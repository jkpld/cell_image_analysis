classdef CellExperiment < handle
    properties
        Experiment_Description(1,1) string
    end
    properties (SetAccess = private)
        Channel_Names(1,:) string
        Channel_TiffImgs(1,:) TiffImg
        Mask(1,1) TiffImg

        Object_Centroid(:,2) single
        Object_Area(:,1) single
        
        DAPI_G1_Area(1,1) struct
        DAPI_G1_Intensity(1,1) struct
        DAPI_G1_Idx(:,1) double 
    end
    properties        
        nucleiPartitioner(1,1) ObjectPartitioner = None
        featureExtractor(1,1) FeatureExtractor

        Surface_Smoothing_Radius(1,1) double = 2 % [mm]
        SurfaceComputation_BlockSize(1,1) double = 1024 % [pixels]
        FeatureComputation_BlockSize(1,1) double = 4096 % [pixels]

        Use_GPU(1,1) logical = true
        Use_Parallel(1,1) logical = true
        Verbose(1,1) logical = true
    end
    properties (SetAccess = private)%, Hidden)
        Threshold_Corrected_Image(1,1) logical
        TiffImg_for_Generating_Mask(1,1) TiffImg
    end

    methods
        function obj = CellExperiment(channel_names, channel_paths)
            if nargin > 0
                if numel(channel_names) ~= numel(channel_paths)
                    error('CellExperiment:badInput','The number of channel_names must match the number of channel_paths.')
                end

                obj.Channel_Names = channel_names;

                if ~contains(obj.Channel_Names, "DAPI")
                    error('CellExperiment:missingDAPI','Missing DAPI channel. A DAPI channel is required for this class.')
                end

                for i = numel(channel_paths):-1:1
                    obj.Channel_TiffImgs(i) = TiffImg(channel_paths(i));
                end
            end
        end

        function Create_Object_Mask(obj, mask_path, thresholdCorrectedImage, minimumObjectSizeForCorrection)

            % Input checking
            narginchk(2,4)

            if nargin < 3
                thresholdCorrectedImage = false;
            end
            if nargin < 4
                minimumObjectSizeForCorrection = 150;
            end

            % Save the thresholdCorrectedImage flag for later use
            obj.Threshold_Corrected_Image = thresholdCorrectedImage;

            % Grab the DAPI TiffImg
            t_ch = obj.Channel_TiffImgs(obj.Channel_Names == "DAPI");
            
            % Make a copy of the image before going farther
            t = copy(t_ch);
            obj.TiffImg_for_Generating_Mask = t; % Save for later reference

            % Set some properties
            t.Surface_Smoothing_Radius = obj.Surface_Smoothing_Radius;
            t.Use_GPU = obj.Use_GPU;
            t.Verbose = obj.Verbose;

            % *Compute threshold*
            t.blockSize = obj.SurfaceComputation_BlockSize;
            t.Compute_Threshold();

            % The nuclei partition could need an Area_Normalizer or
            % Intensity_Normalizer. If so, then we need to compute
            % information about each object to create these normalizers
            % before partitioning. We also need to compute this information
            % if we are using the corrected image to threshold.
            if thresholdCorrectedImage || ...
                    (~isempty(obj.nucleiPartitioner.Area_Normalizer) || ...
                    ~isempty(obj.nucleiPartitioner.Intensity_Normalizer)) || ...
                    obj.nucleiPartitioner.Restricted_Partitioning

                % *Compute background*
                t.Compute_Background(1);

                % *Measure basic props*
                t.blockSize = obj.FeatureComputation_BlockSize;
                x = Measure_BasicProps(t);

                x(x(:,3) < minimumObjectSizeForCorrection,:) = [];
                x = double(x);
                
                % *Compute DAPI forground and stripe correction*
                [flatteningSurface, Xstripe, G1Area] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));

                t.FG_smooth = flatteningSurface;
                t.FG_Xstripe = Xstripe;

                if thresholdCorrectedImage
                    % Use the background and foreground corrected image to
                    % create -potentially better- object masks. This takes a
                    % bit longer.

                    t.FG_smooth.Z = t.FG_smooth.Z./min(t.FG_smooth.Z(:));

                    % *Compute threshold*
                    t.blockSize = obj.SurfaceComputation_BlockSize;
                    t.Surface_Smoothing_Radius = NaN; % smoothing is not necessary for this surface since the median threshold value will be used as a global threshold.
                    t.Compute_Threshold(1,1);
                    t.Surface_Smoothing_Radius = obj.Surface_Smoothing_Radius;

                    % *Measure basic props*
                    t.blockSize = obj.FeatureComputation_BlockSize;
                    x = Measure_BasicProps(t);

                    x(x(:,3) < minimumObjectSizeForCorrection,:) = [];
                    x = double(x); % Need doubles for Compute_DAPI_Corrections

                    % *Compute G1 intensity and area surfaces for normalization*
                    [G1Dapi, ~, G1Area] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));

                    % Set the Intensity_Normalizer of the nucleiPartitioner
                    obj.nucleiPartitioner.Intensity_Normalizer = interpolator2d(G1Dapi.x, G1Dapi.y, G1Dapi.Z);
                end

                % Set the Area_Normalizer of the nucleiPartitioner
                obj.nucleiPartitioner.Area_Normalizer = interpolator2d(G1Area.x, G1Area.y, G1Area.Z);

            end
            % *Write object mask*
            t.blockSize = obj.FeatureComputation_BlockSize;
            Write_Object_Mask(t, mask_path, obj.nucleiPartitioner);

            % *Create the TiffImg object for the object mask*
            obj.Mask = TiffImg(mask_path);
        end

        function Correct_Image_Backgrounds(obj)
            if isempty(obj.Mask)
                error('Correct_Image_Background:noMask','Must first compute the object mask, Compute_Object_Mask, before calling this function.')
            end

            % Correct DAPI first to get the G1 indices
            % shorter name
            t = obj.Channel_TiffImgs(obj.Channel_Names == "DAPI");

            % *Compute background*
            t.blockSize = obj.SurfaceComputation_BlockSize;
            t.Compute_Background(1,'Object_Mask',obj.Mask);

            % *Measure basic props*
            t.blockSize = obj.FeatureComputation_BlockSize;
            x = Measure_BasicProps(t,'DAPI','Object_Mask',obj.Mask);
            obj.Object_Centroid = x(:,1:2);
            obj.Object_Area = x(:,3);
            x = double(x);
            
            % *Compute DAPI forground and stripe correction*
            [flatteningSurface, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));

            t.FG_smooth = flatteningSurface;
            t.FG_Xstripe = Xstripe;
            obj.DAPI_G1_Area = G1Area;
            obj.DAPI_G1_Intensity = flatteningSurface;
            obj.DAPI_G1_Idx = G1_idx;                
                
            
            for i = 1:numel(obj.Channel_Names)
                
                % Skip the DAPI image
                if obj.Channel_Names(i) == "DAPI"
                    continue;
                end
                
                % shorter name
                t = obj.Channel_TiffImgs(i);

                % *Compute background*
                t.blockSize = obj.SurfaceComputation_BlockSize;
                t.Compute_Background(1,'Object_Mask',obj.Mask);

                % *Measure basic props*
                t.blockSize = obj.FeatureComputation_BlockSize;
                I = Measure_Intensity(t,'Object_Mask',obj.Mask);
                I = double(I);
                
                % - Select all G1 nuclei
                % - Fit intensity surface with small-smoothing-radius
                % surface (to get accurate surface)
                % - Fit any x-stripe
                % - Fit the 2% - 4% intensity range with smooth
                % surface; this will be considered the background due
                % to staining inhomogenaity.

                % *Compute X-stripe correction*
                [FG, Xstripe] = Compute_Channel_Corrections(t, x(G1_idx,1), x(G1_idx,2), I(G1_idx));
                t.FG_smooth = FG;
                t.FG_Xstripe = Xstripe;
            end
        end
    end
end
