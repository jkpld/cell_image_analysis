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

        Surface_Smoothing_Radius(1,1) double = NaN % [mm]
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
                    % Create the TiffImg
                    obj.Channel_TiffImgs(i) = TiffImg(channel_paths(i));

                    % Ensure all tileSize's are the same for each image.
                    tileSize = obj.Channel_TiffImgs(i).tileSize;
                    imageSize = obj.Channel_TiffImgs(i).imageSize;
                    if i < numel(channel_paths)
                        if ~isequal(tileSize,prevTileSize)
                            error('CellExperiment:differentTileSizes','The tile size of all channels in a CellExperiment must be the same.')
                        end
                        if ~isequal(imageSize,prevImageSize)
                            error('CellExperiment:differentImageSizes','The image size of all channels in a CellExperiment must be the same.')
                        end
                    end
                    prevTileSize = tileSize;
                    prevImageSize = imageSize;
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
                [FG_f, FG_o, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));

                obj.DAPI_G1_Area = G1Area;
                obj.DAPI_G1_Idx = G1_idx;

                t.FG_offset = FG_o;
                t.FG_factor = FG_f;
                t.FG_stripeX = Xstripe;

                % Set the Area_Normalizer of the nucleiPartitioner
                FG_o_fun = interpolator2d(FG_o.x, FG_o.y, FG_o.Z, false);
                obj.nucleiPartitioner.Intensity_Normalizer = @(I,x,y) I - FG_o_fun(x,y);
                
                if thresholdCorrectedImage
                    % Use the background and foreground corrected image to
                    % create -potentially better- object masks. This takes a
                    % bit longer.

                    minFG = min(t.FG_factor.Z(:));
                    t.FG_factor.Z = t.FG_factor.Z./minFG;

                    % *Compute threshold*
                    t.blockSize = obj.SurfaceComputation_BlockSize;
                    t.Surface_Smoothing_Radius = NaN; % smoothing is not necessary for this surface since the median threshold value will be used as a global threshold.
                    t.Compute_Threshold(1);
                    t.Surface_Smoothing_Radius = obj.Surface_Smoothing_Radius;
                    t.FG_factor.Z = t.FG_factor.Z .* minFG;

                    % Assume already the changes to the objects made by
                    % this thresholding are small. - This might not be
                    % valid
                    %   % *Measure basic props*
                    %   t.blockSize = obj.FeatureComputation_BlockSize;
                    %   x = Measure_BasicProps(t);
                    % 
                    %   x(x(:,3) < minimumObjectSizeForCorrection,:) = [];
                    %   x = double(x); % Need doubles for Compute_DAPI_Corrections
                    % 
                    %   % *Compute G1 intensity and area surfaces for normalization*
                    %   [G1Dapi, ~, ~, G1Area] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));
                    % 
                    %   % Set the Intensity_Normalizer of the nucleiPartitioner
                    %   obj.nucleiPartitioner.Intensity_Normalizer = interpolator2d(G1Dapi.x, G1Dapi.y, G1Dapi.Z);
                end

                % Set the Area_Normalizer of the nucleiPartitioner
                Area_fun = interpolator2d(G1Area.x, G1Area.y, G1Area.Z, false);
                obj.nucleiPartitioner.Area_Normalizer = @(A,x,y) A ./ Area_fun(x,y);
            end
            % *Write object mask*
            t.blockSize = obj.FeatureComputation_BlockSize;
            Write_Object_Mask(t, mask_path, obj.nucleiPartitioner);

            % *Create the TiffImg object for the object mask*
            obj.Mask = TiffImg(mask_path);
        end

        function varargout = Correct_Image_Backgrounds(obj, DEBUG)
            if isEmpty(obj.Mask)
                error('Correct_Image_Background:noMask','Must first compute the object mask, using Compute_Object_Mask(), before calling this function.')
            end

            if nargin < 2
                DEBUG = 0;
            end
            
            % Set Smoothing_Surface_Radius, Use_Parallel, and Verbose;
            % clear any backgrounds and foregrounds.
            setProperties(obj)
            clearCorrections(obj)

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
            
            if DEBUG && nargout > 0
                varargout{1} = x;
            end

            % *Compute DAPI forground and stripe correction*
            [FG_f, FG_o, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));

            obj.DAPI_G1_Area = G1Area;
            obj.DAPI_G1_Idx = G1_idx;
            
            t.FG_offset = FG_o;
            t.FG_factor = FG_f;
            t.FG_stripeX = Xstripe;

            counter = 2;
            for i = 1:numel(obj.Channel_Names)

                % Skip the DAPI image
                name = obj.Channel_Names(i);
                if name == "DAPI"
                    continue;
                end

                % shorter name
                t = obj.Channel_TiffImgs(i);

                % *Compute background*
                t.blockSize = obj.SurfaceComputation_BlockSize;
                t.Compute_Background(1,'Object_Mask',obj.Mask);

                % *Measure basic props*
                t.blockSize = obj.FeatureComputation_BlockSize;
                I = Measure_Intensity(t,obj.Use_Parallel,'Object_Mask',obj.Mask);
                I = double(I);
                
                if DEBUG && nargout > counter-1
                    varargout{counter} = struct('intensity',I,'name',name); %#ok<AGROW>
                    counter = counter + 1;
                end

                % *Compute channel correction*
                [FG, Xstripe] = Compute_Channel_Corrections(t, x(G1_idx,1), x(G1_idx,2), I(G1_idx), [], char(name));
                t.FG_factor = FG;
                t.FG_stripeX = Xstripe;
            end
        end

        function [features, feature_names] = Compute_Object_Features(obj)

            % Check to see if we have a mask
            if isEmpty(obj.Mask)
                error('Compute_Object_Features:noMask','Must first compute the object mask, using Compute_Object_Mask(), before calling this function.');
            end

            % Check to see that FeatureGroups have been added
            if numel(obj.featureExtractor.featureGroups) == 0 || ...
                    all(isempty(obj.featureExtractor.featureGroups))
                error('Compute_Object_Features:noFeatureGroups','No FeatureGroups have been added to the FeatureExtractor.');
            end

            % Set extractor channel names
            obj.featureExtractor.channels = obj.Channel_Names;

            % Number of channels
            N_ch = numel(obj.Channel_TiffImgs);
            try

                % Set block sizes, get correction functions
                CorrectionFunction = cell(1,N_ch);
                postProcIntegratedInt_Offset = cell(1,N_ch);
                
                for i = 1:N_ch
                    t = obj.Channel_TiffImgs(i);
                    t.blockSize = obj.FeatureComputation_BlockSize;
                    
                    CorrectionFunction{i} = generateCorrectionFunction(t);
                    postProcIntegratedInt_Offset{i} = generateFunction(t,t.ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression,false,false);
                end

                % Initialize the feature sizes
                feature_names = [obj.featureExtractor.featureGroups.FeatureNames];
                N_feat = numel(feature_names);
                features = zeros(50000, N_feat, 'single');
                object_count = 1;

                % Number of blocks
                numBlcks = t.numBlcks;

                % Create progress display
                progress = displayProgress(numBlcks(2),'number_of_displays', 15,'active',obj.Verbose, 'name', 'Extracting object features,');
                progress.start();


                % iterate over x blocks
                for blck_x = 1:numBlcks(2)

                    for i = 1:N_ch, obj.Channel_TiffImgs(i).open(); end
                    obj.Mask.open();

                    for blck_y = 1:numBlcks(1)

                        % Get image channels for current block
                        for i = N_ch:-1:1
                            [tmpI,x,y] = getBlock(obj.Channel_TiffImgs(i), blck_x, blck_y);
                            if i == N_ch
                                I = zeros([size(tmpI),3],'like',tmpI);
                            end
                            I(:,:,i) = CorrectionFunction{i}(tmpI,x,y);
                        end

                        % Get mask for current block
                        BW = getBlock(obj.Mask, blck_x, blck_y);

                        % Smooth the image - Assume all channel smoothing
                        % kernels are the same
                        if obj.Use_GPU, I = gpuArray(I); end
                        Is = imfilter(I, obj.Channel_TiffImgs(1).Image_Smooth_Kernel,'symmetric'); clear I;
                        if obj.Use_GPU, I = gather(Is); clear Is; end

                        % Extract object features from curret block
                        x_i = obj.featureExtractor.Compute(BW,I, [t.xEdges(blck_x),t.yEdges(blck_y)]);

                        % Combine with all other object features
                        Ni_obj = size(x_i,1);
                        if Ni_obj + object_count-1 > size(features,1)
                            features = [features; zeros(50000, N_feat,'single')]; %#ok<AGROW>
                        end
                        features(object_count : Ni_obj + object_count - 1,:) = x_i;
                        object_count = object_count + Ni_obj;
                    end

                    % Prevent memory buildup
                    for i = 1:N_ch, obj.Channel_TiffImgs(i).close(); end
                    obj.Mask.close();
                    progress.iteration_end(); % Update progress counter
                end

                % Remove extra allocated features
                features(object_count:end,:) = [];
                
                % Post-process object intensities
                % If an image channel has an
                % ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression,
                % then subract that offset from all integrated intensity
                % features. Additionally, divide this offset by the object
                % area, and subtract this mean intensity correction from
                % all mean intensity features.
                
                % All feature group channels and feature group names
                featChannels = [obj.featureExtractor.featureGroups.Channel];
                featNames = [obj.featureExtractor.featureGroups.GroupName];
                crctn = line('XData',[],'YData',[],'Parent',[]);
                
                % Helpers
                iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}(); % inline if
                curly = @(x, varargin) x{varargin{:}}; % extract cell array element - Use cell array for creating linear sequence of evaluations
                doty = @(x, name) x.(name); % extract object property
                
                channels_withPostProc = find(~cellfun(@isempty, postProcIntegratedInt_Offset));
                for i = channels_withPostProc
                    chGroupIdx = find(obj.Channel_Names(i) == featChannels);
                    
                    % Create helpers that will allow for computing the
                    % intensity correction and mean correction only once.
                    % The intensity correction will be stored as XData in a
                    % line handle, and the mean correction will be stored
                    % as the YData in a line handle. A line handle is used
                    % as a simple way of keeping the data arround.
                    crctn.XData = []; crctn.YData = [];
                    II_offset = @(x) iif(isempty(x.XData), doty(curly({set(x,'XData',postProcIntegratedInt_Offset{i}(obj.Object_Centroid(:,1),obj.Object_Centroid(:,2))),x},2),'XData'), true, x.XData);
                    MuI_offset = @(x) iif(isempty(x.YData), doty(curly({set(x,'YData',II_offset(x)./obj.Object_Area),x},2),'YData'), true, x.YData);
                    
                    % Iterate over the feature groups of the current
                    % channel
                    for grp = chGroupIdx
                        % Apply corrections depending on what type of
                        % feature group it is
                        switch featNames(chGroupIdx)
                            case 'Intensity'
                                II_idx = feature_names == "Intensity_" + obj.Channel_Names(i) + "_IntegratedIntensity";
                                MuI_idx = feature_names == "Intensity_" + obj.Channel_Names(i) + "_MeanIntensity";
                                features(:,II_idx) = features(:,II_idx) - II_offset(crctn);
                                features(:,MuI_idx) = features(:,MuI_idx) - MuI_offset(crctn);
                            case 'RadialIntensity'
                                MuI_idx = startsWith(feature_names, "RadialIntensity_" + obj.Channel_Names(i) + "_Mean");
                                features(:,MuI_idx) = features(:,MuI_idx) - MuI_offset(crctn);
                            case 'Granularity'
                                II_idx = startsWith(feature_names, "Granularity_" + obj.Channel_Names(i));
                                features(:,II_idx) = features(:,II_idx) - II_offset(crctn);
                            otherwise
                                % No other channel needs to be corrected
                        end
                    end
                end
                
                % Post-process Granularity features. We need to divide the
                % granular spectrum 1:L by the 0 element to normalize them.
                % *This needs to happen after the intensity correction,
                % which is why it must be done here.*
                gran_idx = find(featNames == "Granularity");
                toRemove = zeros(numel(gran_idx)); % need to remove one feature for each Granularity group.
                names = {groups.FeatureNames};
                num_xi = cumsum([0,cellfun(@numel, names)]); % group bounds
                
                for i = gran_idx
                    start_idx = num_xi(gran_idx);
                    end_idx = num_xi(gran_idx+1);
                    
                    start_idx = start_idx + 1;
                    toRemove(i) = start_idx;
                    
                    features(:,start_idx+1:end_idx) = features(:,start_idx+1:end_idx) ./ features(:,start_idx);
                end
                
                features(:,toRemove) = [];
                feature_names(toRemove) = [];
                
            catch ME
                for i = 1:N_ch, obj.Channel_TiffImgs(i).close(); end
                obj.Mask.close();
                if obj.Use_GPU
                    gpuDevice([]);
                end
                fprintf('CellExperiment-Compute_Object_Features : entered catch statement\n')
                rethrow(ME)
            end
        end

        function objImgs = Extract_Object_Images(obj, location, channels, patchSize, resizeObject)
            % EXTRACT_OBJECT_IMAGES Use bounding boxes to extract regions
            % of interest (objects) from the channel images.
            %
            % objImgs = Extract_Object_Images(obj, bboxes, channels, patchSize, resizeObject)
            %
            % Input
            %   locations : Either a N x 4 array giving the bounding boxes
            %     of the N objects to extract, or an N x 2 array giving the
            %     centroids of the N objects to extract. 
            %   channels : The names of the channels from which to extract
            %     images.
            %   patchSize : The size of the image to extract (1x2 array),
            %     or the string 'objectSize'. If a 1x2 array [M,L], then
            %     the output is a M x L x K x N image, where K is the
            %     number of channels, and N is the number of objects. If
            %     'objectSize' is selected (and resizeObject = false), then
            %     the ouput is a Nx1 cell array where each element contains
            %     the HxWxK image for each object, where W and H are the
            %     width and hight of the object (from the boundind box).
            %     When the input is a 1x2 array, the objects will be padded
            %     or cropped as necessary. Note, patchSize will be rounded
            %     to odd number.
            %   resizeObject : logical flag. If true, then all objects will
            %     be resized to the the patchSize.
            %
            % Output
            %   objImgs : M x L x K x N array or N x 1 cell array with Hi x
            %     Wi x K array elements
            %
            % Note: Resizing objects is only valid when the bounding boxes
            % are given. Setting patchSize to 'objectSize' is also only
            % valid when the bounding boxes are given.

            % Check inputs
            [validChannel,channelIdx] = ismember([obj.Channel_Names,"Mask"],channels);
            channelIdx = channelIdx(channelIdx>0);
            isMask = channels == "Mask";

            if nargin < 5
                resizeObject = false;
            end
            
            szL = size(location);
            if ~any(validChannel)
                error('Extract_Object_Images:unknownChannel','The channel names given do not match and of the CellExperiment''s Channel_Names.')
            end

            t = obj.Channel_TiffImgs(channelIdx(1));

            if ~(ischar(patchSize) && strcmp(patchSize,'objectSize')) && ~isequal(size(patchSize),[1,2])
                error('Extract_Object_Images:badInput','patchSize must be a [1,2] array giving the size of the extracted images, or the string ''objectSize''.');
            end
            if (resizeObject ~= 0) && (resizeObject ~= 1)
                error('Extract_Object_Images:badInput','resizeObject must be logical.');
            end
            
            if (szL(2) ~= 2 && szL(2) ~= 4) || numel(szL) ~= 2
                error('Extract_Object_Images:badInput','The location input must either be an Nx2 array giving the object centroids or an Nx4 array giving the object bounding boxes.')
            end
            
            if szL(2) == 2
                
                % Round centroids to nearest integer
                location = round(location);
                
                if strcmp(patchSize,'objectSize')
                    error('Extract_Object_Images:badInput','Specifying ''objectSize'' requires that the location input are the bounding boxes (n x 4).')
                end

                if resizeObject
                    error('Extract_Object_Images:badInput','Resizing the objects is only possible when the bounding boxes are given, not the object centroids.')
                end

                if any(any(location<1)) && any(any(location>flip(t.imageSize)))
                    error('Extract_Object_Images:badInput','The object centroids provided are not valid because they go outside of the image bounds.');
                end
            end
            
            if szL(2) == 4 && (any(any(location(:,1:2)<1)) || ...
                    any(any(location(:,1:2)+location(:,3:4)>flip(t.imageSize))))
                error('Extract_Object_Images:badInput','The bounding boxes provided are not valid because they go outside of the image bounds.');
            end

            % If patch size is specified and not resizing objects, then get
            % the new bounding boxes.
            if (~resizeObject && ~ischar(patchSize)) || szL(2)==2
                % Need to create new bounding boxes with W x H = patchSize
                % centered around the given bounding boxes center
                hlfSize = ceil((patchSize-1)/2);
                if szL(2)==2
                    center = round(location);
                else
                    center = ceil(location(:,1:2) + location(:,3:4)/2);
                end
                LT = center - hlfSize;
                RB = center + hlfSize;

                LT(LT(:,1)<1,1) = 1;
                LT(LT(:,2)<1,2) = 1;
                RB(RB(:,1)>t.imageSize(2),1) = t.imageSize(2);
                RB(RB(:,2)>t.imageSize(1),2) = t.imageSize(1);

                bboxes = [LT, RB - LT + 1];
            else
                bboxes = locations;
            end

            % Get tile indices for each bounding box
            obj_tiles = t.computeTileNumbers(bboxes);

            % Get unique list of tiles and the number of bounding boxes
            % each tile is associated with
            N_box = szL(1);
            all_tiles = cell(N_box,1);
            for i = 1:N_box
                [X,Y] = meshgrid(obj_tiles(i,1):obj_tiles(i,3),obj_tiles(i,2):obj_tiles(i,4));

                % compute linear indices - rmember that the obj_tile
                % indices are 0 based
                all_tiles{i} = Y*t.numTiles(2) + X;
            end

            tiles = cellfun(@(x) x(:), all_tiles,'UniformOutput',false);
            tiles = cat(1,tiles{:}); % tiles = m x 1
            [tiles,~,number] = unique(tiles); % tiles = n x 1, number = m x 1 in [1,n]
            number = accumarray(number,1); % number = n x 1

            % Create an array to convert from the tile number to the index
            % of the tile stored in a tile cache
            tile_to_idx = sparse(double(tiles),ones(numel(tiles),1),1:numel(tiles),double(max(tiles)),1,numel(tiles));

            % Get each object image --------------------------------------
            N_ch = numel(channelIdx); % Number of channels
            I = cell(numel(tiles),1); % Create tile cache

            % Create output array
            if ischar(patchSize)
                objImgs = cell(N_box,1);
            else
                objImgs = zeros([patchSize,N_ch,N_box],t.workingClass);
            end

            % Create an array of the TiffImgs
            for i = numel(channelIdx):-1:1
                if isMask(channelIdx(i))
                    tiffImgs(i) = obj.Mask;
                else
                    tiffImgs(i) = obj.Channel_TiffImgs(channelIdx(i));
                end
            end
            
            for i = 1:N_ch, tiffImgs(i).open(); end % Open all images

            tmp_y_tileIdx = 1:t.tileSize(1); % create helper array
            tmp_x_tileIdx = 1:t.tileSize(2); % create helper array

            try
                for i = 1:N_box % iterate over each bounding box

                    % Construct image patch for bounding box -------------
                    tileNums = all_tiles{i}; % Normal objects should have at most 4 tiles
                    tmpI = zeros([size(tileNums).*t.tileSize,N_ch],t.workingClass);

                    for tN_x = 1:size(tileNums,2)
                        for tN_y = 1:size(tileNums,1)
                            tN = tileNums(tN_y,tN_x);

                            % Get the index of this tile in the image tile
                            % cache
                            idx = full(tile_to_idx(tN,1));

                            if isempty(I{idx})
                                % If tile has not been read in before, read
                                % in now. Initialize tile image
                                tileImg = zeros([t.tileSize,N_ch],t.workingClass);
                                for ch_num = 1:N_ch
                                    tmp = tiffImgs(ch_num).getTile(tN); % Read in image tile
                                    tileImg(1:size(tmp,1),1:size(tmp,2),ch_num) = tmp; % Add to tileImg
                                end

                                % Save tile for later if it is needed for
                                % another object.
                                if number(idx) > 1
                                    I{idx} = tileImg;
                                    number(idx) = number(idx) - 1;
                                end
                            else
                                % If tile was read in, then simply use that
                                tileImg = I{idx};
                                number(idx) = number(idx) - 1;

                                % If tile is not needed in future, then
                                % clear it from memory
                                if number(idx) == 0
                                    I{idx} = [];
                                end
                            end

                            % Add the tile to the image patch we need
                            tmpI(tmp_y_tileIdx + (tN_y-1)*t.tileSize(1), tmp_x_tileIdx + (tN_x-1)*t.tileSize(2), :) = tileImg;
                        end
                    end

                    % Extract ROI ----------------------------------------

                    % Offset left-top point of bounding box
                    bbox = bboxes(i,:);
                    bbox(:,1:2) = bbox(:,1:2) - obj_tiles(i,1:2).*t.tileSize;

                    roi = tmpI(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1,:);

                    if ischar(patchSize)
                        objImgs{i} = roi;
                    else
                        if resizeObject
                            roi = imresize(roi,patchSize);
                            objImgs(:,:,:,i) = roi;
                        else
                            % If the object is on the edge of the image,
                            % then the roi might not be the patchSize. In
                            % this case, we must center the roi in a box
                            % with size patchSize.
                            padSize = patchSize - [size(roi,1),size(roi,2)];
                            offset = floor(padSize/2);
                            objImgs(1+offset(1):size(roi,1)+offset(1), 1+offset(2):size(roi,2)+offset(2), :, i) = roi;
                        end
                    end
                end
            catch ME
                for i = 1:N_ch, tiffImgs(i).close(); end % Open all images
                if obj.Use_GPU
                    gpuDevice([])
                end
                fprintf(2,'CellExperiment-Extract_Object_Images, entered catch statement\n')
                rethrow(ME)
            end
            for i = 1:N_ch, tiffImgs(i).close(); end % Open all images
        end
    end

    methods %(Access = private)
        function setProperties(obj)
            for i = 1:numel(obj.Channel_TiffImgs)
                % Set some properties
                obj.Channel_TiffImgs(i).Surface_Smoothing_Radius = obj.Surface_Smoothing_Radius;
                obj.Channel_TiffImgs(i).Use_GPU = obj.Use_GPU;
                obj.Channel_TiffImgs(i).Verbose = obj.Verbose;
            end
        end
        function clearCorrections(obj)
            for i = 1:numel(obj.Channel_TiffImgs)
                % clear backgrounds and foregrounds
                obj.Channel_TiffImgs(i).BG_offset = [];
                obj.Channel_TiffImgs(i).BG_stripeX = [];
                obj.Channel_TiffImgs(i).FG_offset = [];
                obj.Channel_TiffImgs(i).FG_factor = [];
                obj.Channel_TiffImgs(i).FG_stripeX = [];
            end
        end
    end
end
