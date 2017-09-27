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
            if isEmpty(obj.Mask)
                error('Correct_Image_Background:noMask','Must first compute the object mask, using Compute_Object_Mask(), before calling this function.')
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
            %             varargout{1} = x;
            
            % *Compute DAPI forground and stripe correction*
            [flatteningSurface, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(t,x(:,1),x(:,2),x(:,4),x(:,3));
            
            t.FG_smooth = flatteningSurface;
            t.FG_Xstripe = Xstripe;
            obj.DAPI_G1_Area = G1Area;
            obj.DAPI_G1_Intensity = flatteningSurface;
            obj.DAPI_G1_Idx = G1_idx;
            
            %             counter = 2;
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
                %                 varargout{counter} = struct('intensity',I,'name',name);
                %                 counter = counter + 1;
                
                % - Select all G1 nuclei
                % - Fit intensity surface with small-smoothing-radius
                % surface (to get accurate surface)
                % - Fit any x-stripe
                % - Fit the 2% - 4% intensity range with smooth
                % surface; this will be considered the background due
                % to staining inhomogenaity. (this is treated as a
                % foreground when measuring properties - i.e., the
                % intensities are divided by this surface.)
                
                % *Compute X-stripe correction*
                [FG, Xstripe] = Compute_Channel_Corrections(t, x(G1_idx,1), x(G1_idx,2), I(G1_idx), [], char(name));
                t.FG_smooth = FG;
                t.FG_Xstripe = Xstripe;
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
            N_img = numel(obj.Channel_TiffImgs);
            try
                
                % Set block sizes, get background/foreground functions
                hasBG = false(1,N_img);
                hasFG = false(1,N_img);
                BG_fun = cell(1,N_img);
                FG_fun = cell(1,N_img);
                for i = 1:N_img
                    t = obj.Channel_TiffImgs(i);
                    t.blockSize = obj.FeatureComputation_BlockSize;
                    hasBG(i) = ~isempty(t.BG_smooth);
                    hasFG(i) = ~isempty(t.FG_smooth);
                    
                    if hasBG(i)
                        BG_fun{i} = generateFunction(t, t.BG_smooth, t.BG_Xstripe, true);
                    end
                    if hasFG(i)
                        FG_fun{i} = generateFunction(t, t.FG_smooth, t.FG_Xstripe, true);
                    end
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
                    
                    for i = 1:N_img, obj.Channel_TiffImgs(i).open(); end
                    obj.Mask.open();
                    
                    for blck_y = 1:numBlcks(1)
                        
                        % Get image channels for current block
                        for i = N_img:-1:1
                            [tmpI,x,y] = getBlock(obj.Channel_TiffImgs(i), blck_x, blck_y);
                            if i == N_img
                                I = zeros([size(tmpI),3],'like',tmpI);
                            end
                            if hasBG(i)
                                tmpI = tmpI - BG_fun{i}(x,y);
                            end
                            if hasFG(i)
                                tmpI = tmpI ./ FG_fun{i}(x,y);
                            end
                            I(:,:,i) = tmpI;
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
                    for i = 1:N_img, obj.Channel_TiffImgs(i).close(); end
                    obj.Mask.close();
                    progress.iteration_end(); % Update progress counter
                end
                
                % Remove extra allocated features
                features(object_count:end,:) = [];
            catch ME
                for i = 1:N_img, obj.Channel_TiffImgs(i).close(); end
                obj.Mask.close();
                if obj.Use_GPU
                    gpuDevice([]);
                end
                fprintf('CellExperiment-Compute_Object_Features : entered catch statement\n')
                rethrow(ME)
            end
        end
        
        function objImgs = Extract_Object_Images(obj, bboxes, channels, patchSize, resizeObject)
            % EXTRACT_OBJECT_IMAGES Use bounding boxes to extract regions
            % of interest (objects) from the channel images.
            %
            % objImgs = Extract_Object_Images(obj, bboxes, channels, patchSize, resizeObject)
            %
            % Input
            %   bboxes : N x 4 array giving the bounding boxes of the N
            %     objects to extract
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
            
            % Check inputs
            [validChannel,channelIdx] = ismember(obj.Channel_Names,channels);
            channelIdx = channelIdx(channelIdx>0);

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
            if size(bboxes,2) ~= 4 || numel(size(bboxes)) ~= 2 || ...
                    any(bboxes(:,1)+bboxes(:,3)>t.imageSize(2)) || ...
                    any(bboxes(:,2)+bboxes(:,4)>t.imageSize(1)) || ...
                    any(any(bboxes(:,1:2)<1))
                error('Extract_Object_Images:badInput','bboxes must be a n x 4 array, where each row gives an object bounding box.');
            end
            
            
            % If patch size is specified and not resizing objects, then get
            % the new bounding boxes.
            if ~resizeObject && ~ischar(patchSize)
                % Need to create new bounding boxes with W x H = patchSize
                % centered around the given bounding boxes center
                hlfSize = ceil((patchSize-1)/2);
                center = ceil(bboxes(:,1:2) + bboxes(:,3:4)/2);
                LT = center - hlfSize;
                RB = center + hlfSize;
                
                LT(LT(:,1)<1,1) = 1;
                LT(LT(:,2)<1,2) = 1;
                RB(RB(:,1)>t.imageSize(2),1) = t.imageSize(2);
                RB(RB(:,2)>t.imageSize(1),2) = t.imageSize(1);
                
                bboxes = [LT, RB - LT + 1];
            end
            
            % Get tile indices for each bounding box
            obj_tiles = t.computeTileNumbers(bboxes);
            
            % Get unique list of tiles and the number of bounding boxes
            % each tile is associated with
            N_box = size(bboxes,1);
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
            
            for i = channelIdx, obj.Channel_TiffImgs(i).open(); end % Open all images
            
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
                                    tmp = obj.Channel_TiffImgs(channelIdx(ch_num)).getTile(tN); % Read in image tile
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
                for i = channelIdx, obj.Channel_TiffImgs(i).close(); end
                if obj.Use_GPU
                    gpuDevice([])
                end
                fprintf(2,'CellExperiment-Extract_Object_Images, entered catch statement\n')
                rethrow(ME)
            end
            for i = channelIdx, obj.Channel_TiffImgs(i).close(); end
        end
    end
    
    methods (Access = private)
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
                obj.Channel_TiffImgs(i).BG_smooth = [];
                obj.Channel_TiffImgs(i).BG_Xstripe = [];
                obj.Channel_TiffImgs(i).FG_smooth = [];
                obj.Channel_TiffImgs(i).FG_Xstripe = [];
            end
        end
    end
end
