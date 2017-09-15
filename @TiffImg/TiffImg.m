classdef TiffImg < handle
    properties
        blockSize(1,1) double
        Use_GPU(1,1) logical = false;
        Verbose(1,1) logical = false;

        % Iamge_Smooth_Kernel - The Kernel used to smooth the image before
        % all operations (directly after reading the image in.)
        Image_Smooth_Kernel = fspecial('gaussian',7,1);

        % Surface_Smoothing_Radius - The radius (in mm) to use when
        % smoothing the surfaces computed (threshold, background,
        % foreground). If set to NaN, then smoothing will be disabled.
        Surface_Smoothing_Radius(1,1) double = NaN;
    end
    properties (SetAccess = private)
        threshold = []

        BG_smooth = []
        BG_Xstripe = []

        FG_smooth = []
        FG_Xstripe = []

        Sharpness = []
    end
    properties (Dependent)
        BackgroundEvalFun(1,1) function_handle
        ForegroundEvalFun(1,1) function_handle
    end
    properties (SetAccess = private)
        FileName
        FileID
        imageSize
        imageClass
        tileSize
        numTiles
        maxSampleValue
        mmPerPixel
    end
    properties (SetAccess = private, Hidden)
        tiles
        tilesPerBlck
        numBlcks
        xEdges
        yEdges
        xCenters
        yCenters
        tile_y_inds
        tile_x_inds
        blck_y_inds
        blck_x_inds
        
        workingClass
    end
    properties (Access = private, Hidden)
        threshold_fun = []
        Threshold_After_Background(1,1) logical = false;
    end

    methods
        function obj = TiffImg(filename)

            if ( nargin==0 )
                obj.FileName = '';
                obj.FileID = uint64(0);
                return;
            end

            % Get the full file path -------------------------------------
            % ------ ( Taken from matlab Tiff constructor )
            fid = fopen(filename,'r','ieee-le');
            % Ensure file exists.
            if fid == -1
                error(message('MATLAB:imagesci:validate:fileOpen',filename));
            end
            % Get the full path name.
            filename = fopen(fid);
            fclose(fid);
            % -----------------------------------------------------------

            obj.FileID = tifflib('open',filename,'r');
            obj.FileName = filename;

            imageDescription = tifflib('getField',obj.FileID,270);
            aqstnInfo = textscan(imageDescription,'%s','Delimiter','|');
            obj.mmPerPixel = sscanf(aqstnInfo{1}{~cellfun('isempty',strfind(aqstnInfo{1},'MPP'))},'MPP = %f')/1000;

            if ~isfinite(obj.mmPerPixel)
                tifflib('close',obj.FileID);
                obj.FileName = '';
                obj.FileID = uint64(0);
                error('tiffReadHelper:MissingMMPerPixel','The image description does not have a MPP (micron per pixel) field, which is required for this class.')
            end

            sampleFormat = tifflib('getField',obj.FileID,339);
            switch sampleFormat
                case 1
                    numBits = tifflib('getField',obj.FileID,258);
                    obj.maxSampleValue = 2^numBits - 1;
                case 2
                    numBits = tifflib('getField',obj.FileID,258);
                    obj.maxSampleValue = 2^(numBits-1) - 1;
                case 3
                    obj.maxSampleValue = 1;
                otherwise
                    tifflib('close',obj.FileID);
                    obj.FileName = '';
                    obj.FileID = uint64(0);
                    error('tiffReadHelper:unsupporetedFormat','The image has an unsupported format, Supported formats are uint, int, and floating point.')
            end

            tmp = tifflib('readEncodedTile',obj.FileID,1);
            obj.imageClass = class(tmp);
            
            % Set the working class. Always work with a float type to have
            % +- numbers. If image is integers less than 16 bit, (2
            % bytes), then use single precision (4 bytes). If image is
            % integers with more than 16 bits then use double.
            % If the image is logical, then just use logical.
            obj.workingClass = obj.imageClass;
            
            if isinteger(tmp)
                if str2double(obj.imageClass(end-1:end)) > 16
                    obj.workingClass = 'double';
                else
                    obj.workingClass = 'single';
                end
            end

            obj.imageSize = [tifflib('getField',obj.FileID,257), tifflib('getField',obj.FileID,256)];
            obj.tileSize = [tifflib('getField',obj.FileID,323), tifflib('getField',obj.FileID,322)];

            numTiles = flip(ceil(obj.imageSize./obj.tileSize));
            obj.tiles = reshape(0:prod(numTiles)-1,numTiles)';
            numTiles = flip(numTiles);
            obj.numTiles = numTiles;

            obj.blockSize = max(obj.tileSize);

            close(obj); % close the image file.
        end

        function set.blockSize(obj,blockSize)
            if ~isequal(blockSize./obj.tileSize,round(blockSize./obj.tileSize)) %#ok<*MCSUP>
                blockSize = unique(round(blockSize./obj.tileSize).*obj.tileSize);
                warning('blockSize adjusted to be multiple of image tile size. New block size is %d.', blockSize)
            end

            obj.blockSize = blockSize;
            obj.tilesPerBlck = blockSize./obj.tileSize;
            obj.numBlcks = ceil(obj.imageSize/blockSize);

            obj.xEdges = [0:floor(obj.imageSize(2)/blockSize), obj.imageSize(2)/blockSize]*blockSize;
            obj.yEdges = [0:floor(obj.imageSize(1)/blockSize), obj.imageSize(1)/blockSize]*blockSize;

            obj.xCenters = obj.xEdges(1:end-1) + diff(obj.xEdges)/2;
            obj.yCenters = obj.yEdges(1:end-1) + diff(obj.yEdges)/2;

            obj.tile_y_inds = 1:obj.tileSize(1);
            obj.tile_x_inds = 1:obj.tileSize(2);
            obj.blck_y_inds = 1:obj.tilesPerBlck(1);
            obj.blck_x_inds = 1:obj.tilesPerBlck(2);
        end

        function delete(obj)
            obj.close()
        end

        function open(obj)
            if ( obj.FileID == 0 )
                obj.FileID = tifflib('open',obj.FileName,'r');
            end
        end

        function close(obj)
            if ( obj.FileID ~= 0 )
                tifflib('close',obj.FileID)
                obj.FileID = uint64(0);
            end
        end

        function [I,x] = getColumn(obj, col_idx)
            
            % initialize array for image stripe
            I = zeros(tiffImg.imageSize(1),tiffImg.tileSize(2),tiffImg.imageClass);
        
            % Read in image stripe
            counter = 1;
            obj.open()
            for tile_y = obj.tiles(:,col_idx)'
                tmp = tifflib('readEncodedTile',obj.FileID,tile_y);%readEncodedTile(t,cTiles(i));
                [m,n] = size(tmp);
                I(counter:counter+m-1, 1:n) = tmp;
                counter = counter + m;
            end
            obj.close()
            I(:,n+1:end) = [];
            
            % Convert image to working class.
            if ~strcmp(obj.imageClass,'logical')
                I = cast(I,obj.workingClass);
                I = I / cast(obj.maxSampleValue,obj.workingClass);
            end
            
            if nargout > 1
                % Compute x pixels
                x = (1:n) + (col_idx-1)*obj.tileSize(2);
            end
        end
        
        function [I,x,y] = getBlock(obj,blck_x,blck_y)
            tmpB_x_inds = obj.blck_x_inds + (blck_x-1)*obj.tilesPerBlck(2);
            tmpB_y_inds = obj.blck_y_inds + (blck_y-1)*obj.tilesPerBlck(1);

            if blck_x == obj.numBlcks(2)
                tmpB_x_inds(tmpB_x_inds > obj.numTiles(2)) = [];
            end
            if blck_y == obj.numBlcks(1)
                tmpB_y_inds(tmpB_y_inds > obj.numTiles(1)) = [];
            end

            cTiles = obj.tiles(tmpB_y_inds,tmpB_x_inds);
            szcTiles = size(cTiles);

            maxX = 0;
            maxY = 0;

            for t = 1:numel(cTiles)

                [j,k] = ind2sub(szcTiles,t);
                
                I = zeros(obj.blockSize,obj.blockSize,'like',obj.imageClass);

                tmpT_y_inds = obj.tile_y_inds + (j-1)*obj.tileSize(1);
                tmpT_x_inds = obj.tile_x_inds + (k-1)*obj.tileSize(2);

                if size(tmp,2) < obj.tileSize(2)
                    tmpT_x_inds(size(tmp,2)+1:obj.tileSize(2)) = [];
                    maxX = tmpT_x_inds(end)+1;
                end
                if size(tmp,1) < obj.tileSize(1)
                    tmpT_y_inds(size(tmp,1)+1:obj.tileSize(1)) = [];
                    maxY = tmpT_y_inds(end)+1;
                end

                I(tmpT_y_inds, tmpT_x_inds) = tmp;
                clear tmp
            end

            if maxX
                I(:,maxX:obj.blockSize) = [];
            else
                maxX = obj.blockSize+1;
            end
            if maxY
                I(maxY:obj.blockSize,:) = [];
            else
                maxY = obj.blockSize+1;
            end

            % Convert image to working class.
            if ~strcmp(obj.imageClass,'logical')
                I = cast(I,obj.workingClass);
                I = I / cast(obj.maxSampleValue,obj.workingClass);
            end
            
            if nargout > 1
                y = (1:(maxY-1)) + (blck_y-1)*obj.blockSize;
                x = (1:(maxX-1)) + (blck_x-1)*obj.blockSize;
            end
        end

        % Functions defined externally
        Z = smoothSurface(obj,z,type);
        th = Compute_Threshold(obj, removeBackgroundFirst);
        bg = Compute_Background(obj, computeXStripeArtifact);
        fg = Compute_Foreground(obj, computeXStripeArtifact);
        st = Compute_StripeArtifact(obj, BG_or_FG)
        sh = Compute_Sharpness(obj);

        function fun = get.BackgroundEvalFun(obj)
            fun = generateFunction(obj, obj.Background_smooth, obj.Background_Xstripe, false);
        end

        function fun = get.ForegroundEvalFun(obj)
            fun = generateFunction(obj, obj.Foreground_smooth, obj.Foreground_Xstripe, false);
        end
    end

    methods (Access = private)

        function fun = generateFunction(obj, z_smooth, x_stripe, meanOffset)
            % Create a function handle that can be used to evaluate
            % f(x,y) =  z_smooth(x,y) + x_stripe(x)
            %
            % If meanOffset is true, then the mean of z(x,y) will be
            % subracted away.

            
            outClass = obj.workingClass;
            

            if (nargin < 3)
                x_stripe = [];
            end

            if isempty(z_smooth)
                % Just have to evaluate the stripe
                if isempty(x_stripe)
                    warning('generateFunction:noInput','There is nothing to create a function from. First compute a surface (threshold, background, foreground, ...)');
                    fun = @(~,~) 0;
                    return;
                end

                fun = @(x,~) x_stripe(x);
                return;
            else

                if (nargin < 4)
                    meanOffset = false;
                end

                if isstruct(z_smooth)
                    x0 = z_smooth.x(1);
                    y0 = z_smooth.y(1);
                    d = z_smooth.x(2) - z_smooth.x(1); % block size
                    z_smooth = z_smooth.Z;
                else
                    x0 = obj.xCenters(1);
                    y0 = obj.yCenters(1);
                    d = obj.blockSize;
                end

                z_smooth = double(z_smooth); % The fast interpolation function used requires double

                if meanOffset
                    z_smooth = z_smooth - mean(z_smooth(:));
                end

                if isempty(x_stripe)
                    fun = @(x,y) cast(...
                        interp2mex(z_smooth, (x - x0)/d, (y - y0)/d), ...
                        outClass); % cast final result to same class as image
                else
                    x_stripe = double(x_stripe); % Ensure double

                    fun = @(x,y) cast(...
                        x_stripe(x) + interp2mex(z_smooth, (x - x0)/d, (y - y0)/d), ...
                        outClass); % cast final result to same class as image
                end
            end

            % Note: interp2mex only performs nearest neighbor extrapolation
            % outside of the defined limits. Therefore, the half blockSize
            % on the image boarder could have jagged edges.
            %
            % This will only be approximately correct on the right and
            % bottom edges of the image because the size of the right and
            % bottom pixels could be smaller than the rest of them -- due
            % to blocks that are not full.

        end
    end
    
    methods (Static)
        th = otsuthresh_scale(I,scale);
    end
end % class
