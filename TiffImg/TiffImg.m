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

        Background = []
        Background_smooth = []
        Background_stripe = []
        Background_mean = []

        Foreground = []
        Foreground_smooth = []
        Foreground_stripe = []
        Foreground_mean = []

        Sharpness = []
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
    end
    properties (Access = private, Hidden)
        threshold_fun = []
        BG_fun = []
        BG_smooth_fun = []
        BG_stripe_fun = []
        BG_mean = []

        FG_fun = []
        FG_smooth_fun = []
        FG_stripe_fun = []
        FG_mean = []

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
                tmp = tifflib('readEncodedTile',obj.FileID,cTiles(t));%readEncodedTile(t,cTiles(i));
                if t == 1
                    I = zeros(obj.blockSize,obj.blockSize,'like',tmp);
                end

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

            if nargout > 1
                y = (1:(maxY-1)) + (blck_y-1)*obj.blockSize;
                x = (1:(maxX-1)) + (blck_x-1)*obj.blockSize;
            end
        end

        % Functions defined externally
        Z = smoothSurface(obj,z,type);
        th = Compute_Threshold(obj, removeBackgroundFirst);
        bg = Compute_Background(obj);
        fg = Compute_Foreground(obj);
        sh = Compute_Sharpness(obj);

    end

    methods (Access = private)
        
        function fun = generateFunction(obj, z, x_stripe, const)
            % Create a function handel that can be used to evaluate
            % f(x,y) =  z(x,y) + z_stripe(x) + const
            z = double(z); % The fast interpolation function used requires double
            
            if (nargin < 4) || isempty(const)
                const = 0;
            end
            
            const = double(const); % Ensure double
            
            if (nargin < 3) || isempty(x_stripe)
                fun = @(x,y) cast(const + ...
                    interp2mex(z, (x - obj.xCenters(1))/obj.blockSize, (y - obj.yCenters(1))/obj.blockSize), ...
                    obj.imageClass); % cast final result to same class as image
            else
                x_stripe = double(x_stripe); % Ensure double
                
                fun = @(x,y) cast(const + x_stripe(x) + ...
                    interp2mex(z, (x - obj.xCenters(1))/obj.blockSize, (y - obj.yCenters(1))/obj.blockSize), ...
                    obj.imageClass); % cast final result to same class as image
            end
            
            % This only performs nearest neighbor extrapolation outside of
            % the defined limits. Therefore, the half blockSize on the
            % image boarder could have jagged edges.
            %
            % This will only be approximately correct on the right and
            % bottom edges of the image because the size of the right and
            % bottom pixels could be smaller than the rest of them -- due
            % to blocks that are not full.

        end
    end
end % class
