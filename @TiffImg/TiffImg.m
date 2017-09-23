classdef TiffImg < matlab.mixin.Copyable
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
    properties %(SetAccess = private)
        threshold = []

        BG_smooth = []
        BG_Xstripe = []

        FG_smooth = []
        FG_Xstripe = []

        Sharpness = []
    end

    properties (SetAccess = private)
        FileName
        FileID
        
        mmPerPixel
        
        workingClass
        imageClass
        imageSize
        
        tileSize
        tilesPerBlck
        numBlcks
        numTiles
    end
    properties (SetAccess = private, Hidden)
        maxSampleValue
        tiles
        xEdges
        yEdges
        xCenters
        yCenters
        tile_y_inds
        tile_x_inds
        blck_y_inds
        blck_x_inds
    end
    properties (SetAccess = private)
        threshold_fun = []
        Threshold_After_Background(1,1) logical = false;
        Threshold_After_Foreground(1,1) logical = false;
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

            tmp = tifflib('readEncodedTile',obj.FileID,0);
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

        function tf = isEmpty(obj)
            tf = isempty(obj.FileName);
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
            I = zeros(obj.imageSize(1),obj.tileSize(2),obj.imageClass);
        
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
        
        function [I,y] = getRow(obj, row_idx)
            
            % initialize array for image stripe
            I = zeros(obj.tileSize(1),obj.imageSize(2),obj.imageClass);
        
            % Read in image stripe
            counter = 1;
            obj.open()
            for tile_y = obj.tiles(row_idx,:)
                tmp = tifflib('readEncodedTile',obj.FileID,tile_y);%readEncodedTile(t,cTiles(i));
                [m,n] = size(tmp);
                I(1:m, counter:counter+n-1) = tmp;
                counter = counter + n;
            end
            obj.close()
            I(m+1:end,:) = [];
            
            % Convert image to working class.
            if ~strcmp(obj.imageClass,'logical')
                I = cast(I,obj.workingClass);
                I = I / cast(obj.maxSampleValue,obj.workingClass);
            end
            
            if nargout > 1
                % Compute x pixels
                y = (1:m).' + (row_idx-1)*obj.tileSize(1);
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

            I = zeros(obj.blockSize,obj.blockSize,obj.imageClass);
            
            for t = 1:numel(cTiles)

                [j,k] = ind2sub(szcTiles,t);
                tmp = tifflib('readEncodedTile',obj.FileID,cTiles(t));
    
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
                y = (1:(maxY-1)).' + (blck_y-1)*obj.blockSize;
                x = (1:(maxX-1)) + (blck_x-1)*obj.blockSize;
            end
        end
        
        % Functions defined externally
        Z = smoothSurf(obj,z,type);
        th = Compute_Threshold(obj, removeBackgroundFirst, removeForegroundFirst);
        bg = Compute_Background(obj, computeXStripeArtifact, varargin);
        fg = Compute_Foreground(obj, computeXStripeArtifact, varargin);
        st = Compute_StripeArtifact(obj, BG_or_FG, varargin)
        sh = Compute_Sharpness(obj);
        [x, x_names] = Measure_BasicProps(obj,channelName, varargin);
        feature = Measure_Intensity(tiffImg, Use_Parallel, varargin)
        [flatteningSurface, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(tiffImg,x,y,dapi,area,options)
        
        function fig = plot(obj,name)
            % PLOT Plot one of the computed quantities.
            %
            % fig = tiffImg.plot(name)
            %
            % Input 
            %   name : one of the following strings
            %     'threshold'
            %     'background'
            %     'foreground'
            %     'background_stripe'
            %     'foreground_stripe'
            
            valid = {'threshold','background','foreground','stripe_background','stripe_foreground'};
            internal_name = {'threshold','BG_smooth','FG_smooth','BG_Xstripe','FG_Xstripe'};
            idx = find(strncmpi(name, valid, numel(name)));

            if isempty(idx) || (numel(idx) > 1)
                error('plot:unknownName','The "name" must be an unambigious match to one of the following: \n%s.', string(valid).join(', '))
            else
                name = internal_name{idx};
                name_user = valid{idx};
            end
            
            if isempty(obj.(name))
                error('plot:undefinedSurface','The %s has not been computed yet. Please compute it before trying to plot it.', name_user);
            end
            
            fg = figure('visible','off');
            try
                if idx > 3
                    % Plot x stripe
                    x = (1:obj.imageSize(2))*obj.mmPerPixel;
                    y = obj.(name);
                    plot(x,y,'Color','y','LineWidth',1)
                    title(strrep(name_user,'_',' '))
                    xlabel('x / mm')
                    ylabel('intensity')
                    axis tight
                else
                    % Plot surface
                    x = obj.(name).x;
                    y = obj.(name).y;
                    Z = obj.(name).Z;

                    if any(size(Z)==1) % need scattered surface
                        tri = delaunay(x,y);
                        trisurf(tri,x*obj.mmPerPixel,y*obj.mmPerPixel,Z);
                    else
                        surface(x*obj.mmPerPixel, y*obj.mmPerPixel, Z, 'EdgeColor', 'none');
                    end
                    title(name_user)
                    xlabel('x / mm')
                    ylabel('y / mm')
                    zlabel('intensity')
                    colorbar;
                    axis tight
                    try
                        daspect([1,1,range(Z(:))/range(x*obj.mmPerPixel)])
                    catch
                    end
                end
                setTheme(fg,'dark')
            catch ME
                close(fg)
                rethrow(ME)
            end
            fg.Visible = 'on';
            if nargout > 0
                fig = fg;
            end
        end

        function fun = generateFunction(obj, z_smooth, x_stripe, multiplyStripe)
            % Create a function handle that can be used to evaluate
            % f(x,y) =  z_smooth(x,y) + x_stripe(x)
            %
            % If multiplyStripe is true, then the the stripe will be
            % multiplied by the smoothe surface
            % f(x,y) = z_smooth(x,y) .* x_stripe(x)
            %
            % multiplyStripe is false be default
            %
            % This function will take in vectors and output a matrix!!
            % x = [1 x m]
            % y = [n x 1]
            % f = [n x m]

            
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
                    multiplyStripe = false;
                end

                if isstruct(z_smooth)
                    if any(size(z_smooth.Z)==1) % assume data need scattered interpolant
                        fun_scat = scatteredInterpolant(z_smooth.x,z_smooth.y,double(z_smooth.Z));
                        smooth_fun = @(x,y) fun_scat({y,x});
                    else % assume data is gridded
                        xg = z_smooth.x(:);
                        yg = z_smooth.y(:);
                        
                        z_smooth = double(z_smooth.Z);
                        
                        xn = (1:size(z_smooth,2)).';
                        yn = (1:size(z_smooth,1)).';

                        smooth_fun = @(x,y) interp2mex_wExpand(z_smooth, reshape(nakeinterp1(xg,xn,x(:)),size(x)), reshape(nakeinterp1(yg,yn,y(:)),size(y)));
                    end
                else
                    xg = obj.xCenters(:);
                    yg = obj.yCenters(:);
                    xn = (1:size(z_smooth,2)).';
                    yn = (1:size(z_smooth,1)).';
                    
                    z_smooth = double(z_smooth);
                    
                    smooth_fun = @(x,y) interp2mex_wExpand(z_smooth, reshape(nakeinterp1(xg,xn,x(:)),size(x)), reshape(nakeinterp1(yg,yn,y(:)),size(y)));
                end

                if isempty(x_stripe)
                    fun = @(x,y) cast( smooth_fun(x,y), outClass);
                else
                    x_stripe = double(x_stripe); % Ensure double
                    if multiplyStripe
                        fun = @(x,y) cast( x_stripe(x) .* smooth_fun(x,y), outClass); % cast final result to same class as image
                    else
                        fun = @(x,y) cast( x_stripe(x) + smooth_fun(x,y), outClass); % cast final result to same class as image
                    end
                end
            end

            % Note: interp2mex only performs nearest neighbor extrapolation
            % outside of the defined limits. Therefore, the half blockSize
            % on the image boarder could have jagged edges.

        end
    end
    
    methods (Static)
        th = otsuthresh_scale(I,scale);
        [Z,X,Y] = smoothSurface(x,y,z,r,type);
        [S,fun] = decimate_and_smooth(x,y,z,op);
    end
end % class
