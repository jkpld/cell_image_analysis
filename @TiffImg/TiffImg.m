classdef TiffImg < matlab.mixin.Copyable

    properties
        FileName(1,:) char
    end
    properties (SetAccess = private)
        Acquisition_Info(1,1) struct
        mmPerPixel(1,1) double

        workingClass(1,:) char
        imageClass(1,:) char
        imageSize(1,2) double
        stripeWidth(1,1) double
    end

    properties
        blockSize(1,1) double

        % Surface_Smoothing_Radius - The radius (in mm) to use when
        % smoothing the surfaces computed (threshold, background,
        % foreground). If set to NaN, then smoothing will be disabled.
        Surface_Smoothing_Radius(1,1) double = NaN;

        % Iamge_Smooth_Kernel - The Kernel used to smooth the image before
        % all operations (directly after reading the image in.)
        Image_Smooth_Kernel = fspecial('gaussian',7,1);

        % Image_Correction_Expression : The expression used to correct the
        % image. Valid variable names are
        %    S : The input image
        %    BG_o : A background offset (BG_offset) that is added or
        %      subtracted
        %    BG_s : A background x stripe (BG_stripeX) that is multiplied
        %      or divided
        %    FG_o : A foreground offset (FG_offset) that is added or
        %      subtracted
        %    FG_f : A foreground factor (FG_factor) that is multiplied or
        %      divided
        %    FG_s : A foreground x stripe (FG_stripeX) that is multiplied
        %      or divided
        %
        %  If BG_offset or FG_offset have not been assigned to the
        %  TiffImg yet, but they show up in the expression, they are
        %  replaced with 0's.
        %  If BG_factor, FG_factor, BG_stripeX, or FG_stripeX have not
        %  been assigned to the TiffImg yet, but they show up in the
        %  epxression, they are replaced with 1's.
        %
        %   Example expression (this is also the default expression):
        %    (S - BG_o * BG_s + FG_o) / (FG_s * FG_f)
        Image_Correction_Expression(1,1) string {validateExpression(Image_Correction_Expression,"S")} = "S"
    end
    properties (Dependent)
        Current_Image_Correction_Expression
    end
    properties
        ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression(1,1) string {validateExpression} = ""
    end

    properties (SetAccess = private, Hidden)
        tileSize
        tilesPerBlck
        numBlcks
        numTiles
    end

    properties (SetAccess = private, Hidden)
        FileID

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

        Threshold_After_Correction(1,1) logical = false;
        Threshold_CorrectionsExpression;
    end
    properties (Hidden)
        threshold_fun = []
        Threshold_Correction = []
        Secondary_Correction = []
    end
    properties
        threshold = []

        % old props
%         BG_smooth = []
%         BG_Xstripe = [] % this multiplies the BG_offset
%         FG_smooth = []
%         FG_Xstripe = [] % this multiplies the FG_factor

        % new props
%         BG_factor = []
        BG_offset = []
        BG_stripeX = []

        FG_offset = [] % to subtract from the foreground
        FG_factor = [] % to divide the foreground
        FG_stripeX = []


        Sharpness = []

        Use_GPU(1,1) logical = false;
        Verbose(1,1) logical = false;
    end
    properties
        UserData
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

            aqstnInfo = aqstnInfo{1};
            proc_info = regexp(aqstnInfo,'(?<key>[^=].*?) = (?<value>.*)','names');
            notKeyValue = cellfun(@isempty,proc_info);
            info.Extra = aqstnInfo(notKeyValue);
            for i = 1:numel(proc_info)
                if notKeyValue(i)
                    continue;
                else
                    valueNumber = str2double(proc_info{i}.value);
                    if isnan(valueNumber)
                        info.(strrep(proc_info{i}.key,' ','_')) = proc_info{i}.value;
                    else
                        info.(strrep(proc_info{i}.key,' ','_')) = str2double(proc_info{i}.value);
                    end
                end
            end

            if ~isfield(info,'MPP') || ~isfinite(info.MPP)
                tifflib('close',obj.FileID);
                obj.FileName = '';
                obj.FileID = uint64(0);
                error('TiffImg:MissingMMPerPixel','The image description does not have a MPP (micron per pixel) field, which is required for this class.')
            end

            obj.Acquisition_Info = info;
            obj.mmPerPixel = info.MPP/1000;


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
                    error('TiffImg:unsupporetedFormat','The image has an unsupported format, Supported formats are uint, int, and floating point.')
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


            if ~isfield(info,'StripeWidth') || ~isfinite(info.StripeWidth)
                str = ['The image description does not have a StripeWidth ',...
                    'field giving the width, in pixels, of each line scan ',...
                    'of the microscope.\nIf the image has an stripe ',...
                    'artifact along the scan direction, then it could be ',...
                    'important that the blockSize used to compute the ',...
                    'Threshold/Background/... is approximately the same ',...
                    'size as this StripeWidth, or a multiple of it. If it ',...
                    'is not, then artifacts could appear in the ',...
                    'Threshold/Background/...'];
                warning('TiffImg:MissingStripeWidth',str)

                obj.stripeWidth = nan;
                obj.blockSize = max(obj.tileSize);
            else
                obj.stripeWidth = info.StripeWidth;

                % Set the default block size to be the nearest tileSize
                % multiple of the StripeWidth
                obj.blockSize = round(obj.stripeWidth/max(obj.tileSize))*max(obj.tileSize);
            end


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

            if numel(blockSize)==2
                obj.xEdges = [0:floor(obj.imageSize(2)/blockSize(2)), obj.imageSize(2)/blockSize(2)]*blockSize(2);
                obj.yEdges = [0:floor(obj.imageSize(1)/blockSize(1)), obj.imageSize(1)/blockSize(1)]*blockSize(1);
            else
                obj.xEdges = [0:floor(obj.imageSize(2)/blockSize), obj.imageSize(2)/blockSize]*blockSize;
                obj.yEdges = [0:floor(obj.imageSize(1)/blockSize), obj.imageSize(1)/blockSize]*blockSize;
            end


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

        function exprssn = get.Current_Image_Correction_Expression(obj)
            exprssn = processExpression(obj,obj.Image_Correction_Expression,true,"S");
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

        function [I,x,y] = getTile(obj, tile_idx)
            % initialize array for image stripe
            I = zeros(obj.tileSize,obj.imageClass);

            tmp = tifflib('readEncodedTile',obj.FileID,double(tile_idx));

            [m,n] = size(tmp);
            I(1:m, 1:n) = tmp;

            % Convert image to working class.
            if ~strcmp(obj.imageClass,'logical')
                I = cast(I,obj.workingClass);
                I = I / cast(obj.maxSampleValue,obj.workingClass);
            end

            if nargout > 1
                [j,k] = ind2sub(size(obj.tiles),tile_idx+1);
                y = obj.tile_y_inds + (j-1)*obj.tileSize(1);
                x = obj.tile_x_inds + (k-1)*obj.tileSize(2);
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
                %  clear tmp
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

        function tiles_xy = computeTileNumbers(obj,bboxes)
            % COMPUTETILENUMBERS Compute the tile numbers that contain the
            % bounding boxes given.
            %
            % tiles = computeTileNumbers(obj,bboxes)
            %
            % Input
            %   bboxes : Nx4 list of bounding boxes
            %
            % Output
            %   tiles : Nx4 array where the i'th row is [x1,y1,x2,y2].
            %     x1/x2 is the left/right-most tiles. y1/y2 is the
            %     top/bottom-most tiles. These are not the linear indices
            %     that can be used to read in a tile with tifflib; they
            %     must be converted to linear indices. The x/y indices are
            %     0 based.

            % compute x,y points of bounding box corners
            %  -> bounding box = (left, top, width, height)
            bboxes(:,3:4) = bboxes(:,3:4) + bboxes(:,1:2);

            % x,y tile numbers, 0 based
            tiles_xy = ceil(bboxes ./ [obj.tileSize,obj.tileSize]) - 1;
        end

        Z = smoothSurf(obj,z,type);

        th = Compute_Threshold(obj, applyCorrectionsFirst);

        bg = Compute_Background(obj, computeXStripeArtifact, varargin);

        fg = Compute_Foreground(obj, computeXStripeArtifact, varargin);

        st = Compute_StripeArtifact(obj, BG_or_FG, varargin);

        sh = Compute_Sharpness(obj);

        [x, x_names] = Measure_BasicProps(obj,channelName, varargin);

        feature = Measure_Intensity(tiffImg, Use_Parallel, varargin);

        [FG, BG, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(tiffImg,x,y,dapi,area,options);

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

            valid = {'threshold','background','foreground_factor','foreground_offset','sharpness','stripe_background','stripe_foreground'};
            internal_name = {'threshold','BG_offset','FG_factor','FG_offset','Sharpness','BG_stripeX','FG_stripeX'};
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
                if idx > 5
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
                    shading interp
                    title(strrep(name_user,'_',' '))
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

        function fun = generateFunction(obj,expression,removeUndefined,expandInput,requiredVars)
            % GENERATEFUNCTION Create a function that operates on an image
            % based on a string expression
            %
            % fun = generateFunction(obj,expression,expandInput)
            %
            % Input
            %   expression : A string expression, that follows that same
            %     requirments as Image_Correction_Expression
            %   removeUndefined : Logical flag. If true, then any
            %     not-yet-defined variables will be removed from the
            %     expression. If false, then if there are any
            %     not-yet-defined variables in the expression, an error
            %     will be issued.
            %   expandInput : Logical flag. If true, then the input arrays
            %     to evaluate the correction function at, xi and yi, are
            %     expanded to form a matrix of size size(yi,1) x
            %     size(xi,2). If false, then the input arrays are not
            %     expanded, and the output is size(xi). Default value is
            %     true.
            %   requiredVars : (optional) A string array listing inputs
            %     required to show up in the expression. Default is none,
            %     []. These are the first elements to the function handle.
            %     Examples: requiredVars = ["S"] -> @(S,x,y). requiredVars
            %     = ["S","A"] -> @(S,A,x,y). requiredVars = {} -> @(x,y)
            %
            % Ouput
            %   fun : A function handle with three inputs (S,xi,yo). S is
            %     the image to evaluate on, xi and yi are the pixel values
            %     of the image.
            %
            % Example call
            %   img_corrected = fun(img,xi,yi);

            if isempty(expression) || expression == ""
                fun = [];
                return;
            end

            if nargin < 5
                requiredVars = "";
            end

            % Process the expression
            exprssn = processExpression(obj,expression,removeUndefined,requiredVars);

            % Generate the function
            fun = generateFunctionFromExpression(obj, exprssn, expandInput, requiredVars);
        end

        function fun = generateCorrectionFunction(obj, expandInput)
            % GENERATECORRECTIONFUNCITON Create function handle to apply
            % the Current_Image_Correction_Expression
            %
            % fun = generateCorrectionFunction(obj)
            % fun = generateCorrectionFunction(obj, expandInput)
            %
            % Input
            %   expandInput : (optional) Logical flag. If true, then the
            %     input arrays to evaluate the correction function at, xi
            %     and yi, are expanded to form a matrix of size size(yi,1)
            %     x size(xi,2). If false, then the input arrays are not
            %     expanded, and the output is size(xi). Default value is
            %     true.
            %
            % Ouput
            %   fun : A function handle with three inputs (S,xi,yo). S is
            %     the image to evaluate on, xi and yi are the pixel values
            %     of the image.
            %
            % Example call
            %   img_corrected = fun(img,xi,yi);


            if nargin < 2
                expandInput = true;
            else
                expandInput = logical(expandInput(1));
            end

            fun = generateFunctionFromExpression(obj, obj.Current_Image_Correction_Expression, expandInput, "S");
        end

        function clearCorrections(t)
            t.BG_offset = [];
            t.BG_stripeX = [];
            t.FG_offset = [];
            t.FG_factor = [];
            t.FG_stripeX = [];
        end
    end

    methods (Access = private, Hidden)
        function fun = generateFunctionFromExpression(obj, expression, expandInput, requiredVars)
            % GENERATECORRECTIONFUNCITON Create function handle to apply
            % the Current_Image_Correction_Expression
            %
            % fun = generateCorrectionFunction(obj)
            % fun = generateCorrectionFunction(obj, Name1, Value1, ...)
            %
            % Input
            %   expression : Cell array where the first element is the
            %     string expression and the second element is the list of
            %     variables in the expression.
            %   expandInput : (optional) Logical flag. If true, then the
            %     input arrays to evaluate the correction function at, xi
            %     and yi, are expanded to form a matrix of size size(yi,1)
            %     x size(xi,2). If false, then the input arrays are not
            %     expanded, and the output is size(xi). Default value is
            %     true.
            %
            % Ouput
            %   fun : A function handle with three inputs (S,xi,yo). S is
            %     the image to evaluate on, xi and yi are the pixel values
            %     of the image.
            %
            % Example call
            %   img_corrected = fun(img,xi,yi);

            if nargin < 3
                expandInput = true;
            else
                expandInput = logical(expandInput(1));
            end

            currentVars = expression{2};
            currentExpr = expression{1};

            intpl1d = @(x,y,xi) cast(nakeinterp1(double(x(:)), double(y(:)), double(xi(:))),'like',y);

            % Create the functions for the current expression
            for i = 1:numel(currentVars)
                switch currentVars{i}
                    case 'BG_o'
                        BG_o = interpolator2d(obj.BG_offset.x, obj.BG_offset.y, obj.BG_offset.Z, expandInput); %#ok<NASGU>
                        %                     case 'BG_f'
                        %                         BG_f = interpolator2d(obj.BG_factor.x, obj.BG_factor.y, obj.BG_factor.Z, expandInput); %#ok<NASGU>
                    case 'BG_s'
                        xd = (1:obj.imageSize(2)).';
                        if expandInput
                            BG_s = @(x) reshape(intpl1d(xd,obj.BG_stripeX,x),size(x)); %#ok<NASGU>
                        else
                            BG_s = @(x) intpl1d(xd,obj.BG_stripeX,x); %#ok<NASGU>
                        end
                    case 'FG_o'
                        FG_o = interpolator2d(obj.FG_offset.x, obj.FG_offset.y, obj.FG_offset.Z, expandInput); %#ok<NASGU>
                    case 'FG_f'
                        FG_f = interpolator2d(obj.FG_factor.x, obj.FG_factor.y, obj.FG_factor.Z, expandInput); %#ok<NASGU>
                    case 'FG_s'
                        xd = (1:obj.imageSize(2)).';
                        if expandInput
                            FG_s = @(x) reshape(intpl1d(xd,obj.FG_stripeX,x),size(x)); %#ok<NASGU>
                        else
                            FG_s = @(x) intpl1d(xd,obj.FG_stripeX,x); %#ok<NASGU>
                        end
                    case cellstr(requiredVars)
                    otherwise
                        error('createEvaluator:unknownVar','There is an unknown variable. If you see this error, then the createEvaluator() function is out of sync with the validateExpression() and get.Current_Image_Correction_Expression() functions.')
                end
            end

            % Evaluate the current expression and output the resulting
            % function handle.

            fun = eval(currentExpr);
        end

        function exprssn = processExpression(obj,expression,removeUndefined,requiredVars)
            % PROCESSEXPRESSION Process an expression an output a
            % processed (and, if requested, simplified) expression with the
            % list of variables in the expression.
            %
            % fun = processExpression(obj,expression,removeUndefined)
            %
            % Input
            %   expression : A string expression, that follows that same
            %     requirments as Image_Correction_Expression
            %   removeUndefined : Logical flag. If true, then any
            %     not-yet-defined variables will be removed from the
            %     expression. If false, then if there are any
            %     not-yet-defined variables in the expression, an error
            %     will be issued.
            %   requiredVars : (Optional). A string array of variable inputs
            %     that are required to appear as inputs. These are the
            %     first elements to the function handle. Examples:
            %     requiredVars = {"S"} -> @(S,x,y). requiredVars =
            %     {"S","A"} -> @(S,A,x,y). requiredVars = {} -> @(x,y)
            %
            % Ouput
            %   exprssn : A 1 x 2 cell array where the first element is the
            %     processed expression and the second element is the list
            %     of variables.
            %
            %
            % Example:
            %   input expression = "S - BG_o * BG_s"
            %   input removeUndefined = true;
            %
            %   (Assume BG_s has not been defined)
            %
            %   output exprssn = {"@(S,x,y) S - BG_o(x,y)", ["S","BG_o"]}

            if nargin < 4
                requiredVars = "";
            end

            % Validate expression
            validateExpression(expression,requiredVars);
            requiredVars = requiredVars(requiredVars~="");

            exprssn = str2sym(char(expression));

            % Get a list of the variables in the expression
            vars = string(symvar(char(expression)))';

            % If a variable from the expression has not been defined, then
            % throw and error.
            % names = {'BG_offset','BG_factor','BG_stripeX','FG_offset','FG_factor','FG_stripeX'};
            names = {'BG_offset','BG_stripeX','FG_offset','FG_factor','FG_stripeX'};
            for i = 1:numel(vars)
                if contains(vars(i),requiredVars)
                    continue;
                else
                    idx = startsWith(names,vars(i));
                    if isempty(obj.(names{idx}))
                        if removeUndefined
                            % variable name
                            name = names{idx}(1:4);
                            if strcmp(name(4),'o')
                                % If an offset, then replace with 0
                                exprssn = subs(exprssn,name,0);
                            else
                                % If a factor, then replace with 1
                                exprssn = subs(exprssn,name,1);
                            end
                        else
                            error('generateFunction:missingVariable','The property, %s, has not been defined yet.', names{idx})
                        end
                    end
                end
            end

            % Get an updated list of the variables in the expression if
            % necessary
            if removeUndefined
                vars = string(symvar(char(exprssn)))';
            end

            % Convert to symbolic expression and back to string, to make
            % sure that the expression is valid.
            exprssn = string(exprssn);

            % Vectorize the expression
            exprssn = vectorize(exprssn);

            % Iterate over each variable still in the expression and
            % add the function inputs
            for i = 1:numel(vars)
                if contains(vars(i),requiredVars)

                    continue;
                else
                    if endsWith(vars(i),'s')
                        exprssn = strrep(exprssn,vars(i),vars(i) + '(x)');
                    else
                        exprssn = strrep(exprssn,vars(i),vars(i) + '(x,y)');
                    end
                end
            end

            % Tack on the function handle head
            if contains(exprssn,"x")
                hasX = "x";
            else
                hasX = "~";
            end
            if contains(exprssn,"y")
                hasY = "y";
            else
                hasY = "~";
            end
            exprssn = "@(" + join([string(requiredVars),hasX,hasY],',') + ")" + exprssn;

            % Output the expression and variables as a cell array
            exprssn = {exprssn, vars};
        end
    end

    methods (Static)
        th = otsuthresh_scale(I,scale);
        [Z,X,Y] = smoothSurface(x,y,z,r,type);
        [S,fun] = decimate_and_smooth(x,y,z,op);
    end
end % class

function validateExpression(exprssn,requiredVars)
vars = string(symvar(char(exprssn)));

if nargin < 2
    requiredVars = "";
end

if any(requiredVars == "I")
    error('validateExpression:invalidVariableName','The varialbe "I" cannot be used')
end

% validVars = ["S","BG_o","BG_f","BG_s","FG_o","FG_f","FG_s"];
validVars = [requiredVars, "BG_o","BG_s","FG_o","FG_f","FG_s"];
if ~all(ismember(vars,validVars))
    error('validateExpression:invalidVars','Expression can only contain the variables %s.', join(validVars,', '))
end
if nargin > 1 && requiredVars ~= ""
    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars(i),vars)
            error('validateExpression:missingS','The variable %s must appear in the expression!', requiredVars(i))
        end
    end
end
end
