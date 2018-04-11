classdef Location < FeatureGroup
    % LOCATION Compute the centroid and bounding box for each
    % object
    %
    % Properties :
    %
    % Channel - Not used. This FeatureGroup works directly on a mask.
    % requiredOptions - Structure with the following fields
    %   Use_GPU : logical flag determining if a GPU is used to speed up the
    %     computation.
    %
    % Methods :
    %
    % Compute(~, L) - Compute the centroid and bounding box for each object
    %   in label matrix, L.
    %   x = Location.Compute([], L)
    
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = Location(channel,options)
            obj.GroupName = class(obj);

            % Default options
            defaultOptions = struct('Use_GPU', 0);
            obj.requiredOptions = fieldnames(defaultOptions)';
            
            if nargin < 2
                options = defaultOptions;
            end
            obj.Options = options;
            
            if nargin >= 1
                if ~strcmp(channel, '') || ~isempty(channel)
                    warning('Location:channelNotUsed', 'The Location FeatureGroup does not operate on an image channel. Channel given will be ignored.')
                end
            end
        end
        
        function names = get.FeatureNames(obj)

            featList = [ ...
                "Centroid_X", ...
                "Centroid_Y", ...
                "BBox_" + ["L","T","W","H"]];
            
            names = obj.GroupName + "_" + featList;
        end
        
        function x = Compute(obj, ~, L)
            % LOCATION.COMPUTE Return the centroid and bounding box for
            % each object
            %
            % x = Location.Compute([],L)
            %
            % Input
            %   L : label matrix
            %
            % Output
            %   x : N x 6 array. N is the number of objects. The columns
            %     correspond to the features given by FeatureNames.
            
            % James Kapaldo
            
            useGPU = obj.Options.Use_GPU;
            
            % Initialize linear indices
            linIdx = (single(1):single(numel(L)))';
            
            % Get number of rows in image
            nRow = size(L,1);
            
            % Remove background pixels
            FG = L>0;
            L = single(L(FG));
            linIdx = linIdx(FG);
            
            % If using GPU, then send the arrays there
            if useGPU
                L = gpuArray(L);
                linIdx = gpuArray(linIdx);
                nRow = gpuArray(nRow);
            end
            
            N_obj = max(L); % number of objects
            
            % Row and column indices
            r = rem(linIdx-1, nRow) + 1;
            c = (linIdx - r)/nRow + 1; clear linIdx nRow;
            
            if useGPU
                % Pixels per object
                N = accumarray(L, gpuArray.ones(size(L),'single'), [N_obj, 1], @sum, single(1));
                
                % Initialize features array
                x = gpuArray.zeros(N_obj, 6, 'single');
            else
                % Pixels per object
                % Using sparse on CPU is faster than accumarray
                N = single(full(sparse(double(L),ones(size(L)), ones(size(L)), double(N_obj), 1)));
                N(N==0) = 1;
                
                % Initialize features array
                x = zeros(N_obj, 6, 'single');
            end
            
            % Centroids
            x(:,1) = accumarray(L, c) ./ N;
            x(:,2) = accumarray(L, r) ./ N; clear N;
            
            % Bounding Box
            minC = accumarray(L, c, [N_obj, 1], @min);
            maxC = accumarray(L, c, [N_obj, 1], @max); clear c;
            minR = accumarray(L, r, [N_obj, 1], @min);
            maxR = accumarray(L, r, [N_obj, 1], @max);
            
            x(:,3) = minC; % left corner
            x(:,4) = minR; % top corner
            x(:,5) = maxC - minC + 1; % width
            x(:,6) = maxR - minR + 1; % height
            
            if useGPU
                x = gather(x);
            end
        end
    end
end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
