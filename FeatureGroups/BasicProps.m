classdef BasicProps < FeatureGroup
    % BASICPROPS Compute the basic properties (centroid, integrated
    % intensity, and area) for each object in an image
    %
    % Properties :
    %
    % Channel - The channel from which the features will be extracted
    % requiredOptions - Structure with the following fields
    %   Use_GPU : (logical flag) Determine if to use gpu to speed up.
    %
    % Methods :
    %
    % Compute(I, L) - Compute the integrated intensity from image, I, for each
    % object in label matrix, L, and the area and centroid of each object.
    %   x = BasicProps.Compute(I, L)
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = BasicProps(channel,options)
            obj.GroupName = class(obj);
            
            % Default options
            defaultOptions = struct('Use_GPU',0);
            obj.requiredOptions = fieldnames(defaultOptions)';
            
            if nargin < 2
                options = defaultOptions;
            end
            obj.Options = options;
            
            if nargin >= 1
                obj.Channel = channel;
            end
        end
        
        function names = get.FeatureNames(obj)
            names = obj.GroupName + "_" + ...
                ["Area", ...
                 "Centroid_" + ["X","Y"], ...
                 "Intensity_" + obj.Channel + "_Integrated"];
        end
        
        function x = Compute(obj, I, L)
            % BASICPROPS.COMPUTE Return the area, centroid, and integrated
            % intensity for each object in an image I labeled with matrix
            % L.
            %
            % x = BasicProps.Compute(I, L)
            %
            % Input
            %  I : image
            %  L : object label matrix
            %
            % Output
            %  x : N x 4 array where N is the number of objects, and each
            %    column corresponds to FeatureNames
            
            % James Kapaldo
            
            useGPU = obj.Options.Use_GPU;
            
            FG = L~=0;
            
            nRow = size(L,1);
            linIdx = single(1:nmel(L))';
            
            I = single(I(FG));
            L = single(L(FG));
            linIdx = linIdx(FG);
            
            if useGPU
                I = gpuArray(I);
                L = gpuArray(L);
                linIdx = gpuArray(linIdx);
                nRow = gpuArray(nRow);
            end
            
            N_obj = max(L);
            
            % Row and column indices
            y = rem(linIdx-1, nRow) + 1;
            x = (linIdx - y)/nRow + 1; clear linIdx nRow;
            
            % Pixels per object
            if useGPU
                N = accumarray(L, gpuArray.ones(size(L),'single'), [N_obj, 1], @sum, single(1));
            else
                % Using sparse on CPU is faster than accumarray
                N = single(full(sparse(double(L),ones(size(L)), ones(size(L)), double(N_obj), 1)));
                N(N==0) = 1;
            end
            
            % Centroids
            x_mu = accumarray(L, x) ./ N;
            y_mu = accumarray(L, y) ./ N;
            
            % Integrated intensity
            S = accumarray(L, I, [N_obj, 1], @sum);
            
            x = [N, x_mu, y_mu, S];
            
            if useGPU
                x = gather(x);
            end
        end
    end
end