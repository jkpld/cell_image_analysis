classdef BasicProps < FeatureGroup
    % BASICPROPS Compute the basic properties (the integrated intensity and
    % area) for each object in an image
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
    % object in label matrix, L, and the area of each object.
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
            names = obj.GroupName + ...
                ["Area", obj.Channel + "_IntegratedIntensity"];
        end
        
        function x = Compute(obj, I, L)
            % BASICPROPS.COMPUTE Return the area and integrated intensity
            % for each object in an image I labeled with matrix L.
            %
            % x = BasicProps.Compute(I, L)
            %
            % Input
            %  I : image
            %  L : object label matrix
            %
            % Output
            %  x : N x 2 array where N is the number of objects, and each
            %    column corresponds to FeatureNames
            
            % James Kapaldo
            
            useGPU = obj.Options.Use_GPU;
            
            FG = L~=0;
            
            I = single(I(FG));
            L = single(L(FG));
            
            if useGPU
                I = gpuArray(I);
                L = gpuArray(L);
            end
            
            N_obj = max(L);
            
            % Pixels per object
            if useGPU
                N = accumarray(L, gpuArray.ones(size(L),'single'), [N_obj, 1], @sum, single(1));
            else
                % Using sparse on CPU is faster than accumarray
                N = single(full(sparse(double(L),ones(size(L)), ones(size(L)), double(N_obj), 1)));
                N(N==0) = 1;
            end
            
            % Integrated intensity
            S = accumarray(L, I, [N_obj, 1], @sum);
            
            x = [N, S];
            
            if useGPU
                x = gather(x);
            end
        end
    end
end