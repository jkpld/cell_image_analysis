classdef Intensity < FeatureGroup
    % INTENSITY Compute the intensity features for each object in an image.
    %
    % Properties :
    %
    % Channel - The channel from which the features will be extracted
    % requiredOptions - Structure with the following fields
    %   Use_GPU : (logical flag) Determine if to use gpu to speed up.
    %
    % Methods :
    %
    % Compute(I, L) - Compute intensity features from image, I, for each
    % object in label matrix, L.
    %   x = Intensity.Compute(I, L)
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = Intensity(channel,options)
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
            featList = ["Integrated", ...
                "Mean", ...
                "Std", ...
                "Skewness", ...
                "Kurtosis"];
            
            names = obj.GroupName + "_" + obj.Channel + "_" + featList;
        end
        
        function x = Compute(obj, I, L)
            % INTENSITY.COMPUTE Return the integrated, mean, standard
            % deviation, skewness, and kurtosis of the intensities for each
            % object in an image I labeled with matrix L.
            %
            % x = Intensity.Compute(I, L)
            %
            % Input
            %  I : image
            %  L : object label matrix
            %
            % Output
            %  x : N x 5 array where N is the number of objects, and each
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
            
            % NOTE : accumarray on gpu does not accept @mean, @std, or the
            % others; thus, they need to be calculated out by hand. This is
            % actually faster than using the other commands anyway, because
            % many of the commands calculate the same thing.
            
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
            
            % Mean
            mu = S./N;
            
            % expand out the object mean to each pixel
            I_mu = mu(L);
            
            % Standard deviation
            m2_i = (I - I_mu).^2;
            m2 = accumarray(L,m2_i, [N_obj, 1], @sum); clear m2_i
            s = sqrt(m2./(N-1));
            m2 = m2 ./ N;
            
            % Skewness
            m3_i = (I - I_mu).^3;
            m3 = accumarray(L,m3_i, [N_obj, 1], @sum) ./ N; clear m3_i
            skw = sqrt((N-1)./N) .* N./(N-2) .* m3 ./ m2.^(1.5); % unbiased skewness
            skw(N<3) = NaN;
            
            % Kurtosis
            m4_i = (I - I_mu).^4;
            m4 = accumarray(L,m4_i, [N_obj, 1], @sum) ./ N; clear m4_i
            krt = ((N-1)./((N-2).*(N-3)) .* ( (N+1) .* m4 ./ m2.^2 - 3*(N-1) ));
            krt(N<4) = NaN;
            
            x = [S, mu, s, skw, krt];
            
            if useGPU
                x = gather(x);
            end
            
        end
    end
end