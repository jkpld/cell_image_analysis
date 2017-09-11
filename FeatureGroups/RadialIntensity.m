classdef RadialIntensity < FeatureGroup
    % RADIALINTENSITY Compute the mean and standard deviation of object
    % intensity values in concentric rings.
    %
    % Properties :
    %
    % Channel - The channel from which the features will be extracted
    % requiredOptions - Structure with the following fields
    %   Number_of_Rings : The number of rings to compute
    %   Ring_Width : The width of the rings
    %   Use_GPU : logical flag determining if a GPU is used to speed up the
    %     computation.
    %
    % Methods :
    %
    % Compute(I, L) - Compute the mean and std of the intensity from image,
    %   I, for concentric rings from each object in label matrix, L.
    %   x = RadialIntensity.Compute(I, L)
    
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = RadialIntensity(channel,options)
            obj.GroupName = class(obj);

            % Default options
            defaultOptions = struct(...
                'Number_of_Rings', 1, ...
                'Ring_Width', 1,...
                'Use_GPU', 0);
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
            NR = obj.Options.Number_of_Rings;
            W = obj.Options.Ring_Width;
            featList = (["Mean_", "Std_"] + "W" + string(W) + "_") ...
                + ["Center", "R" + string(NR:-1:1)]';
            featList = featList(:)'           ;
            
            names = obj.GroupName + "_" + obj.Channel + "_" + featList;
        end
        
        function x = Compute(obj, I, L)
            % RADIALINTENSITY.COMPUTE Return the mean and standard
            % deviation of object intensity values in concentric rings.
            %
            % x = RadialIntensity.Compute(I, L)
            %
            % Input
            %  I : image
            %  L : label matrix
            %
            % Output
            %  x : N x 2*(Number_of_Rings + 1) array giving the mean and
            %    std intensity value for each ring. N is the number of
            %    objects. (The +1 is explained below.) The columns
            %    correspond to the FeatureNames
            %
            % Note 1. This function computes Number_of_Rings rings with
            % width Ring_Width starting from the boundary of the objects.
            % The mean and standard devation of any region in the center of
            % the object that is not included in a ring is also returned;
            % this is why there are (Number_of_Rings + 1) columns in the
            % output.
            %
            % Note 2. Ring 1 is the outermost ring
            %
            % Note 3. If an object is small, then the inner rings and
            % central region might be empty. Any empty inner ring or
            % central region will have a mean and standard deviation of 0
            % returned.
            %
            % Method: This function computes the pixels belonging to each
            % ring by thresholding the distance transform of the objects.
            
            % James Kapaldo
            
            N_R = 1 + obj.Options.Number_of_Rings; % number of rings; +1 will hold any leftover region at the center of the object
            W = obj.Options.Ring_Width; % width of the N_R rings
            useGPU = obj.Options.Use_GPU;
            
            BW = L>0; % compute binary mask
            BWi = find(BW);
            % Note: below I clear all variables that are no longer needed. This can be
            % helpful when working on a GPU.
            I = single(I(BWi));
            L = single(L(BWi));
            
            % Compute distance transform. This is used for computing the pixels in each
            % ring.
            %    NOTE : There is some error with computing the dist. trans. on the GPU
            %    with large images (1000x1000 works but 1500x1500 does not). Therefore,
            %    compute it on the CPU, even though this is much slower.
            %    TODO : Try to get the distance transform to work on the GPU.
            D_full = bwdist(~BW);
            D = single(D_full(BWi)); clear D_full;
            
            if useGPU
                % If using gpu, then send over the arrays.
                I = gpuArray(I);
                L = gpuArray(L);
                D = gpuArray(D);
            end
            
            N_obj = max(L); % The number of objects
            N_D = numel(D); % The number of object pixels
            
            % Create a temporary array that holds the logical indices for layer
            if useGPU
                layer_idx = [gpuArray.false(N_D,1), D > ( N_R-1:-1:1 ) * W, gpuArray.true(N_D,1)]; clear D
            else
                layer_idx = [false(N_D,1), D > ( N_R-1:-1:1 ) * W, true(N_D,1)]; clear D
            end
            
            % The pixels belonging to each ring are the difference between the upper
            % layer and the current layer.
            ring_idx = layer_idx(:,2:end) & ~layer_idx(:,1:end-1); clear layer_idx;
            
            % Convert the logical indices to linear indices
            ring_linidx_full = ring_idx .* (1:N_D)';
            
            % Compute the ring number of each index
            ring_num = ring_idx .* (1:N_R);
            
            % Remove all 0 indices
            ring_linidx = ring_linidx_full(ring_linidx_full ~= 0); clear ring_linidx_full;
            
            % Compute the indices for the final output (first column is the
            % object number, second column is the ring number. (Note that
            % the first ring is on the inside of the object and the last
            % ring is on the outside of the object.)
            inds = [ L(ring_linidx), ring_num(ring_idx)]; clear ring_idx ring_num L;
            
            % Compute the values of each pixel. Note vals should be the
            % same size as I, just re-ordered. (The same goes for L above.)
            vals = I(ring_linidx); clear ring_linidx I;
            
            % Now, compute the mean and standard deviation. Note that the
            % mean and standard deviation functions are not available for
            % accumarray on a GPU; thus, I calculate them the long way.
            
            % Compute the number of pixels in each ring of each object. The
            % default value is 1, to prevent the NaN when computing the
            % mean.
            
            if useGPU
                N = accumarray(inds, gpuArray.ones(numel(vals),1,'single'), [N_obj, N_R], @sum, single(1));
            else
                % N = accumarray(inds, ones(numel(vals),1,'single'), [N_obj, N_R], @sum, single(1));
                N = single(full(sparse(double(inds(:,1)), double(inds(:,2)), ones(numel(vals),1), double(N_obj), N_R)));
                N(N==0) = 1;
                
                % For some reason using sparse instead of accumarray here
                % is faster for computing N, but it is slower for computing
                % the other quantities that use accumarray below.
            end
            
            % Compute the sum of the intensities in each ring of each
            % object
            S = accumarray(inds, vals, [N_obj, N_R], @sum);
            
            % Compute the mean intensity in each ring of each object
            Mu = S ./ N;
                
            % Compute the square difference between the intensity and mean
            % intensity
            Var_i = (vals - Mu(inds(:,1) + (inds(:,2)-1)*N_obj)).^2; clear vals

            % Compute the variance in each ring for each object
            Var = accumarray(inds, Var_i, [N_obj, N_R], @sum); clear Var_i inds

            % Compute the standard deviation
            Sigma = sqrt(Var./(N-1)); clear Var

            % Replace any Inf's with 0.
            Sigma(N==1) = 0;
            
            if useGPU
                % If using gpu, then gather the results
                Mu = gather(Mu);
                Sigma = gather(Sigma);
            end
            x = [Mu, Sigma];
        end
    end
end