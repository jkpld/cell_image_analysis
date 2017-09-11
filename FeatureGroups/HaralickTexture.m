classdef HaralickTexture < FeatureGroup
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = HaralickTexture(channel,options)
            obj.GroupName = class(obj);
            
            % Default options
            defaultOptions = struct('PixelOffset',1,'Use_GPU',0);
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
            D = obj.Options.PixelOffset;
            featList = [...
                "Contrast", ...
                "Corr", ...
                "DiffEntropy", ...
                "DiffVar", ...
                "Energy", ...
                "Entropy", ...
                "InfoMeasCorr1", ...
                "InfoMeasCorr2", ...
                "InvDiffMoment", ...
                "SumAvg", ...
                "SumEntropy", ...
                "SumSquaresVar", ...
                "SumVar"]' + "_" + string(D(:)');
            
            featList = featList(:)';
            
            names = obj.GroupName + "_" + obj.Channel + "_" + featList;
        end
        
        function x = Compute(obj, I, L)
            % HARALICKTEXTURE.COMPUTE Return the 13 Haralick features from
            % the gray co-occurance matrix for each object in an image.
            %
            % x = HaralickTexture.Compute(obj, I, L)
            %
            % Input
            %   I : image - single/double in range [0,1]
            %   L : label matrix
            %
            % Output
            %   x : N x 13*N_Offset array. N is the number of objects, 13
            %     is the number of features per offset distance, and
            %     N_Offset = size(Offset) is the number of offset
            %     distances. The features are measured along [0,45,90,135]
            %     degrees and averaged. Each column of the returned
            %     features corresponds to the names given by FeatureNames.
            %
            % Note : Computes textures over four directions, assumes
            % symmetric co-occurance matrix, and averages for the four
            % directions.
            %
            % Note : These are the same 13 Haralick features that
            % CellProfiler outputs. The results match cellprofiler;
            % however, the information features match CellProfiler 1, but
            % not CellProfiler 2.
            
            % James Kapaldo
            
            D = obj.Options.PixelOffset(:);
            useGPU = obj.Options.Use_GPU;
            N_D = numel(D);
            
            
            NL = 8; % Use 8 levels.
            
            L = single(L); % Make L the same class as I
            N_obj = max(L(:)); % Number of objects
            
            % Direction offsets
            offsets = [0, -1, -1, -1; 1, 1, 0, -1];
            offsets = reshape(offsets.*permute(D,[3,2,1]),2,4*N_D)';
            
            x = zeros(N_obj,13*N_D); % Initialize features matrix
            
            % If using GPU, then send the arrays to the GPU
            if useGPU
                I = gpuArray(I);
                L = gpuArray(L);
                offsets = gpuArray(offsets);
            end
            
            % Compute the per object min and max to scale the intensity
            % levels
            BG = L==0; % background mask
            FG = ~BG; % forground mask
            Lf = L(FG); % forground labels
            If = I(FG); clear FG; % forground image
            
            objMin = [single(0); accumarray(Lf,If,[N_obj,1],@min)]; 
            objMax = [single(1); accumarray(Lf,If,[N_obj,1],@max) + 1e-6]; clear If Lf; % this small offset ensures we do not have any NL+1 values after scaling (below in the for-loop)
            objRange = objMax - objMin; clear objMax;
            
            % Scale intensities of each object to be in the range [1,NL].
            I = floor(NL*(I - objMin(L+1))./objRange(L+1) + 1); clear objMin objRange
            
            L(BG) = NaN; clear BG; % Set all non object pixels to NaN
            
            % Pad the arrays for circular shifting
            maxD = max(D);
            I = padarray(I,[maxD,maxD], 0);
            L = padarray(L,[maxD,maxD], NaN);
            
            
            % Compute the GLCM matrix and then the Harlick features for
            % each direction
            for i = 1:size(offsets,1)
                
                % Shift the images by the offset
                I2 = circshift(I,offsets(i,:));
                L2 = circshift(L,offsets(i,:));
                
                % Find the valid overlap regions
                valid = find(L==L2); % find() is much faster than the logical alternative.
                
                % Get the intensities and object number for the valid
                % region
                Inds = [I(valid), I2(valid), L(valid)];
                
                % Compute the GLCM for each object
                GLCM = accumarray(Inds, 1, [NL, NL, N_obj]);
                
                % Make the GLCM symmetric
                GLCM = GLCM + permute(GLCM,[2,1,3]);
                
                if useGPU
                    % If using a GPU, then bring the GLCM back to the CPU.
                    % The GLCMFeatures function below is optimized for CPU,
                    % not GPU.
                    GLCM = gather(GLCM);
                end
                
                % Add the features
                feat_idx = (1:13) + floor((i-1)/4)*13;
                x(:,feat_idx) = x(:,feat_idx) + GLCMFeatures(GLCM);
            end
            
            % Get the mean feature values over each direction
            x = x/4;
            
            % Set any NaN values to 0
            x(isnan(x)) = 0;
            
        end
    end
end