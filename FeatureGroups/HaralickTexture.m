classdef HaralickTexture < FeatureGroup
    properties (Dependent, SetAccess = protected)
        FeatureNames
    end
    methods
        function obj = HaralickTexture(channel,options)
            obj.GroupName = class(obj);
            
            % Default options
            defaultOptions = struct('Offset',1,'Use_GPU',0);
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
            D = obj.Options.Offset;
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
            
            D = obj.Options.Offset(:);
            useGPU = obj.Options.Use_GPU;
            N_D = numel(D);
            
            
            NL = 8; % Use 8 levels.
            
            L = cast(L,'like',I); % Make L the same class as I
            numObjs = max(L(:)); % Number of objects
            
            % Direction offsets
            offsets = [0, -1, -1, -1; 1, 1, 0, -1];
            offsets = reshape(offsets.*permute(D,[3,2,1]),2,4*N_D)';
            
            x = zeros(numObjs,13*N_D); % Initialize features matrix
            
            % If using GPU, then send the arrays to the GPU
            if useGPU
                I = gpuArray(I);
                L = gpuArray(L);
                offsets = gpuArray(offsets);
            end
            
            % Compute the per object min and max to scale the intensity
            % levels
            BG = L==0; % background mask
            nBG = ~BG; % forground mask
            Lf = L(nBG); % forground labels
            If = I(nBG); % forground image
            objMin = accumarray(Lf,If,[numObjs,1],@min); % object min
            objMax = accumarray(Lf,If,[numObjs,1],@max) + 1e-6; % object max : this small offset ensures we do not have any NL+1 values after scaling (below in the for-loop)
            objRange = objMax - objMin; % object range
            
            L(BG) = NaN; % Set all non object pixels to NaN
            
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
                valid = L==L2;
                
                % Get the intensities and object number for the valid
                % region
                Inds = [I(valid), I2(valid), L(valid)];
                
                % Scale the intensities (per object) to be in the range
                % [1,NL]
                Inds(:,1:2) = floor(NL*(Inds(:,1:2)-objMin(Inds(:,3)))./objRange(Inds(:,3)) + 1);
                
                % Compute the GLCM for each object
                GLCM = accumarray(Inds, 1, [NL, NL, numObjs]);
                
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