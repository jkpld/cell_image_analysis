function features = extractGLCMFeatures(I,L,options)
% EXTRACTGLCMFEATURES Compute the 13 Haralick features from the gray
% co-occurance matrix for each object in an image.
%
%   - Assume single or double input with range [0,1]
%   - Computes textures over four directions, assumes symmetric, and
%   averages for the four directions.
%   
% Input
%   I : image
%   L : label matrix
%   options : structure with the following fields
%     D : offset distances
%     useGPU : Use GPU to speed up computation. (optional, defualt is
%       false) Using a GPU more than doubles the speed on my computer.
%
% Output
%   features : N x 13*ND array. N is the number of objects, 13 is the
%     number of features per offset distance,  and ND is the number of
%     offset distances.  The features are measured along [0,45,90,135]
%     degrees and averaged. 
%         contrast 
%         correlation 
%         differenceEntropy 
%         differenceVariance 
%         energy 
%         entropy 
%         informationMeasureOfCorrelation1 
%         informationMeasureOfCorrelation2 
%         inverseDifferenceMoment 
%         sumAverage 
%         sumEntropy 
%         sumOfSquaresVariance 
%         sumVariance 
%
% These are the same 13 Haralick features that CellProfiler outputs. The
% results match cellprofiler; hwoever, the information features matche
% CellProfiler 1, but not CellProfiler 2.

% James Kapaldo

D = options.Offset(:);
N_D = numel(D);
useGPU = options.Use_GPU;

NL = 8; % Use 8 levels.

L = cast(L,'like',I); % Make L the same class as I
numObjs = max(L(:)); % Number of objects

% Direction offsets
offsets = [0, -1, -1, -1; 1, 1, 0, -1];
offsets = reshape(offsets.*permute(D,[3,2,1]),2,4*N_D)';

features = zeros(numObjs,13*N_D); % Initialize features matrix

% If using GPU, then send the arrays to the GPU
if useGPU
    I = gpuArray(I);
    L = gpuArray(L);
    offsets = gpuArray(offsets);
end

% Compute the per object min and max to scale the intensity levels
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


% Compute the GLCM matrix and then the Harlick features for each direction
for i = 1:size(offsets,1)

    % Shift the images by the offset
    I2 = circshift(I,offsets(i,:));
    L2 = circshift(L,offsets(i,:));
    
    % Find the valid overlap regions
    valid = L==L2;

    % Get the intensities and object number for the valid region
    Inds = [I(valid), I2(valid), L(valid)];
    
    % Scale the intensities (per object) to be in the range [1,NL]
    Inds(:,1:2) = floor(NL*(Inds(:,1:2)-objMin(Inds(:,3)))./objRange(Inds(:,3)) + 1);
    
    % Compute the GLCM for each object
    GLCM = accumarray(Inds, 1, [NL, NL, numObjs]);
    
    % Make the GLCM symmetric
    GLCM = GLCM + permute(GLCM,[2,1,3]);
    
    if useGPU
        % If using a GPU, then bring the GLCM back to the CPU. The
        % GLCMFeatures function below is optimized for CPU, not GPU.
        GLCM = gather(GLCM);
    end

    % Add the features
    feat_idx = (1:13) + floor((i-1)/4)*13;
    features(:,feat_idx) = features(:,feat_idx) + GLCMFeatures(GLCM);
end

% Get the mean feature values over each direction
features = features/4;

% Set any NaN values to 0
features(isnan(features)) = 0;

end