function features = extractGLCMFeatures(I,L,D,useGPU)
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
%   D : offset distance
%   useGPU : Use GPU to speed up computation. (optional, defualt is false)
%
% Output
%   features : N x 13, N is the number of objects :: N = max(L(:)) and 13
%   is the number of features. The features are measured along
%   [0,45,90,135] degrees and averaged.
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
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973


% II = imread('image_67_normalized.tif');
% % BW = imread('mask_67.tif');
% BW = BW_salrGeo;
% BW = bwareaopen(BW,100);
% I = single(I)/single(intmax(class(II)));
% L = bwlabel(BW);
% D = 7;

if nargin < 4
    useGPU = false;
end

NL = 8; % Use 8 levels.

numObjs = max(L(:));
L = cast(L,'like',I);

offsets = [0,D; -D,D; -D,0; -D,-D];
features = zeros(numObjs,13);

if useGPU
    I = gpuArray(I);
    L = gpuArray(L);
    offsets = gpuArray(offsets);
end

BG = L==0;
nBG = ~BG;
objMin = accumarray(L(nBG),I(nBG),[numObjs,1],@min);
objMax = accumarray(L(nBG),I(nBG),[numObjs,1],@max) + 1e-6;
objRange = objMax - objMin;

% I = floor(I*NL + 1); % scale image to have NL levels

% I(I > NL) = NL;
% I(I < 1) = 1;


L(BG) = NaN;
I = padarray(I,[D,D], 0);
L = padarray(L,[D,D], NaN);


for i = 1:size(offsets,1)

    I2 = circshift(I,offsets(i,:));
    L2 = circshift(L,offsets(i,:));
    
    valid = L==L2;
    
    Inds = [I(valid), I2(valid), L(valid)];
    Inds(:,1:2) = floor(NL*(Inds(:,1:2)-objMin(Inds(:,3)))./objRange(Inds(:,3)) + 1);
    
    GLCM = accumarray(Inds, 1, [NL, NL, numObjs]);
    GLCM = GLCM + permute(GLCM,[2,1,3]);
    
    if useGPU
        GLCM = gather(GLCM);
    end

    features = features + GLCMFeatures(GLCM);
end

features = features/4;