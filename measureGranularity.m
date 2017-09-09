function G = measureGranularity(I,mask,options)
% MEASUREGRANULARITY Compute the granular spectrum for each object in mask
% using the image I.
%
% G = measureGranularity(I,mask,options)
%
% Input
%   I : image, single or double between [0,1]
%   maks : object mask
%   options : structure with the following fields
%     SubSampleSize - the factor to re-size the image by before computing
%       the spectrum. 
%     BackgroundSampleSize - the factor to re-size the image before
%       computing the image background
%     ElementSize - The size of the primary elements
%     GranularSpectrumLength - The length of the granular spectrum
%     Use_GPU : (boolean) Determine if to use gpu to speed up. (speed up
%       factor is ~3 on my computer)
%
% Output
%   G : granular spectrum. N x GranularSpectrumLength, where N is the
%     number of objects. The object number corresponds to the object number
%     given by bwlabel/bwconncomp/...

% James Kapaldo

SubSampleSize = options.SubSampleSize;
BackgroundSampleSize = options.BackgroundSampleSize;
ElementSize = options.ElementSize;
GranularSpectrumLength = options.GranularSpectrumLength;
Use_GPU = options.Use_GPU;

sz = size(I);

mask = mask>0;

%% Subsample

if Use_GPU
    I = gpuArray(I);
    mask = gpuArray(mask);
end

if SubSampleSize ~= 1
    I = imresize(I, SubSampleSize);
    mask = imresize(single(mask), SubSampleSize)>0.5;
end

%% Background correct

% Resize image and mask
if BackgroundSampleSize ~= 1
    Ism = imresize(I, BackgroundSampleSize);
    mask_sm = imresize(single(mask), BackgroundSampleSize)>0.5;
else
    Ism = I;
    mask_sm = mask;
end

% Compute background
% - gpu erosion/dilation require uint8 and resizing requires single/double
%   so need to cast to different types
% - if using a gpu, then it could be necessary to conserve memory, so
% remove temporary variables when they are no longer needed.
se = strel('disk',ElementSize);
Ism8 = im2uint8(Ism);                   clear Ism
BGsm8 = imopen_mask(Ism8,se,~mask_sm);  clear Ism8 mask_sm
BGsm = single(BGsm8);                   clear BGsm8
BG = imresize(BGsm,sz);                 clear BGsm

% Set image to uint8 to make the same scale as BG and subtract BG
I = im2uint8(I) - uint8(BG);            clear BG

%% Label objects 
L = bwlabel(mask);

L = L(:);
FG = L>0;
L = L(FG);

%% Compute granular spectrum

% just need the ~mask from now on.
notmask = ~mask;                        clear mask;

ero = I; % copy the image to hold the erosions
N = max(L); % get the number of objects

% Create structure element
se = gpuArray([0 1 0; 1 1 1; 0 1 0]);

% Starting object mean intensity
Is = single(I); % need single array if on gpu
startmean = accumarray(L,Is(FG),[N,1],@sum);    clear Is
currentmean = startmean;

% initialize granular spectrum
if Use_GPU
    G = gpuArray.zeros(N,GranularSpectrumLength,'single');
else
    G = zeros(N,GranularSpectrumLength,'single');
end

% compute granular spectrum on gpu
for i = 1:GranularSpectrumLength
    prevmean = currentmean;
    ero = imerode_mask(ero,se,notmask);
    rec = imreconstruct(ero,I,4);
    
    if Use_GPU
        % If using a GPU then we need single arrays for accumarray.
        recS = single(rec);
        currentmean = accumarray(L,recS(FG),[N,1],@sum);
        clear recS
    else
        currentmean = accumarray(L,rec(FG),[N,1],@sum);
    end

    G(:,i) = prevmean - currentmean;
    clear rec
end

% normalize granular spectrum
G = 100*(G ./ startmean);

if Use_GPU
    G = gather(G);
end

end

function I = imopen_mask(I,se,notMask)

I = imerode_mask(I,se,notMask);
I = imdilate_mask(I,se,notMask);

end

function ero = imerode_mask(I,se,notMask)

% notMask = ~mask;
maskedI = I;
maskedI(notMask) = cast(Inf,'like',I); % Set not-mask to Inf so that it is not considered in the background correction
ero = imerode(maskedI,se);
ero(notMask) = I(notMask);

end

function dil = imdilate_mask(I,se,notMask)

maskedI = I;
maskedI(notMask) = cast(0,'like',I); % Set not-mask to Inf so that it is not considered in the background correction
dil = imdilate(maskedI,se);
dil(notMask) = I(notMask);

end
