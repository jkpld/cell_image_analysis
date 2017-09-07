function G = measureGranularity(I,mask,options)

SubSampleSize = options.SubSampleSize;
BackgroundSampleSize = options.BackgroundSampleSize;
ElementSize = options.ElementSize;
GranularSpectrumLength = options.GranularSpectrumLength;

I = double(I);
mask = mask>0;

%% Subsample
if SubSampleSize ~= 1
    I = imresize(I,SubSampleSize,'bilinear');
    mask = imresize(mask,SubSampleSize,'bilinear')>0.9;
end

%% Background correct

% Create background structuring element.
%   -> instead of resizing the image, use a larger structering element
% se = strel('disk',round(ElementSize*(1/BackgroundSampleSize)));

% I_BG = imopen_mask(I,se,mask);
% I_FG = I - I_BG;
% I(mask) = I_FG(mask);


% Resize image and mask
if BackgroundSampleSize ~= 1
    Ism = imresize(I,BackgroundSampleSize,'bilinear');
    mask_sm = imresize(mask,BackgroundSampleSize,'bilinear')>0.5;
else
    Ism = I;
    mask_sm = mask;
end
mask_sm_d = gpuArray(mask_sm);

% Remove background
se = strel('disk',ElementSize);
minIsm = min(Ism(:));
maxIsm = max(Ism(:));
Ism = uint8(255*(Ism-minIsm)/(maxIsm-minIsm));
Ism_d = gpuArray(Ism);
Ism_BG_d = imopen_mask(Ism_d,se,~mask_sm_d);
I_FG = Ism - gather(Ism_BG_d); %*(maxIsm-minIsm)/255 + minIsm;
clear Ism_d Ism_BG_d mask_sm_d

% Un-resize
I = imresize(I_FG,size(I),'bilinear');

%% Convert image to uint8

% I(~mask) = 0;% = I.*mask;
% I = I - min(I(mask));
% I = uint8(255*I/max(I(mask)));


%% Label objects

% Label objects
% CC = bwconncomp(mask);
% L = labelmatrix(CC);
L = bwlabel(mask);

L = L(:);
FG = L>0;
L = L(FG);


%% Compute granular spectrum



% Send image to GPU
I_d = gpuArray(I);
ero_d = gpuArray(I);
notmask_d = gpuArray(~mask);
% L_d = gpuArray(L);
% FG_d = gpuArray(FG);
% N_d = gpuArray(CC.NumObjects);

% I_d = I;
% ero_d = I;
% ero_d(~mask) = 0;
% mask_d = mask;

N = max(L);

% Create structure element
se_d = gpuArray([0 1 0; 1 1 1; 0 1 0]);%strel('diamond',1);
% se = strel('diamond',1);

% Starting object mean intensity
startmean = accumarray(L,I(FG),[N,1],@sum);
currentmean = startmean;

% initialize granular spectrum
G = zeros(N,GranularSpectrumLength);

% compute granular spectrum on gpu
for i = 1:GranularSpectrumLength
    prevmean = currentmean;
    ero_d = imerode_mask(ero_d,se_d,notmask_d);
    rec_d = imreconstruct(ero_d,I_d,4);
    rec = gather(rec_d);
%     rec = rec_d;

    currentmean = accumarray(L,rec(FG),[N,1],@sum);
%     currentmean_d = accumarray(L_d,single(rec_d(FG_d)),[N_d,1],@sum);
%     G(:,i) = prevmean - gather(currentmean_d);
    G(:,i) = prevmean - currentmean;
end

% clean up memory
clear ero_d I_d rec_d se_d notmask_d

% normalize granular spectrum
G = 100*(G ./ startmean);


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
