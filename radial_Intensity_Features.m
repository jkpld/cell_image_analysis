function [Mu, Sigma] = radial_Intensity_Features(BW, I, options)
% RADIAL_INTENSITY_Features Compute the mean and standard deviation of
% object intensity values in concentric rings.
%
% Mu = radial_Intensity_Features(BW, I, options)
% [Mu, Sigma] = radial_Intensity_Features(BW, I, options)
%
% Input
%  BW : binary mask of objects
%  I : image
%  options : structure with the following fields
%    Number_of_Rings : The number of rings to compute
%    Ring_Width : The width of the rings
%    Use_GPU : boolian flag determining if a GPU is used to speed up the
%    computation.
%
% Output
%  Mu : N x (Number_of_Rings + 1) array giving the mean intensity value for
%    each ring. N is the number of objects. The object number corresponds
%    to the numbering by bwlabel/bwconncomp/... (The +1 is explained
%    below.)
%  Sigma : N x (Number_of_Rings + 1) array giving the standard deviation of
%    the intensity values for each ring. N is the number of objects.
%
% Note 1. This function computes Number_of_Rings rings with width
% Ring_Width starting from the boundary of the objects. The mean and
% standard devation of any region in the center of the object that is not
% included in a ring is also returned; this is why there are
% (Number_of_Rings + 1) columns in the output.
%
% Note 2. The last column of Mu and Sigma corresponds the the outermost
% ring. The first column corresponds to the central region of the object.
%
% Note 3. If an object is small, then the inner rings and central region
% might be empty. Any empty inner ring or central region will have a mean
% and standard deviation of 0 returned.
%
% Method: This function computes the pixels belonging to each ring by
% thresholding the distance transform of the objects. 

% James Kapaldo

N_R = 1 + options.Number_of_Rings; % number of rings; +1 will hold any leftover region at the center of the object
W = options.Ring_Width; % width of the N_R rings
useGPU = options.Use_GPU;


% Note: below I clear all variables that are no longer needed. This can be
% helpful when working on a GPU.
I = double(I(BW));

% Compute distance transform. This is used for computing the pixels in each
% ring. 
%    NOTE : There is some error with computing the dist. trans. on the GPU
%    with large images (1000x1000 works but 1500x1500 does not). Therefore,
%    compute it on the CPU, even though this is much slower. 
%    TODO : Try to get the distance transform to work on the GPU.
D_full = bwdist(~BW);
D = double(D_full(BW)); clear D_full;

if useGPU
    % If using gpu, then send over the arrays.
    I = gpuArray(I); 
    BW = gpuArray(BW);
    D = gpuArray(D);
end



% Compute the label matrix to identify each ring with a specific object
L_full = double(bwlabel(BW));
L = L_full(BW); clear L_full BW;

N_obj = max(L); % The number of objects
N_D = numel(D); % The number of object pixels 

% Create a temporary array that holds the logical indices for layer
if useGPU
    layer_idx = [gpuArray.false(N_D,1), D > ( N_R-1:-1:1 ) * W, gpuArray.true(N_D,1)]; %clear D
else
    layer_idx = [false(N_D,1), D > ( N_R-1:-1:1 ) * W, true(N_D,1)]; %clear D
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

% Compute the indices for the final output (first column is the object
% number, second column is the ring number. (Note that the first ring is on
% the inside of the object and the last ring is on the outside of the
% object.)
inds = [ L(ring_linidx), ring_num(ring_idx)]; clear ring_idx ring_num L;

% Compute the values of each pixel. Note vals should be the same size as I,
% just re-ordered. (The same goes for L above.)
vals = I(ring_linidx); clear ring_linidx I;

% Now, compute the mean and standard deviation. Note that the mean and
% standard deviation functions are not available for accumarray on a GPU;
% thus, I calculate them the long way.

% Compute the number of pixels in each ring of each object. The default
% value is 1, to prevent the NaN when computing the mean.

if useGPU
    N = accumarray(inds, gpuArray.ones(numel(vals),1,'double'), [N_obj, N_R], @sum, double(1));
else
%     N = accumarray(inds, ones(numel(vals),1,'double'), [N_obj, N_R], @sum, double(1));
    N = double(full(sparse(double(inds(:,1)), double(inds(:,2)), ones(numel(vals),1), double(N_obj), N_R)));
    N(N==0) = 1;
    
    % For some reason using sparse instead of accumarray here is faster for
    % computing N, but it is slower for computing the other quantities that
    % use accumarray below.
end

% Compute the sum of the intensities in each ring of each object
S = accumarray(inds, vals, [N_obj, N_R], @sum);

% Compute the mean intensity in each ring of each object
Mu = S ./ N;

% If the standard deviation was requested, then compute it
if nargout > 1
    
    % Compute the square difference between the intensity and mean
    % intensity
    Var_i = (vals - Mu(inds(:,1) + (inds(:,2)-1)*N_obj)).^2; clear vals
    
    % Compute the variance in each ring for each object
    Var = accumarray(inds, Var_i, [N_obj, N_R], @sum); clear Var_i inds
    
    % Compute the standard deviation
    Sigma = sqrt(Var./(N-1)); clear Var
    
    % Replace any Inf's with 0.
    Sigma(N==1) = 0;
    
end

if useGPU
    % If using gpu, then gather the results
    Mu = gather(Mu);
%     N = gather(N);
%     D = gather(D);
    if nargout > 1
        Sigma = gather(Sigma);
    end
end

end