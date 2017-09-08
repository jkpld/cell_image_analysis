function x = intensity_features(L, I, options)
% INTENSITY_FEATRUES Return the integrated, mean,
% standard deviation, skewness, and kurtosis of the intensities for each
% object in an image I labeled with matrix L.
%
% x = intensity_features(I, L, options)
%
% Input 
%  L : object label matrix
%  I : image
%  options : stucture with the following fields
%    Use_GPU : logical flag. If true, then a GPU will be used to speed up
%    the computation
%
% Output
%  x : N x 5 array where N is the number of objects. The columns are --
%    (integral | mean | standard deviation | skewness | kurtosis )

% James Kapaldo

useGPU = options.Use_GPU;

FG = L~=0;

I = single(I(FG));
L = single(L(FG));

if useGPU
    I = gpuArray(I);
    L = gpuArray(L);
end

N_obj = max(L);

% NOTE : accumarray on gpu does not accept @mean, @std, or the others;
% thus, they need to be calculated out by hand. This is actually faster
% than using the other commands anyway, because many of the commands
% calculate the same thing.

% Pixels per object
if useGPU
    N = accumarray(L, gpuArray.ones(size(L),'single'), [N_obj, 1], @sum, single(1));
else
    % Using sparse on CPU is faster than accumarray
    N = single(full(sparse(double(L),ones(size(L)), ones(size(L)), double(N_obj), 1)));
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