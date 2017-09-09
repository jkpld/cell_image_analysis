function [Z, Z_shapeReproduction] = fourierDescriptors(B,N1,N2,useParallel)
% FOURIERDESCRIPTORS Compute fourier descriptors for a set of boundaries
%
% [Z, Z_shapeReproduction] = fourierDescriptors(B,N1,N2,useParallel)
%
% Input 
%   B : cell array of object boundaries (as returned by bwboundaries).
%   N1 : The number of elements that all boundaries will first be resized
%     to.
%   N2 : The number of fourier descripters returned.
%   useParallel : Logical flag. If true, then the boundaries will be
%     resized in parallel.
%
% Output
%   Z : N x (N2-2) array of translation, rotation, and scale invariant
%     fourier descriptors. Use these for as features for any machine
%     learning. N is the number of boundary elements. There are N2-2
%     columns because the first and last columns of these translation and
%     scale invariant features is always 0 and 1 respectively.
%   Z_shapeReproduction : N x N2 array of fourier descriptors that maintain
%     all translation, rotation, and scale information. Use these if you
%     want to recover the original objects exact position.
%
% - output variables are class Single
%
% This function will first resize each boundary to have N1 elements (so
% that the frequency components returned for each boundary are the same).
% The fourier transform will be computed Row i of matrix Z will correspond
% to boundary B{i}
%
% Note: N2 should be even.
% Note: The boundaries must have clockwise orientation.

% James Kapaldo

s = zeros(N1,numel(B));
if useParallel
    parfor i = 1:numel(B)
        m = size(B{i},1);
        xG = (1:m)';
        xi = linspace(1,m,N1)';
        s(:,i) = nakeinterp1(xG,B{i}(:,1),xi) + 1i * nakeinterp1(xG,B{i}(:,2),xi);
    end
else
    for i = 1:numel(B)
        m = size(B{i},1);
        xG = (1:m)';
        xi = linspace(1,m,N1)';
        s(:,i) = nakeinterp1(xG,B{i}(:,1),xi) + 1i * nakeinterp1(xG,B{i}(:,2),xi);
    end
end

s = single(s);

Z = fft(s);

% Extract only the number of components requested
d = round(N2/2);
Z((d+1):(N1-d),:) = [];

if nargout > 1
    % If the descriptors capbable of reproducing the shape are asked for,
    % then create these
    % Transpose so that each row is an obect and each column a feature.
    Z_shapeReproduction = Z.';
end

% Create the scale, rotation, translation invariant descriptors for
% learning.
Z = abs(Z); % rotation invariance
Z(1,:) = 0; % translation invariance
Z = Z ./ Z(end,:); % scale invariance

% remove constant quantities, which provide zero information about the
% shape
Z([1,end],:) = []; 

% Transpose so that each row is an obect and each column a feature.
Z = Z.';



