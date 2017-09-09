function B = invFourierDescriptors(z,N1,N2)
% INVFOURIERDESCRIPTORS Compute the original shape from the fourier
% descriptors.
%
% B = invFourierDescriptors(z,N1,N2)
%
% Compute the the shape using N2 fourier descriptors. B will be returned as
% a complex matrix. Each column will represent a boundary. The real part is
% x, the complex part is y.
%
% N1 = size(z,1);

% Undo the transpose at the end of fourierDescriptors()
z = z.';
d = round(N2/2);

z = [z(1:d,:); zeros(N1-N2,size(z,2)); z(d+1:end,:)];

B = ifft(z);
