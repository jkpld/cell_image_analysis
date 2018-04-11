function B = invFourierDescriptors(z,N1,N2)
% INVFOURIERDESCRIPTORS Compute the original shape from the fourier
% descriptors.
%
% B = invFourierDescriptors(z,N1,N2)
%
% Compute the the shape using N2 fourier descriptors. B will be returned as
% a complex matrix. Each column will represent a boundary. The real part is
% x, the complex part is y.


% Undo the transpose at the end of fourierDescriptors()
z = z.';

if size(z,1) == N2 - 2
    % The features given were for machine learning and are translation and
    % scale invarient. Add back in the columns.
    z = [zeros(1,size(z,2)); z; ones(1,size(z,2))];
end

d = round(N2/2);

z = [z(1:d,:); zeros(N1-N2,size(z,2)); z(d+1:end,:)];

B = ifft(z);

%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
