function [out] = GLCMFeatures(glcm)
% GLCMFEATURESJK Compute the 13 Haralick features for each page of the glcm
% matrix. (These are the 13 texture features used by CellProfiler.)
%
% Features computed 
% contrast [1]
% correlation [1]
% differenceEntropy [1]
% differenceVariance [1]
% energy [1]
% entropy [1]
% informationMeasureOfCorrelation1 [1]
% informationMeasureOfCorrelation2 [1]
% inverseDifferenceMoment [1]
% sumAverage [1]
% sumEntropy [1]
% sumOfSquaresVariance [1]
% sumVariance [1]
%
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973
%
%
% File started from
% https://www.mathworks.com/matlabcentral/fileexchange/55034-glcmfeatures-glcm-
%
% Modified by James Kapaldo to only output the 13 given by CellProfiler and
% to speed it up without the forloop over pages.


if nargin == 0
    return
elseif (nargin > 1) 
    error('Too many or too few input arguments')
else
    if ((size(glcm,1) <= 1) || (size(glcm,2) <= 1))
        error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcm,1) ~= size(glcm,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end  
end



% Get size of GLCM
nGrayLevels = size(glcm,1);
nglcm = size(glcm,3);

% Normalize the GLCMs
glcm = glcm./sumUp(glcm);


% Create indices for vectorising code:
I = (1:nGrayLevels).*ones(nGrayLevels,1);
J = I';

tmpOffsetInd = permute((0:nglcm-1)*nGrayLevels,[1,3,2]);

uX = sumUp( I.*glcm );
uY = sumUp( J.*glcm );
sX = sumUp( (I-uX).^2.*glcm );
sY = sumUp( (J-uY).^2.*glcm );

contrast                =  sumUp( glcm .* (I-J).^2           );
correlation             = (sumUp(I.*J.*glcm) - uX.*uY)./sqrt(sX.*sY);
energy                  =  sumUp( glcm .^ 2                  );
entropy                 = -sumUp( glcm .* log(glcm)          );
inverseDifferenceMoment =  sumUp( glcm ./ ( 1 + (I-J).^2)    );
sumOfSquaresVariance    =  sumUp( glcm .* (I-uX).^2          );

tmp1 = (I+J) + permute((0:nglcm-1)*(2*nGrayLevels-1),[1,3,2]) - 1;
tmp2 = abs(I-J) + tmpOffsetInd + 1;

pXplusY  = reshape(accumarray(tmp1(:),glcm(:)), [2*nGrayLevels-1, nglcm]); 
pXminusY = reshape(accumarray(tmp2(:),glcm(:)), [nGrayLevels, nglcm]);

idx1 = (2:2*nGrayLevels)';
idx2 = (0:nGrayLevels-1)';

sumAverage          =  sum( idx1 .* pXplusY                                                      , 1, 'omitnan');
sumEntropy          = -sum( pXplusY .* log(pXplusY)                                              , 1, 'omitnan');
differenceEntropy   = -sum( pXminusY .* log(pXminusY)                                            , 1, 'omitnan');
differenceVariance  =  sum( (idx2 - permute(sumUp( abs(I - J) .* glcm ),[1,3,2])).^2 .* pXminusY , 1, 'omitnan');
sumVariance         =  sum( (idx1 - sumAverage).^2 .* pXplusY                                    , 1, 'omitnan');

pX = sum(glcm,2);
pY = permute(sum(glcm,1), [2,1,3]);

tmpI = I + tmpOffsetInd;
tmpJ = J + tmpOffsetInd;
tmp = reshape(pX(tmpI(:)).*pY(tmpJ(:)), [nGrayLevels,nGrayLevels,nglcm]);

HXY1 = -sumUp( glcm .* log(tmp)   );
HXY2 = -sumUp( tmp  .* log(tmp)   );
HX =   -sum(   pX   .* log(pX)    , 1, 'omitnan');
HY =   -sum(   pY   .* log(pY)    , 1, 'omitnan');

informationMeasureOfCorrelation1 = (entropy - HXY1) ./ max(HX,HY);
informationMeasureOfCorrelation2 = sqrt(1 - exp(-2.*(HXY2 - entropy)));

out = [contrast(:),...
correlation(:),...
differenceEntropy(:),...
differenceVariance(:),...
energy(:),...
entropy(:),...
informationMeasureOfCorrelation1(:),...
informationMeasureOfCorrelation2(:),...
inverseDifferenceMoment(:),...
sumAverage(:),...
sumEntropy(:),...
sumOfSquaresVariance(:),...
sumVariance(:)];


end

function x = sumUp(X)
    x = sum(sum( X, 1,'omitnan'), 2,'omitnan');
end
