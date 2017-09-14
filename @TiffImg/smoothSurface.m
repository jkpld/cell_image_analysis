function Z = smoothSurface(obj,z,type)

% If smoothing is not turned on, then output the input
if isnan(obj.Surface_Smoothing_Radius)
    Z = z;
    return;
end

if nargin<3
    type = 'LowessBisquare';
end

[Y,X] = ndgrid(obj.yCenters,obj.xCenters);
smoothRadius = (obj.Surface_Smoothing_Radius / (obj.blockSize*obj.mmPerPixel) )^2;

switch type

    case 'Lowess'
        
        Z = curvefit.LowessFit;
        Z.Span = pi*smoothRadius^2 / numel(z);
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));
    case 'LowessLAR'
        
        Z = curvefit.LowessFit;
        Z.Robust = 'LAR';
        Z.Span = pi*smoothRadius^2 / numel(z);
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));
    case 'LowessBisquare'
        
        Z = curvefit.LowessFit;
        Z.Robust = 'Bisquare';
        Z.Span = pi*smoothRadius^2 / numel(z);
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));
    case 'median'
        r = round(smoothRadius);
        D = strel('disk',r);
        D = D.Neighborhood;
        r = (size(D,1)-1)/2 + 1;
        order = (r-1)*size(D,1) + r;
        Z = ordfilt2(z,order,D,'symmetric');
        
        Z = smoothSurface(obj,Z,smoothRadius,'Lowess');
        
    otherwise
        
        Z = curvefit.LowessFit;
        Z.Span = pi*smoothRadius^2 / numel(z);
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));

end