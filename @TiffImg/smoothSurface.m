function [Z,X,Y] = smoothSurface(x,y,z,smoothRadius,type)

% If the number of elements in x is not the same as the number of elements
% in Z, then assume we need to create a grid.
if numel(x) == numel(z)
    Y = y;
    X = x;
else
    [Y,X] = ndgrid(y,x);
end

% If smoothing is not turned on, then output the input
if ~isfinite(smoothRadius)
    Z = z;
    return;
end

% If less than 5 points, cannot smooth;
if numel(z) < 5
    Z = z;
    return;
end

% If the surface is an integer, then convert it to single.
if isinteger(z)
    z_class = class(z);
    z = single(z);
end

if nargin<5
    type = 'LowessBisquare';
end

span = max(min(1,pi*smoothRadius^2 / numel(z)),0);

switch type

    case 'Lowess'
        
        Z = curvefit.LowessFit;
        Z.Span = span;
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));
        
    case 'LowessLAR'
        
        Z = curvefit.LowessFit;
        Z.Robust = 'LAR';
        Z.Span = span;
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));
        
    case 'LowessBisquare'
        
        Z = curvefit.LowessFit;
        Z.Robust = 'Bisquare';
        Z.Span = span;
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));
        
    case 'median'
        if any(size(z)==1)
            error('smoothSurface:medianUnsupported','Median filtering is only supported for gridded data. If your data is gridded, then consider reshaping it to a matrix before passing to this function.')
        end
        
        r = round(smoothRadius);
        D = strel('disk',r);
        D = D.Neighborhood;
        r = (size(D,1)-1)/2 + 1;
        order = (r-1)*size(D,1) + r;
        Z = ordfilt2(z,order,D,'symmetric');        
        
    otherwise
        
        Z = curvefit.LowessFit;
        Z.Span = span;
        Z = fit(Z, [Y(:),X(:)], z(:));
        Z = reshape(evaluate(Z, [Y(:),X(:)]),size(X));

end

if isinteger(z)
    Z = cast(Z,z_class);
end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
