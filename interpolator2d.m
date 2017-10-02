function fun = interpolator2d(varargin)

%
% fun = interpolator2d(Z)
% fun = interpolator2d(Z, expand)
% fun = interpolator2d(x, y, Z)
% fun = interpolator2d(x, y, Z, expand)
%
% Z : 2d matrix size (n x m)
% x : vector, numel(x) = m
% y : vector, numel(y) = n
% expand : logical. if true, then fun(xi,yi) will be matrix (size(yi,1) x
%   size(xi,2)), otherwise will be array (numel(xi) x 1). (default value is
%   true)
%
% output class is same as Z



switch nargin
    case {1,2}
        Z = varargin{1};
        interpolateInputPoints = false;
    case {3,4}
        x = varargin{1};
        y = varargin{2};
        Z = varargin{3};
        if numel(x) ~= size(Z,2) || numel(y) ~= size(Z,1)
            error('interpolator2d:badVectorSizes', 'Vector x and y should have size(Z,2) and size(Z,1) number of elements, respectively.')
        end
        x = x(:);
        y = y(:);
        interpolateInputPoints = true;
end
createGrid = true;
if nargin == 2 || nargin == 4
    createGrid = logical(varargin{end}(1));
end


if interpolateInputPoints    
    xg = (1:size(Z,2))';
    yg = (1:size(Z,1))';
    
    xi = @(xi) reshape(nakeinterp1(x, xg, xi(:)), size(xi));
    yi = @(yi) reshape(nakeinterp1(y, yg, yi(:)), size(yi));
else
    xi = @(xi) xi;
    yi = @(yi) yi;
end

if createGrid
    fun = @(x,y) cast(interp2mex_wExpand(double(Z), xi(double(x)), yi(double(y))),'like',Z);
else
    fun = @(x,y) cast(interp2mex(double(Z), xi(double(x(:))), yi(double(y(:)))),'like',Z);
end


end