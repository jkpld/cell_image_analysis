function fun = interpolator2d(varargin)

charInput = cellfun(@ischar, varargin);
numCharInput = sum(charInput);

createGrid = contains(varargin(charInput), 'expand');
interpolateInputPoints = true;

if nargin - numCharInput == 1
    Z = varargin{1};
    interpolateInputPoints = false;
    
elseif nargin - numCharInput == 3
    x = varargin{1};
    y = varargin{2};
    Z = varargin{3};
    if numel(x) ~= size(Z,2) || numel(y) ~= size(Z,1)
        error('interpolator2d:badVectorSizes', 'Vector x and y should have size(Z,2) and size(Z,1) number of elements, respectively.')
    end
    x = x(:);
    y = y(:);
else
    error('interpolator2d:badInput','Function input can be matrix Z; or, it can be vector x, vector y, matrix Z.')
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