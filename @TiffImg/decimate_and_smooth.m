function [S,fun] = decimate_and_smooth(x,y,z,varargin)
% DECIMATE_AND_SMOOTH Decimate scatter point data, smooth the surface, and
% return a function handle that can be used to evaluate the surface.
%
% fun = decimate_and_smooth(x,y,z,options)
%
% Input
%   x : spatial x data
%   y : spatial y data
%   z : value data
%   options : parameter value pairs.
%     smoothingRadius : the radius (in the same units as the data) which
%       the surface will be smoothed over. If value is not finite, then the
%       surface will not be smoothed.
%     smoothType : The type of surface smoothing to perform {'Lowess',
%       'LowessLAR', 'LowessBisquare', 'median'}. Default is
%       'LowessBisquare'.
%
%     Any/All options accepted by decimateData() function.
%
% Output 
%   fun : function handle that takes in x and y and evaluts the surface.
%     z = fun(x,y);

% James Kapaldo

% Parse inputs
p = inputParser;
p.FunctionName = 'decimate_and_smooth';
p.KeepUnmatched = 1;

addParameter(p,'smoothingRadius',nan, @(t) validateattributes(t,{'numeric'},{'scalar'}))
addParameter(p,'smoothType','LowessBisquare', @(t) ischar(t))
parse(p,varargin{:})
sR = p.Results.smoothingRadius;
sType = p.Results.smoothType;
options = p.Unmatched; % options to pass through to decimateData

% offset the hexagons so that the first one starts in the corner of the
% data (not centered at 0 0).
if strcmp(options.gridType,'hexagonal')
    Xoffset = 0.9*options.binSize(1)/2;
    Yoffset = 0.9*prod(options.binSize)/2;
    
    x = x - Xoffset;
    y = y - Yoffset;
end

[dcmt, options] = decimateData(x,y,z,options);

if strcmp(options.gridType,'hexagonal')
    dcmt.X = dcmt.X + Xoffset;
    dcmt.Y = dcmt.Y + Yoffset;
end

if isfinite(sR)
    switch options.gridType
        case 'rectangular'
            % mean bin size
            gridArea = prod(options.binSize);
        case 'hexagonal'
            % approximate squewed hexagon with inscribed ellipse
            gridArea = pi * (options.binSize(1)/2)^2 * options.binSize(2);     
    end

    % Get the number of grid points in the smoothing radius
    sR = sqrt((pi*sR^2)/gridArea)/pi;

    [Z,X,Y] = TiffImg.smoothSurface(dcmt.X,dcmt.Y,dcmt.Z,sR,sType);
    
else
    Z = dcmt.Z;
    X = dcmt.X;
    Y = dcmt.Y;
end

S.X = X;
S.Y = Y;
S.Z = Z;

if nargout > 1
    if ~any(size(Z)==1) % assume the data is gridded
        f = griddedInterpolant(Y',X',Z);
        fun = @(x,y) f(y,x);
    else
        fun = scatteredInterpolant(X,Y,Z);
    end
end

end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
