function [flatteningSurface, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(tiffImg,x,y,dapi,area,options)
% COMPUTE_DAPI_CORRECTIONS Compute the smooth surface and X stripe artifact
% such that when the DAPI data is divided by both of these corrections the
% G1 and G2 bands are both flat and located at 1 and 2, respectively, and
% there is not stripe artifact along the x direction.
%
% [flatteningSurface, Xstripe, G1_areaNorm, G1_idx] = Compute_DAPI_Corrections(x,y,dapi,mmPerPixel)
%
% Input 
%   x : spatial x data [pixels]
%   y : spatial y data [pixels]
%   dapi : integrated dapi intensity for each object
%   area : area of each object
%   options : (optional) structure with options for decimate_smooth()
%
% Output
%   flatteningSurface : struct with fields x, y, and Z. The smooth surface
%     that, when the dapi data is divided by, positions the G1 and G2 bands
%     in the correct positions.
%   Xstripe : The X stripe artifact (1 x imageSize(2)) that, when the dapi
%     data is divided by, removes the x stripe artifact from the dapi data.
%   G1Area : struct with fields x, y, and Z. The smooth surface that gives
%     the G1 area.
%   G1_idx : Logical index array giving the G1 nuclei of the input data.

% James Kapaldo

% Nuclei density
rho = numel(x)/prod(tiffImg.imageSize*tiffImg.mmPerPixel); % nuclei/mm^2

if nargin < 6
    nucleiPerBin = 300; % emperical
    bin = sqrt( (2/sqrt(3)) * nucleiPerBin / rho );
    smooth = 2.5*bin;
    options = struct( ...
        'reductionMethod', 'mode', ...
        'gridType', 'hexagonal', ...
        'cleanHexagonData', true, ...
        'binSize', [bin/tiffImg.mmPerPixel,1], ... % 1.5/mmPerPixel
        'smoothingRadius',smooth/tiffImg.mmPerPixel, ... % 4/mmPerPixel
        'defaultValue',1);
end

DEBUG = 1;

% Initial correction. ----------------------------------------------------

% Flatten the g1 band and position it approximately at 1 (assuming the G1
% band is the mode)

[G1_1,DAPI_mode] = TiffImg.decimate_and_smooth(x, y, dapi, options);
dapi_c = dapi ./ DAPI_mode(x,y);

if DEBUG
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('after initial g1 flattening')
    setTheme(gcf,'dark')
    axis tight
    zlim([0,3])
    view(0,0)
end

% Remove the stripe artifact with two iterations of fitting. -------------

% Select the g1 band. Compute the median dapi value along small x-slices to
% extract the stripe. Divide the stripe away.

idx = dapi_c > 0.6 & dapi_c < 1.4; % Get G1 band
G1_stripe1 = decimateData(x(idx),ones(sum(idx),1),dapi_c(idx),'binSize',[100,100],'defaultValue',1); % Get median dapi value from bins 100 pixels wide along x direction
dapi_c = dapi_c ./ nakeinterp1(G1_stripe1.X(:,1), G1_stripe1.Z(:,1), x); % Divide out the median value

idx = dapi_c > 0.7 & dapi_c < 1.3; % Get G1 band
G1_stripe2 = decimateData(x(idx),ones(sum(idx),1),dapi_c(idx),'binSize',[100,100],'defaultValue',1);  % Get median dapi value from bins 100 pixels wide along x direction
dapi_c = dapi_c ./ nakeinterp1(G1_stripe2.X(:,1), G1_stripe2.Z(:,1), x); % Divide out the median value

if DEBUG
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('after stripe correction')
    setTheme(gcf,'dark')
    axis tight
    zlim([0,3])
    view(0,0)
end

% Flatten the G1 band, again. --------------------------------------------

% Select the g1 band. Fit the mode with a smooth surface. Divide the
% surface away.

idx = dapi_c > 0.7 & dapi_c < 1.3;
[G1_2, DAPI_G1_mode] = TiffImg.decimate_and_smooth(x(idx), y(idx), dapi_c(idx), options);
dapi_c = dapi_c ./ DAPI_G1_mode(x,y);

% Flatten the G2 band. ---------------------------------------------------

% Select the g2 band. Fit the mode with a smooth surface. Divide the
% surface away.

options.reductionMethod = 'median'; % Use median instead of mode as the distribution is not that narrow.
options.defaultValue = 2;

idx = dapi_c > 1.7 & dapi_c < 2.5;
[G2_1, DAPI_G2_mode] = TiffImg.decimate_and_smooth(x(idx), y(idx), dapi_c(idx), options);
dapi_c = dapi_c ./ (DAPI_G2_mode(x,y)/2); % dapi_c = ((dapi_c-1) ./ (DAPI_G2_mode(x,y)/2)) + 1;

% Remove any stripe artifact from the G2 band. ---------------------------

% Select the g2 band. Compute the median dapi value along small x-slices to
% extract the stripe. Divide the stripe away.

idx = dapi_c > 1.7 & dapi_c < 2.3;
G2_stripe = decimateData(x(idx),ones(sum(idx),1),dapi_c(idx),'binSize',[100,100],'defaultValue',2);
dapi_c = dapi_c ./ (nakeinterp1(G2_stripe.X(:,1), G2_stripe.Z(:,1), x)/2) ;

idx = dapi_c > 0.7 & dapi_c < 1.3; % Get G1 band
G1_stripe3 = decimateData(x(idx),ones(sum(idx),1),dapi_c(idx),'binSize',[100,100],'defaultValue',1);  % Get median dapi value from bins 100 pixels wide along x direction
dapi_c = dapi_c ./ nakeinterp1(G1_stripe3.X(:,1), G1_stripe3.Z(:,1), x); % Divide out the median value


if DEBUG
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('Data after correction')
    setTheme(gcf,'dark')
    zlim([0,3])

    figure
    ax = axes;
    hold on
    h = histogram(ax,dapi_c,'BinLimits',[0,3]);
end

% Combine all results into a single smooth surface and stripe. -----------

% Interpolate the stripe to all for each x-pixel of the image.
Xstripe = nakeinterp1(G1_stripe1.X(:,1), G1_stripe1.Z(:,1), (1:tiffImg.imageSize(2)).') ...
    .* nakeinterp1(G1_stripe2.X(:,1), G1_stripe2.Z(:,1), (1:tiffImg.imageSize(2)).') ...
    .* nakeinterp1(G1_stripe3.X(:,1), G1_stripe3.Z(:,1), (1:tiffImg.imageSize(2)).') ...
    .* nakeinterp1(G2_stripe.X(:,1), G2_stripe.Z(:,1)/2, (1:tiffImg.imageSize(2)).');

Xstripe = Xstripe'; %should be row.

% Construct the flattening surface

% Each flattening surface could have a different number of points. If this
% is the case, then they need to be interpolated onto a grid of the same
% size before combining them.

% Grid to re-interpolate results on - use the grid for creating the
% background.
xg = tiffImg.BG_smooth.x;
yg = tiffImg.BG_smooth.y;

if isequal(G1_1.X,G1_2.X) && isequal(G1_2.X,G2_1.X)
    Z = G1_1.Z .* G1_2.Z .* G2_1.Z/2;
    X = G1_1.X;
    Y = G1_1.Y;
    
    % This surface would need to be evaluated using scattered interpolants;
    % however, this is quite slow to evaluate, so reinterplate onto the same
    % square grid
    fun = scatteredInterpolant(double(X),double(Y),double(Z));
    Z = fun({xg,yg})';
else
    fun = scatteredInterpolant(double(G1_1.X),double(G1_1.Y),double(G1_1.Z));
    
    G1_1Z = fun({xg,yg})';
    
    fun = scatteredInterpolant(double(G1_2.X),double(G1_2.Y),double(G1_2.Z));
    G1_2Z = fun({xg,yg})';
    
    fun = scatteredInterpolant(double(G2_1.X),double(G2_1.Y),double(G2_1.Z));
    G2_1Z = fun({xg,yg})';
    
    Z = G1_1Z .* G1_2Z .* G2_1Z/2;
end

flatteningSurface.Z = Z;
flatteningSurface.x = xg;
flatteningSurface.y = yg;

% Compute the G1 area surface --------------------------------------------
if nargout > 2
    nucleiPerBin = 800; % emperical
    bin = sqrt( (2/sqrt(3)) * nucleiPerBin / rho );
    smooth = 2.5*bin;
    idx = dapi_c > 0.7 & dapi_c < 1.3;
    options.defaultValue = nan;
    options.binSize = [bin/tiffImg.mmPerPixel,1];%[2.5/tiffImg.mmPerPixel,1];
    options.smoothingRadius = smooth/tiffImg.mmPerPixel;%6/tiffImg.mmPerPixel;
    options.reductionMethod = 'median';
    G1_area  = TiffImg.decimate_and_smooth(x(idx), y(idx), area(idx), options);
    
    % interpolate onto the same square grid as the threshold
    fun = scatteredInterpolant(double(G1_area.X), double(G1_area.Y), double(G1_area.Z));
    Z = fun({xg,yg})';
    
    figure
    
    % Plot surface
    tri = delaunay(G1_area.X,G1_area.Y);
    trisurf(tri,G1_area.X,G1_area.Y,G1_area.Z);
    line(x(idx),y(idx),area(idx),'marker','.','linestyle','none','color','g','markersize',1)
    
%     surface(xg,yg,Z)
    title('Area surface')
    setTheme(gcf,'dark')
    axis tight
    
    G1Area.Z = Z;
    G1Area.x = xg;
    G1Area.y = yg;
end

if nargout > 3
    G1_idx = idx;
end

if DEBUG
    tmp = @(x,y) interp2mex(flatteningSurface.Z, nakeinterp1(xg(:),(1:size(flatteningSurface.Z,2))',x), nakeinterp1(yg(:),(1:size(flatteningSurface.Z,1))',y));%griddedInterpolant({yg,xg},Z,'linear','nearest');%
    tmp_s = @(x) nakeinterp1((1:tiffImg.imageSize(2)).', Xstripe(:), x);
    
    dapi_c = dapi ./ tmp(x,y) ./ tmp_s(x);
    
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('Data when applying correction all at once')
    setTheme(gcf,'dark')
    axis tight
    zlim([0,3])
    
    histogram(ax,dapi_c,'BinEdges',h.BinEdges,'FaceColor','r')
end
end