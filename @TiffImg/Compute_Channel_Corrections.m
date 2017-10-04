function [flatteningSurface, Xstripe] = Compute_Channel_Corrections(tiffImg, x, y, I, options, name)
% COMPUTE_CHANNEL_CORRECTIONS Compute a smooth surface that hopefully
% compensates for uneven staining, and compute the Xstripe artifact.
%
% [flatteningSurface, Xstripe] = Compute_Channel_Corrections(tiffImg, x, y, I, options)
%
% Input 
%   x : spatial x data [pixels]
%   y : spatial y data [pixels]
%   I : integrated intensity for each object
%   options : (optional) structure with options for decimate_smooth()
%
% Output
%   flatteningSurface : struct with fields x, y, and Z. The smooth surface
%     that ties to compensate for uneven staining (see more below)
%   Xstripe : The X stripe artifact (1 x imageSize(2)) that, when the
%     channel data is divided by, removes the x stripe artifact.
%
% Note : The flatteningSurface is computed by using the 2% to 4% intensity
% data with large bin sizes.
%
% Note : It is suggested that input data be only from G1 nuclei. This could
% lead to a better correction.

% James Kapaldo

% Nuclei density
rho = numel(x)/prod(tiffImg.imageSize*tiffImg.mmPerPixel); % nuclei/mm^2

if nargin < 5 || isempty(options)
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

if nargin < 6
    name = 'unknown';
else
    name = char(name);
end

DEBUG = 1;

if DEBUG
    figure %#ok<*UNRCH>
    line(x,y,I,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    title(['Initial data : channel ' name])
    view(0,0)
    setTheme(gcf,'dark')
end

% Remove the stripe artifact. --------------------------------------------

% First flatten the intensity data.
[~, flattener] = TiffImg.decimate_and_smooth(x, y, I, options);
I_c = I ./ flattener(x,y);

if DEBUG
    figure
    line(x,y,I_c,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    title(['Flattened data before x stripe correction : channel ' name])
    view(0,0)
    setTheme(gcf,'dark')
end

% Now compute the stripe with two iterations
idx = I_c > 0.6 & I_c < 1.4; % Get data around flattened region
Xstripe1 = decimateData(x(idx),ones(sum(idx),1),I_c(idx),'binSize',[100,100],'defaultValue',1); % Get median dapi value from bins 100 pixels wide along x direction
Xstripe1.Z = highpass(Xstripe1.Z(:,2),3,100*tiffImg.mmPerPixel) + 1;
I_c = I_c ./ nakeinterp1(Xstripe1.X(:,1), Xstripe1.Z(:,1), x); % Divide out the median value

idx = I_c > 0.7 & I_c < 1.3; % Get data around flattened region
Xstripe2 = decimateData(x(idx),ones(sum(idx),1),I_c(idx),'binSize',[100,100],'defaultValue',1); % Get median dapi value from bins 100 pixels wide along x direction
Xstripe2.Z = highpass(Xstripe2.Z(:,2),3,100*tiffImg.mmPerPixel) + 1;

% Combine stripe corrections
xg = (1:tiffImg.imageSize(2)).';
Xstripe = nakeinterp1(Xstripe1.X(:,1), Xstripe1.Z(:,1), xg) ...
    .* nakeinterp1(Xstripe2.X(:,1), Xstripe2.Z(:,1), xg);

Xstripe = Xstripe'; %should be row.

I_c = I ./ nakeinterp1(xg, Xstripe, x);

if DEBUG    
    figure
    line(x,y,I_c,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    title(['After x stripe correction (without flattening) : channel ' name])
    view(0,0)
    setTheme(gcf,'dark')
    
    figure
    line(xg*tiffImg.mmPerPixel,Xstripe,'color','y')
    title(['XStripe : channel ' name])
    setTheme(gcf,'dark')
end

% Compute flattening surface. --------------------------------------------
nucleiPerBin = 800; % emperical
bin = sqrt( (2/sqrt(3)) * nucleiPerBin / rho );
smooth = 1.5*bin;
    
options.binSize = [bin/tiffImg.mmPerPixel,1]; % [2.5/tiffImg.mmPerPixel,1];
options.smoothingRadius = smooth/tiffImg.mmPerPixel; % 6/tiffImg.mmPerPixel;
options.defaultValue = nan;

% Decimate data using the mean of the data in the 2-4% range
    function v = reductionMethod(x)
        prct = prctile(x,[2,4]);
        v = mean(x(x>=prct(1) & x<=prct(2)));
    end
options.reductionMethod = @reductionMethod;

S = TiffImg.decimate_and_smooth(x, y, I_c, options);

% This surface would need to be evaluated using scattered interpolants;
% however, this is quite slow to evaluate, so reinterplate onto the same
% square grid used to define the background surface.
fun = scatteredInterpolant(double(S.X),double(S.Y),double(S.Z));
xg = tiffImg.BG_offset.x;
yg = tiffImg.BG_offset.y;
Z = fun({xg,yg})';

flatteningSurface.Z = Z;
flatteningSurface.x = xg;
flatteningSurface.y = yg;

if DEBUG
    figure
    line(x,y,I_c,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    surface(flatteningSurface.x,flatteningSurface.y,flatteningSurface.Z)
    title(['Channel staining correction (2-4% value, after stripe crrctn) : channel ' name])
    setTheme(gcf,'dark')
end

end