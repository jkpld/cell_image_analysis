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
else
    bin = 1.5;
    smooth = 2.5*bin;
end

if nargin < 6
    name = 'unknown';
else
    name = char(name);
end

DEBUG = 0;


% Remove the stripe artifact. --------------------------------------------

% First flatten the intensity data.
options.smoothingRadius = 5*bin/tiffImg.mmPerPixel;
[fltnr, flattener] = TiffImg.decimate_and_smooth(x, y, I, options);
I_c = I ./ flattener(x,y);

if DEBUG
    figure %#ok<*UNRCH>
    tri = delaunay(fltnr.X,fltnr.Y);
    trisurf(tri,fltnr.X,fltnr.Y,fltnr.Z);
    line(x,y,I,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    title(['Initial data with flattening surface: channel ' name])
    view(0,0)
    setTheme(gcf,'dark')
    axis tight
    zlim([0,100])
    
    figure
    line(x,y,I_c,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    title(['Flattened data before x stripe correction : channel ' name])
    view(0,0)
    setTheme(gcf,'dark')
end

% Now compute the stripe with two iterations

[Xstripe1, stripe_x, ~, fsample] = Fit_Stripe_Artifact(x*tiffImg.mmPerPixel,I_c,'Threshold',1e-2,'generatePlots',DEBUG);

if isempty(fsample)
    % No stripe was detected
    Xstripe = [];
    I_c = I;
else
    stripe_x = stripe_x/tiffImg.mmPerPixel;

    size(Xstripe1)
    Xstripe1 = Xstripe1(:);
    stripe_x = stripe_x(:);
    if isnan(tiffImg.stripeWidth)
        maxPeriod = 3;
    else
        maxPeriod = 3*tiffImg.stripeWidth*tiffImg.mmPerPixel;
    end
    Xstripe1 = highpass(Xstripe1,maxPeriod,1/fsample) + 1;

    I_c = I ./ nakeinterp1(stripe_x, Xstripe1, x);
    
    xg = (1:tiffImg.imageSize(2)).';
    Xstripe = nakeinterp1(stripe_x, Xstripe1, xg);
    Xstripe = Xstripe';
end

% 
% % idx = I_c > 0.6 & I_c < 1.4; % Get data around flattened region
% idx = I_c > 0 & I_c < 4; % Get data around flattened region
% Xstripe1 = decimateData(x(idx),ones(size(x(idx))),I_c(idx),'binSize',[200,100],'defaultValue',1,'reductionMethod','mad'); % Get median dapi value from bins 100 pixels wide along x direction
% % Xstripe1 = decimateData(x(idx),ones(sum(idx),1),I_c(idx),'binSize',[100,100],'defaultValue',1,'reductionMethod','std'); % Get median dapi value from bins 100 pixels wide along x direction
% % Xstripe1.Z = highpass(Xstripe1.Z(:,2),3,100*tiffImg.mmPerPixel) + 1;
% I_c = I_c ./ nakeinterp1(Xstripe1.X(:,1), Xstripe1.Z(:,1), x); % Divide out the median value
% figure
% plot(Xstripe1.X(:,2),Xstripe1.Z(:,2))
% idx = I_c > 0.7 & I_c < 1.3; % Get data around flattened region
% % Xstripe2 = decimateData(x(idx),ones(sum(idx),1),I_c(idx),'binSize',[100,100],'defaultValue',1); % Get median dapi value from bins 100 pixels wide along x direction
% % Xstripe2.Z = highpass(Xstripe2.Z(:,2),3,100*tiffImg.mmPerPixel) + 1;
% 
% % Combine stripe corrections
% xg = (1:tiffImg.imageSize(2)).';
% Xstripe = nakeinterp1(Xstripe1.X(:,1), Xstripe1.Z(:,1), xg) ...
%     .* nakeinterp1(Xstripe2.X(:,1), Xstripe2.Z(:,1), xg);
% 
% Xstripe = Xstripe'; %should be row.
% 
% I_c = I ./ nakeinterp1(xg, Xstripe, x);
% 

% Compute flattening surface. --------------------------------------------
% Nuclei density
p2to4 = prctile(I_c,[2,4]);

rho = sum(I_c>=p2to4(1) & I_c<=p2to4(2))/prod(tiffImg.imageSize*tiffImg.mmPerPixel); % nuclei/mm^2

nucleiPerBin = 200; % emperical
bin = sqrt( (2/sqrt(3)) * nucleiPerBin / rho );
if bin > 3
    bin = 3;
    nucleiPerBin = bin^2 * rho / (2/sqrt(3));
end

smooth = 2.5*bin;
    
options.binSize = [bin/tiffImg.mmPerPixel,1]; % [2.5/tiffImg.mmPerPixel,1];
options.smoothingRadius = smooth/tiffImg.mmPerPixel; % 6/tiffImg.mmPerPixel;
options.defaultValue = nan;

% Decimate data using the mean of the data in the 2-4% range
    function v = reductionMethod(x)
        prct = prctile(x(x>0),[2,4]);
        idx = x>=prct(1) & x<=prct(2);
        
        v = mean(x(idx));
%         v = median(x(idx));
%         [numel(idx), sum(idx), v, mean(x(idx)), median(x(idx))]
        if sum(idx) < 0.5*nucleiPerBin
            v = nan;
        end
        
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
    tri = delaunay(S.X,S.Y);
    ts = trisurf(tri,S.X,S.Y,S.Z);
    line(x,y,I_c,'Marker','.','MarkerSize',1,'LineStyle','none','Color','g')
%     surface(flatteningSurface.x,flatteningSurface.y,flatteningSurface.Z)
    title(['Channel staining correction (2-4% value, after stripe crrctn) : channel ' name])
    setTheme(gcf,'dark')
    shading interp
    ts.EdgeColor = 'k';
    axis tight
    zlim([0.9*min(S.Z(:)),1.1*max(S.Z(:))])

end

end