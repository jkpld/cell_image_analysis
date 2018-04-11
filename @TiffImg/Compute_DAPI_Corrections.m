function [FG_f, FG_o, Xstripe, G1Area, G1_idx] = Compute_DAPI_Corrections(tiffImg,x,y,dapi,area,options)
% COMPUTE_DAPI_CORRECTIONS Compute corrections that : will flatten the
% DAPI G1 band and position the G1 band at 1 and the G2 band at 2, remove
% the x-stripe artifact, and flatten/normalize the G1 area.
%
% [FG, BG, Xstripe, G1area, G1_idx] = Compute_DAPI_Corrections(tiffImg,x,y,dapi,area,options)
%
% Input 
%   x : spatial x data [pixels]
%   y : spatial y data [pixels]
%   dapi : integrated dapi intensity for each object
%   area : area of each object
%   options : (optional) structure with options for decimate_smooth()
%
% Output
%   FG_f : struct with fields x, y, and Z. This is a surface that the DAPI
%     channel should be divided by.
%   FG_o : struct with fields x, y, and Z. This is the surface that should
%     be added to the DAPI channel *after* dividing by FG.
%   Xstripe : The X stripe artifact (1 x imageSize(2)) that, when the dapi
%     data is divided by, removes the x stripe artifact from the dapi data.
%   G1Area : struct with fields x, y, and Z. The smooth surface that gives
%     the G1 area.
%   G1_idx : Logical index array giving the G1 nuclei of the input data.

% James Kapaldo

mmPP = tiffImg.mmPerPixel;

% Nuclei density
imageArea = prod(tiffImg.imageSize*mmPP);
rho = numel(x)/imageArea; % nuclei/mm^2

if nargin < 6
    nucleiPerBin = 300; % emperical
    bin = sqrt( (2/sqrt(3)) * nucleiPerBin / rho );
    smooth = 5*bin;%2.5*bin
    options = struct( ...
        'reductionMethod', 'mode', ...
        'gridType', 'hexagonal', ...
        'cleanHexagonData', true, ...
        'binSize', [bin/mmPP,1], ... % 1.5/mmPerPixel
        'smoothingRadius',smooth/mmPP, ... % 4/mmPerPixel
        'defaultValue',1);
else
    bin = options.binSize(1)*mmPP;
end

DEBUG = 0;

%% Initial correction. ----------------------------------------------------

% Flatten the g1 band with a large smoothing radius. This will be used to
% remove outliers.
options.smoothingRadius = 5*bin/mmPP;
[G1_0,DAPI_mode] = TiffImg.decimate_and_smooth(x, y, dapi, options);
dapi_c = dapi ./ DAPI_mode(x,y);

% Flatten the g1 band and position it approximately at 1 (assuming the G1
% band is the mode)
options.smoothingRadius = 2.5*bin/mmPP;
idx = dapi_c > 0.5 & dapi_c < 1.5;
[G1_1,DAPI_mode] = TiffImg.decimate_and_smooth(x(idx), y(idx), dapi(idx), options);
dapi_c = dapi ./ DAPI_mode(x,y);

if DEBUG
    figure
    hold on
    tri = delaunay(G1_1.X,G1_1.Y);
    trisurf(tri,G1_1.X,G1_1.Y,G1_1.Z);
    tri = delaunay(G1_0.X,G1_0.Y);
    trisurf(tri,G1_0.X,G1_0.Y,G1_0.Z);
    line(x,y,dapi,'marker','.','linestyle','none','color','g','markersize',1)
    title('original data with fitted surface')
    setTheme(gcf,'dark')
    drawnow
    
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('after initial g1 flattening')
    setTheme(gcf,'dark')
    axis tight
    zlim([0,3])
    view(0,0)
end

%% Remove the stripe artifact. ------------------------------------------
rho = numel(x)/(range(x)*mmPP);
% initBin = max(200/rho,0.005)
initBin = 0.005*(1 + 9*(200/rho>0.005));
[G1_stripe, G1_stripeX, ~, fsample] = Fit_Stripe_Artifact(x*mmPP,dapi_c,'Threshold',1e-2,'generatePlots',DEBUG,'initialBinSize',initBin);

if isempty(fsample)
    % No stripe was detected
    Xstripe = [];
else
    % Covert locations back to pixels
    G1_stripeX = G1_stripeX/mmPP;

    % Make columns
    G1_stripe = G1_stripe(:);
    G1_stripeX = G1_stripeX(:);
    
    % Get maximum period for highpass filter
    if isnan(tiffImg.stripeWidth)
        maxPeriod = 3;
    else
        maxPeriod = 3*tiffImg.stripeWidth*mmPP;
    end

    % Sent the computed stripe through a highpass filter and then recenter
    % at 1.
    G1_stripe = highpass(G1_stripe,maxPeriod,1/fsample) + 1;
    
    % Apply stripe correction
    dapi_c = dapi_c ./ nakeinterp1(G1_stripeX(:), G1_stripe(:), x);

    % Interpolate the stripe to all for each x-pixel of the image.
    Xstripe = nakeinterp1(G1_stripeX(:), G1_stripe(:), (1:tiffImg.imageSize(2)).');
    Xstripe = Xstripe'; %should be row.
end

%% Flatten the G1 band, again. --------------------------------------------

% Select the g1 band. Fit the mode with a smooth surface. Divide the
% surface away.
options.reductionMethod = 'median';

idx = dapi_c > 0.7 & dapi_c < 1.3;
[G1_2, DAPI_G1_mode] = TiffImg.decimate_and_smooth(x(idx), y(idx), dapi_c(idx), options);
dapi_c = dapi_c ./ DAPI_G1_mode(x,y);

%% Flatten the G2 band. ---------------------------------------------------
% Flattening the G2 band can be more difficult as it could be much less
% dence than the G1 band[1]. Thus, the fitting uses a more envolved
% reductionMethod() for the decimation.
%
% [1] : The cells I wrote this code for spend about 80% of there time, or
% much more, in G1 phase.

% Approximately select the G2 band
idx = dapi_c > 1.8 & dapi_c < 2.8;

% Get the bin size for use with the G2 band
nucleiPerBin = 300; % emperical
bin = sqrt( (2/sqrt(3)) * nucleiPerBin * imageArea / sum(idx) );
% fprintf('dapi bin size : %f\n', bin)

% Compute the histogram of the G2 band accross the entire image. The
% maximum value should give the approximate location of the G2 band. Note,
% that we can expect the G2 band to be mostly flat because it should
% ideally be exactly twice the G1 band, and since we have already flattened
% the G1 band, the G2 band should ideally be flat. In reality, the G2 band
% can have small changes accross the image for the same reason the G2 band
% is not located at 2 after correction the G1 band.
[nn,edg] = histcounts(dapi_c(idx),'BinWidth',0.05);
cnt = edg(1:end-1) + diff(edg)/2;
% flip the count and centers so that maximum returns the largest center if
% two bins have the same count.
nn = flip(nn);
cnt = flip(cnt);
[~,nI] = max(nn);
cnt = cnt(nI);
% This center will be used to create a weighting gaussian in the
% reductionMethod below.

% figure
% ax1 = axes;
% hold on

    function v = reductionMethod(x)
        % Idea : 
        % * Compute the histogram of the points.      
        % * Weight the counts using a gaussian with a FWHM of 1 centered on
        % the previously computed average G2 center. (This weighing helps
        % select the best value when several bins have almost the same
        % count accros the G2 band.)
        % * Do not accept the value if it is the smallest possible bin
        % center.
        % * Do not consider the image region if it only had a single bin
        % * Do not consider the image region if the maximum count value is
        % about the minimum requested density.
        
        % Bin the data of the image region
        [n,edgs] = histcounts(x,'BinWidth',0.05);
        cnts = edgs(1:end-1) + diff(edgs)/2;
        
        % Weight the data with the gaussian
        n = n.*exp(-(cnts-cnt).^2/(2*(1/2.355)^2)); % weight the counts towards the approximate value
        
        [n,sI] = sort(n,'descend');
        cnts = cnts(sI);
        
        % Do not use the result if there was only one bin
        if numel(n)==1
            v = nan;
            return;
        end
        % Do not consider the smallest value
        if sI(1)==1 
            n(1) = [];
            cnts(1) = [];
        end
        
        % If the distribution is close to uniform, then set the value to
        % the default value
        if n(1) < 1.2*nucleiPerBin/numel(n)
            v = nan;
            return;
        else
            v = cnts(1);
        end
        if isempty(v)
            v = nan;
        end

%         if v < 2.1
%             plot(cnts,n,'parent',ax1)
%         end
%         [v, n(1), 1.2*nucleiPerBin/numel(n)]
    end

% Set the options
smooth = 2.5*bin;
options.reductionMethod = @reductionMethod; % Use median instead of mode as the distribution is not that narrow.
options.defaultValue = nan;
options.binSize = [bin/mmPP,1]; % [2.5/mmPP,1];
options.smoothingRadius = smooth/mmPP; % 6/mmPP;

% Compute the G2 band
[G2_1, DAPI_G2_mode] = TiffImg.decimate_and_smooth(x(idx), y(idx), dapi_c(idx), options);

if DEBUG
    figure
    tri = delaunay(G2_1.X,G2_1.Y);
    tsh = trisurf(tri,G2_1.X,G2_1.Y,G2_1.Z);
    shading interp
    tsh.EdgeColor = 'k';

    line(x, y, dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('G2 surface fit')
    setTheme(gcf,'dark')
    axis tight
    axis square
    zlim([0.5,2.8])
    view(0,90)
end


% This is only an approximate correction that will be used for displaying
% results if debugging
dapi_c = ((dapi_c-1) ./ (DAPI_G2_mode(x,y)/2)) + 1; % dapi_c = dapi_c ./ (DAPI_G2_mode(x,y)/2); % 

% We want to have the G1 band at 1 and the G2 band at 2; therefore, the
% correction surface is
%
% (dapi_norm - 1) = (dapi - g1) * (2 - 1) / (g2 - g1)
% dapi_norm = dapi / (g2 - g1) - g1 / (g2-g1) + 1
%
% FG_f <-- (g2 - g1)
% FG_o <--1 - g1/(g2-g1)


if DEBUG
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('Data after correction')
    setTheme(gcf,'dark')
    zlim([0,3])

    figure
    ax = axes;
    hold on
    h = histogram(ax,dapi_c,'BinLimits',[0,6]);
    grid on
    setTheme(gcf,'dark')
end

% Combine all results into a single smooth surface and stripe. -----------



% Construct the flattening surface

% Each flattening surface could have a different number of points. If this
% is the case, then they need to be interpolated onto a grid of the same
% size before combining them.

% Grid to re-interpolate results on - use the grid for creating the
% background.
xg = tiffImg.BG_offset.x;
yg = tiffImg.BG_offset.y;

if isequal(G1_1.X,G1_2.X) && isequal(G1_2.X,G2_1.X)
    
    % Combine the two G1 surfaces
    Z = G1_1.Z .* G1_2.Z;
    X = G1_1.X;
    Y = G1_1.Y;
    
    % Scattered interpolants are quite slow to evaluate, so reinterplate
    % onto the same square grid
    fun = scatteredInterpolant(double(X),double(Y),double(Z));
    G1z = fun({xg,yg})';
    
    fun = scatteredInterpolant(double(X),double(Y),double(G2_1.Z));
    G2z = fun({xg,yg})';
else
    fun = scatteredInterpolant(double(G1_1.X),double(G1_1.Y),double(G1_1.Z));
    G1z1 = fun({xg,yg})';
    
    fun = scatteredInterpolant(double(G1_2.X),double(G1_2.Y),double(G1_2.Z));
    G1z2 = fun({xg,yg})';
    
    fun = scatteredInterpolant(double(G2_1.X),double(G2_1.Y),double(G2_1.Z));
    G2z = fun({xg,yg})';
    
    G1z = G1z1 .* G1z2;
end

% G2z is currently relative to G1z, make it so that it is not
G2z = G2z.*G1z;

FG_f.Z = (G2z-G1z);
FG_f.x = xg;
FG_f.y = yg;

FG_o = FG_f;
FG_o.Z = 1 - G1z./(G2z-G1z);

%% Compute the G1 area surface --------------------------------------------
if nargout > 3
    nucleiPerBin = 800; % emperical
    bin = sqrt( (2/sqrt(3)) * nucleiPerBin / rho );
    smooth = 2.5*bin;
    idx = dapi_c > 0.7 & dapi_c < 1.3;
    options.defaultValue = nan;
    options.binSize = [bin/mmPP,1]; % [2.5/mmPP,1];
    options.smoothingRadius = smooth/mmPP; % 6/mmPP;
    options.reductionMethod = 'median';
    [G1_area, G1_area_fun]  = TiffImg.decimate_and_smooth(x(idx), y(idx), area(idx), options);

    area_c = area./G1_area_fun(x,y);
    [G1_stripe, G1_stripeX, ~, fsample] = Fit_Stripe_Artifact(x(idx)*mmPP,area_c(idx),'Threshold',1e-2,'generatePlots',DEBUG);
    
    % interpolate onto the same square grid as the threshold
    fun = scatteredInterpolant(double(G1_area.X), double(G1_area.Y), double(G1_area.Z));
    Z = fun({xg,yg})';
    
    if DEBUG
        figure
        tri = delaunay(G1_area.X,G1_area.Y);
        trisurf(tri,G1_area.X,G1_area.Y,G1_area.Z);
        line(x(idx),y(idx),area(idx),'marker','.','linestyle','none','color','g','markersize',1)
        title('Area surface')
        setTheme(gcf,'dark')
        axis tight
    end
    
    G1Area.Z = Z;
    G1Area.x = xg;
    G1Area.y = yg;
end

if nargout > 4
    G1_idx = dapi_c > 0.7 & dapi_c < 1.3;
end

if DEBUG
        
    FG_f_fun = @(x,y) interp2mex(FG_f.Z, nakeinterp1(xg(:),(1:size(FG_f.Z,2))',x), nakeinterp1(yg(:),(1:size(FG_f.Z,1))',y));%griddedInterpolant({yg,xg},Z,'linear','nearest');%
    FG_o_fun = @(x,y) interp2mex(FG_o.Z, nakeinterp1(xg(:),(1:size(FG_o.Z,2))',x), nakeinterp1(yg(:),(1:size(FG_o.Z,1))',y));%griddedInterpolant({yg,xg},Z,'linear','nearest');%
    
    if isempty(Xstripe)
        dapi_c = (dapi) ./ (FG_f_fun(x,y)) + FG_o_fun(x,y);
    else
        FG_s = @(x) nakeinterp1((1:tiffImg.imageSize(2)).', Xstripe(:), x);
        dapi_c = (dapi) ./ (FG_s(x).*FG_f_fun(x,y)) + FG_o_fun(x,y);
    end
    
    figure
    line(x,y,dapi_c,'marker','.','linestyle','none','color','g','markersize',1)
    title('Data when applying correction all at once')
    setTheme(gcf,'dark')
    axis tight
    zlim([0,3])
    
    histogram(ax,dapi_c,'BinEdges',h.BinEdges,'FaceColor','r')
end
end



% Remove any stripe artifact from the G2 band. ---------------------------

% Select the g2 band. Compute the median dapi value along small x-slices to
% extract the stripe. Divide the stripe away.

% idx = dapi_c > 1.7 & dapi_c < 2.3;
% G2_stripe = decimateData(x(idx),ones(sum(idx),1),dapi_c(idx),'binSize',[stripeBinWidth,100],'defaultValue',2);
% dapi_c = dapi_c ./ (nakeinterp1(G2_stripe.X(:,1), G2_stripe.Z(:,1), x)/2) ;

% idx = dapi_c > 0.7 & dapi_c < 1.3; % Get G1 band
% G1_stripe3 = decimateData(x(idx),ones(sum(idx),1),dapi_c(idx),'binSize',[stripeBinWidth,100],'defaultValue',1);  % Get median dapi value from bins 100 pixels wide along x direction
% G1_stripe3.Z = highpass(G1_stripe3.Z(:,2),2.5,100*mmPP) + 1;
% dapi_c = dapi_c ./ nakeinterp1(G1_stripe3.X(:,1), G1_stripe3.Z(:,1), x); % Divide out the median value
%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
