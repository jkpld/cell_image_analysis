function [stripe, xCents, z, fsample] = Fit_Stripe_Artifact(x,z,varargin)
% FIT_STRIPE_ARTIFACT will fit the stripe artifact along the slow
% axis in two steps. First it will fit it with using a very small binSize
% for x, then it will determine a new bin size based on the frequency
% spectrum of the first fit. Using the new bin size it will refit the
% stripe artifact. This is necessary because if the sample frequency is
% too much larger (binsize to small) than necessary, you end up fitting
% noise and this results in a distinct artifact.
%
% This function should only be called after first fitting the data to a
% surface and dividing the data by the surface to flatten out the data --
% i.e. bg_fitG1Band should be called and the result used to flatten the
% data.
%
% Input parameters
% x - Spatial x data. Data should be in [mm]
%
% z - Data to fix optical correction by, for example DAPI or GFP intensity.
% In general, each color channel of the image needs to be corrected.
%
% Optional parameters:
% initialBinSize - Size of the initial bin used to discretize x data. Unit
% should be mm. Default is 0.005 mm. 
%
% z_bounds - Bounds of the z data to use when fitting the optical
% distortion. Default is [0.65, 1.5]. 
%
% threshold - Threshold value of when to stop fitting the optical
% distortion. Default is 1e-4. --- Note that decreasing this value could
% lead to an artifact where the histogram of z shows a distinct minimum at
% the upper sample range.
%
% generatePlots - logical flag. if true, then some plots will be made, if
% false, then no plots will be made.
%
% Output parameters:
% z - corrected data
%
% distortion - the amplitude of the stripe artifact
%
% xGrid - the x values used for calculated the distortion.
%
% fsample - the sample frequency of the final fit. The grid spacing is
% 1/fsample
%
% JamesKapaldo
% August 10, 2016

% Make sure inputs are columns.
if isrow(x)
    x = x(:);
end
if isrow(z)
    z = z(:);
end

% Get input values
p = inputParser;
p.FunctionName = 'bg_fitOpticalDistortion';
addParameter(p,'initialBinSize',0.005,@(t) numel(t) == 1 && t > 0 && t < max(x)-min(x));
addParameter(p,'z_bounds',[0.65, 1.5],@(t) numel(t) == 2 && all(t > 0 & t < max(z)-min(z)) && t(1) < t(2));
addParameter(p,'threshold',1e-4,@(t) numel(t) == 1 && t > 0);
addParameter(p,'generatePlots',false,@(t) islogical(t) || (t==0 || t == 1) );

parse(p,varargin{:})

voxelSizes = p.Results.initialBinSize; % Sample period
zB = p.Results.z_bounds;
threshold = p.Results.threshold;
generatePlots = p.Results.generatePlots;

boundEdges = prctile(x,[0,100]);

xEdges = (boundEdges(1)-voxelSizes(1)) : voxelSizes(1)  : (boundEdges(2) + voxelSizes(1));

% The number of xCents must be odd for frequency calculation
if mod(numel(xEdges),2) == 0
    xEdges(end+1) = xEdges(end)+voxelSizes;
end

xCents = xEdges(1:end-1) + diff(xEdges)/2;

fsample = 1/voxelSizes; % sample frequency
L = numel(xCents);
freq = fsample*(0:(L/2))/L;


stripe = ones(numel(xCents),1);

%     figure
%     axhist = axes;
%     hold on

G1_err = inf;
counter = 1;
origZ = z;
while G1_err > threshold
    
    if counter > 10
        warning('Optical distortion correction reached maximum iteration limit.')
        break;
    end
    
    G1_ind = z > zB(1) & z < zB(2);% + rand(1)*0.2 - 0.1;
    G1_z = z(G1_ind);
    bin = discretize(x(G1_ind),xEdges);
    
    optDist = accumarray(bin(all(bin>0,2),:),G1_z(all(bin>0,2)),[numel(xEdges)-1,1],@median,1);

    stripe = stripe.*optDist;
    crrctn = interp1(xCents,optDist,x,'linear');
    crrctn(isnan(crrctn) | isinf(crrctn) | crrctn == 0) = 1;
    
    newz = z./crrctn;
    
    G1_err = sqrt(sum((newz(G1_ind)-z(G1_ind)).^2))/sum(G1_ind);
    
    if generatePlots
        fprintf('difference (%d) = %0.5g \n', counter, G1_err)
    end
    
    z = newz;
    
    counter = counter + 1;
end

stripe(isnan(stripe) | isinf(stripe) | stripe == 0) = 1;

Y = fft(stripe);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%     figure
%     line(freq, P1)
%     ylim([0,0.1])


% Now repeate with proper sample frequency -------------------------------

z = origZ;

% fsample = 2*(max(freq(P1>0.01)) + 40); % old sample frequency - this is not a harmonic of the stripe frequency
fsample = 3*max(freq(P1>median(P1(2:end))+2*std(P1(2:end)))); % new sample frequency - this should be a harmonic of the stripe frequency

voxelSizes = 1/fsample;
xEdges = (boundEdges(1)-voxelSizes(1)) : voxelSizes(1)  : (boundEdges(2) + voxelSizes(1));

if mod(numel(xEdges),2) == 0
    xEdges(end+1) = xEdges(end)+voxelSizes;
end

xCents = xEdges(1:end-1) + diff(xEdges)/2;
stripe = ones(numel(xCents),1);

G1_err = inf;
counter = 1;

while G1_err > threshold
    
    if counter > 10
        warning('Optical distortion correction reached maximum iteration limit.')
        break;
    end
    
    G1_ind = z > zB(1) & z < zB(2);% + rand(1)*0.2 - 0.1;
    G1_z = z(G1_ind);
    bin = discretize(x(G1_ind),xEdges);
    
    optDist = accumarray(bin(all(bin>0,2),:),G1_z(all(bin>0,2)),[numel(xEdges)-1,1],@median,1);

    stripe = stripe.*optDist;
    crrctn = interp1(xCents,optDist,x,'linear');
    crrctn(isnan(crrctn) | isinf(crrctn) | crrctn == 0) = 1;
    
    newz = z./crrctn;
    
    G1_err = sqrt(sum((newz(G1_ind)-z(G1_ind)).^2))/sum(G1_ind);
    
    if generatePlots
        fprintf('difference (%d) = %0.5g \n', counter, G1_err)
    end
    
    z = newz;
    
    counter = counter + 1;
end

stripe(isnan(stripe) | isinf(stripe) | stripe == 0) = 1;

if generatePlots
    figure;
    subplot(2,1,1)
    line(x,origZ,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    line(xCents,stripe,'Color','b','LineWidth',1)
    drawnow;
    axis tight
    ylim([0.5,2.5]);
    
    subplot(2,1,2)
    line(x,z,'Marker','.','MarkerSize',1,'LineStyle','none','Color','y')
    drawnow;
    setTheme(gcf,'dark')
    axis tight
    ylim([0.5,2.5]);
end

end