pth = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_GFP.tif';
pth_n = 'K:\RadiationLab\Cell_Incubation\Cancer\p24_GFP_small.tif';
%%
t = TiffImg(pth);
%%
internal_Write_Small_Image(t, pth_n, [10,18], [10,18])
disp here
%%

% Create an instance of the CellExperiment class. We have to image channels DAPI, and
% GFP. Note, a DAPI channel is required for this class.
% cE = CellExperiment(channel_names,channel_paths)
cE = CellExperiment(["DAPI","GFP"], ["K:\RadiationLab\cell_image_analysis\example\DAPI.tif","K:\RadiationLab\cell_image_analysis\example\GFP.tif"]);
%%
% Setup image correction functions for each channel.
% The corrections applied to each channel can be controlled by inputing
% algebraic expressions.
%
% Set DAPI correction:
cE.Channel_TiffImgs(1).Image_Correction_Expression = "(S-BG_o*BG_s)/(FG_s*FG_f)";
% This subtracts away the background and divides by the foreground correction
% S : our image
% BG_o : image background
% BG_s : image background stripe - the Aperio slide reading can have a stripe artifact
%  along the slow scan axis. This stripe artifact is computed relative to the image
%  background, therefore the total background is BG_o * BG_s.
% FG_f : foreground correction factor - for the DAPI channel, this factor will flatten the
%  G1 and G2 bands across the image.
% FG_s : foreground stripe (same as BG_s)
%
% After extracting the channel intensities for each object, we can apply an additional
% correction to the intensities.
cE.Channel_TiffImgs(1).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "FG_o";
% Here we subtract away a forground ofset correction (FG_o) which will place the DAPI G1
% band at value 1 and the DAPI G2 band at value 2.

% Set GFP image corrections
cE.Channel_TiffImgs(2).Image_Correction_Expression = "(S-BG_o*BG_s)/FG_s";
cE.Channel_TiffImgs(2).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "-FG_o";
% Setup object segmentation
% Object segmentation is performed using the DAPI channel. Segmentation has two steps.
% - Segment all nuclei using adaptive log-weighted otsu thresholding
% - De-clump (or partition) all nuclei that are overlapping with each other.
%
% Nuclei partitioning (step two) is handled by adding an ObjectPartitioner class object
% to the CellExperiment. There are currently three ObjectPartitions
%  1. None : no object partitioning is performed (default)
%  2. DTWatershed : Watershed partitioning is used on the distance transform of the
%     objects.
%  3. SALRGeoPartition : Use SALR clustering to locate nuclei centers and then use
%     a seed-point based geometric partitioning. This method requires that SALR clustering
%     and geometric partitioning toolboxes have been installed.
partitioner = DTWatershed;
partitioner.Use_Parallel = true; % Set to use parallel processing

% Setup feature extraction
% Features are extracted using a FeatureExtractor class. Different FeatureGroups are added
% to the FeatureExtractor. Here I extract, location and shape features from the objects
% and intensity, radial intensity, haralick texture, and granularity features from each
% channel.

% Create FeatureGroups
defaultOptions = struct(...
    'Use_GPU', 1, ...
    'Use_Parallel', 0, ...
    'Number_Boundary_Vertices', 100, ...
    'Number_Fourier_Descriptors', 20, ...
    'Number_of_Rings', 3, ...
    'Ring_Width', 4, ...
    'PixelOffset', [3,7], ...
    'SubSampleSize', 1, ...
    'BackgroundSampleSize', 1, ...
    'ElementSize', 10, ...
    'GranularSpectrumLength', 7);

warning('off','FeatureGroup:extraOptions')

clear featGroups
featGroups(10) = EmptyFeature;
featGroups(1) = Location('', defaultOptions);
featGroups(2) = Shape('', defaultOptions);
featGroups(3) = Intensity('DAPI', defaultOptions);
featGroups(4) = Intensity('GFP', defaultOptions);
featGroups(5) = RadialIntensity('DAPI', defaultOptions);
featGroups(6) = RadialIntensity('GFP', defaultOptions);
featGroups(7) = HaralickTexture('DAPI', defaultOptions);
featGroups(8) = HaralickTexture('GFP', defaultOptions);
featGroups(9) = Granularity('DAPI', defaultOptions);
featGroups(10) = Granularity('GFP', defaultOptions);
featGroups(featGroups.isempty()) = [];
% Set the partitioner and featureExtractor
cE.nucleiPartitioner = partitioner;
cE.featureExtractor.featureGroups = featGroups;

%% With everything setup, process the images.

% Processes
% Compute object masks
pth = 'K:\RadiationLab\cell_image_analysis\example\Mask.tif';
Create_Object_Mask(cE, pth);
%
% Background correct the channels
cE.Surface_Smoothing_Radius = nan;
Correct_Image_Backgrounds(cE,true);
%
% Extract features
fprintf('Computing features...\n')
pth = 'K:\RadiationLab\cell_image_analysis\example\FeaturesFile.mat'; % path to save features
Compute_Object_Features(cE, pth);

%% Look at some cells in M-phase

% Load up the features
[X,Xn] = Load_Features(cE);

% Normalize area
a = X(:,Xn=="Shape_Area") / median(cE.DAPI_G1_Area.Z(:));
d = X(:,Xn=="Intensity_DAPI_Integrated");            
inRange = find(d<2.1 & d>1.9 & a<1.1 & a>0.9);% & (X(:,1)*mpp>4 & X(:,1)*mpp<14) & (X(:,2)*mpp>4 & X(:,2)*mpp<14);

idx = inRange(randperm(numel(inRange),10));

objImgs = cE.Extract_Object_Images(X(idx,1:2), ["DAPI","GFP","Mask"], [71,71]);


viewNucleiImages(objImgs,["DAPI","GFP","Mask"])