%% Process example image.
% The image (~4x4 mm^2) is a smal portion of a full slide image with ~20000
% cells. In this example, the nuclei will be segmented and partitioned, the
% images will be corrected, and features extracted. Finally, we will look
% at a summary of the computed features, and plot several nuclei likely to
% be in the M-phase.

%% Load images

% Create an instance of the CellExperiment class. We have to image channels DAPI, and
% GFP. Note, a DAPI channel is required for this class.
% cE = CellExperiment(channel_names,channel_paths)
cE = CellExperiment(["DAPI","GFP"], ["example\ca_DAPI.tif","example\ca_GFP.tif"]);

%% Setup processing, object partitioner, and features to extract

% Set image correction functions for each channel.
% The corrections applied to each channel can be controlled by inputing
% algebraic expressions.

% Set DAPI correction:
cE.Channel_TiffImgs(1).Image_Correction_Expression = "(S-BG_o*BG_s)/(FG_s*FG_f)";
% This subtracts away the background and divides by the foreground correction
% S : our image
% BG_o : image background
% BG_s : image background stripe - the Aperio slide reading can have a
%  stripe artifactalong the slow scan axis. This stripe artifact is
%  computed relative to the image background, therefore the total
%  background is BG_o * BG_s.
% FG_f : foreground correction factor - for the DAPI channel, this factor
%  will flatten the G1 and G2 bands across the image.
% FG_s : foreground stripe (same as BG_s)
%
% After extracting the channel intensities for each object, we can apply an additional
% correction to the intensities.
cE.Channel_TiffImgs(1).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "FG_o";
% Here we add away a forground offset (FG_o) which will place the DAPI G1
% band at value 1 and the DAPI G2 band at value 2.

% Set GFP image corrections
cE.Channel_TiffImgs(2).Image_Correction_Expression = "(S-BG_o*BG_s)/FG_s";
cE.Channel_TiffImgs(2).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "-FG_o";

% Setup object segmentation
% Object segmentation is performed using the DAPI channel. Segmentation has
% two steps.
% - Segment all nuclei using adaptive log-weighted otsu thresholding
% - De-clump (or partition) all nuclei that are overlapping with each other.
%
% In this example, objects will be partitioned using watershed on the
% distance transform.
partitioner = DTWatershed;
partitioner.Use_Parallel = true; % Set to use parallel processing

% Setup feature extraction
% Features are extracted using a FeatureExtractor class. Different
% FeatureGroups are added to the FeatureExtractor. Here I extract, location
% and shape features from the objects and intensity, radial intensity,
% haralick texture, and granularity features from each channel.

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

%% Process the iamge
% With everything setup, process the images.

% Compute object masks
pth = 'example\Mask.tif';
Create_Object_Mask(cE, pth);

% Background correct the channels
cE.Surface_Smoothing_Radius = nan;
Correct_Image_Backgrounds(cE,true);

% Extract features
pth = 'example\FeaturesFile.mat'; % path to save features
Compute_Object_Features(cE, pth);

% Create summary
% - Plot the distribution of DAPI intensity vs. normalized area - you can
% qualitatively extract cell-phase from this distribution.
% - Plot the x-y distributions of the DAPI and GFP intensities.
% - Show a summary of some statistics for each feature extracted.
cE.Create_Summary()

%% Look at some cells likely in M-phase

% Load up the features
[X,Xn] = Load_Features(cE);

a = X(:,Xn=="Shape_Area") / median(cE.DAPI_G1_Area.Z(:)); % normalized area
d = X(:,Xn=="Intensity_DAPI_Integrated"); % dapi

inRange = find(d<2.3 & d>2 & a<1.15 & a>0.9); % region with high probability of finding cells in M-phase
idx = inRange(randperm(numel(inRange),10)); % Randomly select 10 cells

% Extract the cell images
objImgs = cE.Extract_Object_Images(X(idx,1:2), ["DAPI","GFP","Mask"], [71,71]);

viewNucleiImages(objImgs,["DAPI","GFP","Mask"])
