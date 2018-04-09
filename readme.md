# Cell image analysis
A library of functions and classes to speed up the analysis of large images containing stained (immunofluorescence) cell nuclei (images ~20 mm x 20 mm with ~200,000 cells). There are methods for segmenting nuclei, image background and artifact correction, and feature extraction. Methods use GPU acceleration for image processing and features of each nuclei can be extracted in parallel on a GPU.

![GitHub Logo](/docs/processing_img1.png)

A high level overview of the classes and methods can be seen below. Note that several of the functions are specific for images where the nuclei DNA is stained, such as when using DAPI or PI. Further several methods are specific to processing images created using whole slide fluorescence readers.

All code is written for use in Matlab. The tiff images are assumed to have been created with a Leica Aperio fluorescence slide scanner. These tiff images have their acquisition info stored in the `ImageDescription` tag (id: 270) of the tiff images. A requried piece of information in this information is `MPP` which gives the spatial scale of the image in microns-per-pixels

## Example usage

```Matlab
% Create an instance of the CellExperiment class. We have to image channels DAPI, and
% GFP. Note, a DAPI channel is required for this class.
% cE = CellExperiment(channel_names,channel_paths)
cE = CellExperiment(['DAPI','GFP'], ['/example/DAPI.tiff','/example/GFP.tiff'])

%% Setup image correction functions for each channel.
% The corrections applied to each channel can be controlled by inputing
% algebraic expressions.
%
% Set DAPI correction:
eE.Channel_TiffImgs(1).Image_Correction_Expression = "(S-BG_o*BG_s)/(FG_s*FG_f)";
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
eE.Channel_TiffImgs(1).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "FG_o";
% Here we subtract away a forground ofset correction (FG_o) which will place the DAPI G1
% band at value 1 and the DAPI G2 band at value 2.

% Set GFP image corrections
eE.Channel_TiffImgs(2).Image_Correction_Expression = "(S-BG_o*BG_s)/FG_s";
eE.Channel_TiffImgs(2).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "-FG_o";

%% Setup object segmentation
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

%% Setup feature extraction
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

% With everything setup, process the images.

%% Processes
% Compute object masks
Create_Object_Mask(cE, imgPth(trials(i),"Mask"));

% Background correct the channels
cE.Surface_Smoothing_Radius = nan;
Correct_Image_Backgrounds(cE,true);

% Extract features
pth = 'FeaturesFile.mat'; % path to save features
Compute_Object_Features(cE, pth);
```
<!--
## Class `TiffImg`
This class is a wrapper for the `tifflib` library. This class is much faster at repeated reading of small portions of large tiff files than Matlab's built in `Tiff` class.

### High level methods

##### Load tif image
```Matlab
t = TiffImg('path-to-img.tiff');
```

##### Threshold:
```Matlab
% Compute threshold for object segmentation
th = t.Compute_Threshold();
```
- This uses an adaptive log-weighted otsu-thresholding similar to CellProfiler software.
- The threshold is computed for each block of the image (specified by t.blockSize)
- When applying the threshold to the image, the computed threshold values from each block
  are linearly interpolated to get the threshold value at each pixel.
- The computed threshold is stored in `t.threshold` and is optionally

##### Compute background
```Matlab
% Compute image background
bg = t.Compute_Background();
```
- Background is determined by 1) creating a mask with all pixels less than the threshold 2) eroding this mask with a small diamond (radius 2) 3) computing the median image intensity of all mask pixels from each block of the image.
- Background is stored in `t.BG_offset`

The images aquired with the Aperio slide reader can have a stripe artifact along the slow scan axis. This stripe artifact can be determined using
```Matlab
% Compute background stripe
bg_s = t.Compute_StripeArtifact('background');
```
- The background stripe is computed relative to the background. That is, after computing the background, the image is _divided_ by the background, and then the background stripe is computed.
- The background stripe is stored in `t.BG_stripeX`

The background and background stripe can be computed with a single function call using
```Matlab
% Compute background and background stripe
bg = t.Compute_Background(true);
```

-->
