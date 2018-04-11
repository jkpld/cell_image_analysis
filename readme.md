# Cell image analysis
A library of functions and classes to speed up the analysis of large images containing stained (immunofluorescence) cell nuclei (images ~20 mm x 20 mm with ~200,000 cells). There are methods for segmenting nuclei, image background and artifact correction, and feature extraction. Methods use GPU acceleration for image processing and features of each nuclei can be extracted in parallel on a GPU.

![image_processing](/docs/processing_img1.png)


## Requirements
* Matlab version R2017b or higher.
* The Image Processing Toolbox.
* The Statistics and Machine Learning Toolbox.
* The Parallel Computing Toolbox - _This toolbox is required to use parallel computing and GPU acceleration._

#### Assumptions
Several of the functions are specific for images where the nuclei DNA is stained, such as when using DAPI or PI; and the `CellExperiment` class requires a channel labeled DAPI. Further several methods are specific to processing images created using whole slide fluorescence readers.

The tiff images are assumed to have been created with a Leica Aperio fluorescence slide scanner. These tiff images have their acquisition info stored in the `ImageDescription` tag (id: 270) of the tiff images. A required piece of information in this information is `MPP` which gives the spatial scale of the image in microns-per-pixels.

## Installation
1. Download the repository.
2. Open Matlab on your computer, and change the current folder to the downloaded repository.
3. Run `setup.m`. This function will
    * Add the necessary files to your Matlab path.
    * Attempt to compile two C functions that will increase the program speed. If a C compiler is not found on your computer, then the code will still run, but it could be slower.
4. Start by browsing the example located in the `example/` folder, and described below.

#### Installation Notes
* The GitHub repository could take a few minutes to download as the repository is >100 MB. The majority of this space comes from the example images included.
* `setup.m` will not permanently add the files to your Matlab path. If you want the files to remain on your Matlab path, the run the command `savepath` in the Matlab command window after running `setup.m`..


## Example usage
The example image (~4x4 mm<sup>2</sup>) is a small portion of a full slide image with ~20000 cells. In this example, the nuclei will be segmented and partitioned, the images will be corrected, and features extracted. Finally, we will look at a summary of the computed features, and plot several nuclei likely to be in the M-phase.

### Load images

Create an instance of the CellExperiment class. We have to image channels DAPI, and GFP. Note, a DAPI channel is required for this class.
```Matlab
cE = CellExperiment(channel_names,channel_paths)
cE = CellExperiment(["DAPI","GFP"], ["example\ca_DAPI.tif","example\ca_GFP.tif"]);
```
### Setup processing, object partitioner, and features to extract

Set image correction functions for each channel.
The corrections applied to each channel can be controlled by inputing
algebraic expressions.

#### Set DAPI correction:
```Matlab
cE.Channel_TiffImgs(1).Image_Correction_Expression = "(S-BG_o*BG_s)/(FG_s*FG_f)";
```
This subtracts away the background and divides by the foreground correction

* `S` : our image
* `BG_o` : image background
* `BG_s` : image background stripe - the Aperio slide reading can have a stripe artifactalong the slow scan axis. This stripe artifact is computed relative to the image background, therefore the total background is BG_o * BG_s.
* `FG_f` : foreground correction factor - for the DAPI channel, this factor will flatten the G1 and G2 bands across the image.
* `FG_s` : foreground stripe (same as BG_s)

After extracting the channel intensities for each object, we can apply an additional
correction to the intensities.
```Matlab
cE.Channel_TiffImgs(1).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "FG_o";
```
Here we add away a foreground offset (FG_o) which will place the DAPI G1 band at value 1 and the DAPI G2 band at value 2.

#### Set GFP image corrections
```Matlab
cE.Channel_TiffImgs(2).Image_Correction_Expression = "(S-BG_o*BG_s)/FG_s";
cE.Channel_TiffImgs(2).ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression = "-FG_o";
```

### Setup object segmentation
Object segmentation is performed using the DAPI channel. Segmentation has two steps.
- Segment all nuclei using adaptive log-weighted otsu thresholding
- De-clump (or partition) all nuclei that are overlapping with each other.

In this example, objects will be partitioned using watershed on the distance transform.
```Matlab
partitioner = DTWatershed;
partitioner.Use_Parallel = true; Set to use parallel processing
```

### Setup feature extraction
Features are extracted using a `FeatureExtractor` class. Different `FeatureGroups` are added to the `FeatureExtractor`. Here I extract, location and shape features from the objects and intensity, radial intensity, haralick texture, and granularity features from each channel.

#### Create `FeatureGroups`
```Matlab
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
```

### Set the partitioner and featureExtractor
```Matlab
cE.nucleiPartitioner = partitioner;
cE.featureExtractor.featureGroups = featGroups;
```

### Process the iamge
With everything setup, process the images.

Compute object masks
```Matlab
pth = 'example\Mask.tif';
Create_Object_Mask(cE, pth);
```
Background correct the channels
```Matlab
cE.Surface_Smoothing_Radius = nan;
Correct_Image_Backgrounds(cE,true);
```

Extract features
```Matlab
pth = 'K:\RadiationLab\cell_image_analysis\example\FeaturesFile.mat'; path to save features
Compute_Object_Features(cE, pth);
```

### Create summary
- Plot the distribution of DAPI intensity vs. normalized area (shown below) - you can
qualitatively extract cell-phase from this distribution.
- Plot the x-y distributions of the DAPI and GFP intensities.
- Show a summary of some statistics for each feature extracted.

```Matlab
cE.Create_Summary()
```
![area_vs_dapi](/docs/area_vs_dapi_distribution.PNG)

### Look at some cells likely in M-phase

Load up the features
```Matlab
[X,Xn] = Load_Features(cE);
```

Get the normalized nuclei area and DAPI intensity
```Matlab
a = X(:,Xn=="Shape_Area") / median(cE.DAPI_G1_Area.Z(:)); % normalized area
d = X(:,Xn=="Intensity_DAPI_Integrated"); % dapi
```

Select a random 10 nuclei that are likely in M-phase
```Matlab
inRange = find(d<2.3 & d>2 & a<1.15 & a>0.9); % region with high probability of finding cells in M-phase
idx = inRange(randperm(numel(inRange),10)); % Randomly select 10 cells
```

Extract the nuclei images
```Matlab
objImgs = cE.Extract_Object_Images(X(idx,1:2), ["DAPI","GFP","Mask"], [71,71]);
```

View the nuclei
```Matlab
viewNucleiImages(objImgs,["DAPI","GFP","Mask"])
```
![example_nuclei](/docs/example_nuclei.PNG)
