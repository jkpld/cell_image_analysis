# Class `TiffImg`
This class is a wrapper for the `tifflib` library. This class is much faster at repeated reading of small portions of large tiff files than Matlab's built in `Tiff` class. Additionally, the processing can be accelerated by using a GPU.

## List of properties
See the class definition for a list of all properties. Here a truncated list of properties is shown.
* `FileName` - Path to the .tif image
* `Acquisition_Info` - Structure containing the information stored in the `ImageDescription` tag of the .tif image. It is assumed to have a `key = value` structure, with entries separated by at `|`.
* `mmPerPixel` - milimeters per pixel. This value is derived from the `MPP` key in the `ImageDescription`; this key is required.
* `blockSize` - The image block size used when processing the image. This size must be a multiple of the .tif image tile size.
* `Surface_Smoothing_Radius` - The radius (in mm) to use when smoothing the surfaces computed (threshold, background, foreground). If set to NaN, then smoothing will be disabled.
* `Image_Smooth_Kernel` - The Kernel used to smooth the image before all operations (directly after reading the image in.)
* `Image_Correction_Expression` - The expression used to correct the image. Valid variable names are
  - `S` : The input image
  - `BG_o` : A background offset (`BG_offset`) that is added or subtracted
  - `BG_s` : A background x stripe (`BG_stripeX`) that is multiplied or divided
  - `FG_o` : A foreground offset (`FG_offset`) that is added or subtracted
  - `FG_f` : A foreground factor (`FG_factor`) that is multiplied or divided
  - `FG_s` : A foreground x stripe (`FG_stripeX`) that is multiplied or divided

  If `BG_offset` or `FG_offset` have not been assigned to the `TiffImg` yet, but they show up in the expression, they are replaced with 0's.

  If `BG_factor`, `FG_factor`, `BG_stripeX`, or `FG_stripeX` have not been assigned to the `TiffImg` yet, but they show up in the expression, they are replaced with 1's.

  Example expression:
        (S - BG_o * BG_s + FG_o) / (FG_s * FG_f)

* `ObjectIntegratedIntensity_FeaturePostProcess_OffsetExpression` - An expression, as in `Image_Correction_Expression`, that corrections the object intensities after being computed.
* `Use_GPU` - Logical flag determining if a GPU should be used to accelerate computation
* `Verbose` - Logical flag determining if updates are displayed while computing
* `threshold` - Structure (fields: X, Y, Z) storing the `threshold`
* `BG_offset` - Structure (fields: X, Y, Z) storing the `BG_offset`
* `BG_stripeX` - Row array storing the `BG_stripeX`
* `FG_offset` - Structure (fields: X, Y, Z) storing the `FG_offset`
* `FG_factor` - Structure (fields: X, Y, Z) storing the `FG_factor`
* `FG_stripeX` - Row array storing the `FG_stripeX`
* `Sharpness` - Structure (fields: X, Y, Z) storing the `Sharpness`


## List of high level methods
* `TiffImg` - Create object instance
```Matlab
t = TiffImg('path-to-img.tiff');
```

* `Compute_Threshold` - Compute image threshold
```Matlab
th = t.Compute_Threshold()
```
 - Threshold is computed in each image block using adaptive log-weighted Otsu thresholding
 - Threshold is stored in `t.threshold`.
 - When the threshold is applied it is linearly interpolated to get the threshold value at each pixel of the image.


* `Compute_Background` - Compute image background
```Matlab
bg = t.Compute_Background();
```
 - Background is determined by
   1. Creating a mask with all pixels less than the threshold
   2. Eroding this mask with a small diamond (radius 2)
   3. Computing the median image intensity of all mask pixels from each block of the image.
 - Background is stored in `t.BG_offset`


 * `Compute_Foreground` - Compute image foreground
 ```Matlab
 fg = t.Compute_Foreground();
 ```
  - Foreground is determined by
    1. Creating a mask with all pixels greater than the threshold
    2. Eroding this mask with a small diamond (radius 2)
    3. Computing the median image intensity of all mask pixels from each block of the image.
  - Foreground is stored in `t.FG_factor`


  * `Compute_StripeArtifact` - Compute stripe artifact in image background or foreground. When using the Aperio slide reader, there is often a strip artifact that can be seen along the slow scan axis.
  ```Matlab
  st = t.Compute_StripeArtifact(BG_or_FG)
  ```
   - Stripe is computed by
     1. Compute the background or foreground mask (as in `Compute_Background`/`Compute_Foreground`)
     2. Compute the median masked image intensity along entire rows of the image. (This assumes that the rows are the fast axis of the scanner.)


  * `Compute_Sharpness` - Compute image sharpness
  ```Matlab
  fg = t.Compute_Sharpness();
  ```
   - Sharpness is computed in each image block as the number of pixels with a fourier transform magnitude greater than the maximum fourier transform magnitude over 1000, and normalized by the square root of the number of pixels.


* `Measure_BasicProps` - Measure the basic properties (location, intensity, area) for each object
```Matlab
[feature,feature_names] = t.Measure_BasicProps(channelName)
```


* `Measure_Intensity` - Measure only the object intensities
```Matlab
feature = t.Measure_Intensity(useParallel)
```


* `Write_Object_Mask` - Write the object mask to file, and optionally apply partitioning before writing
```Matlab
t.Write_Object_Mask(fileName, partitioner)
```
- `partitioner` is an instance of class `ObjectPartitioner`


* `Compute_Channel_Corrections` - Compute the intensity and stripe corrections based on object intensities for general image channels.
```Matlab
[flatteningSurf, xStripe] = t.Compute_Channel_Corrections(x, y, I);
```
 - The background computation (`flatteningSurf`) uses the bottom 2%-4% of the object intensities.
 - The stripe is computed relative to the object intensities (thus, it should be multiplied or divided).


* `Compute_DAPI_Corrections` - Compute the intensity and strope corrections based on object intensities for the DAPI channel (or any channel that stains all of the nuclei DNA, e.g. PI)
```Matlab
[FG_f, FG_o, Xstripe, G1Area, G1_idx] = t.Compute_DAPI_Corrections(x,y,I,A)
```
 - Fit the G1 band using the median intensity, and then fit the G2 band
 - Compute corrections to flatten the G1 and G2 bands and position them at intensity values of 1 and 2, respectively.
 - The stripe is computed relative to the object intensities (thus, it should be multiplied or divided).


* `generateFunction` - Create a simple function that operates on an image based on a string expression.
```Matlab
fun = t.generateFunction(expression,removeUndefined,expandInput,requiredVars)
% Ex.
% expression = "(S - BG_o*BG_s)/(FG_f*FG_s)"
% fun = t.generateFunction(expression,false,true,["S"])
% fun
%   Function handle
%   @(S,x,y) (S - BG_o(x,y) .* BG_s(x)) / (FG_f(x,y) .* FB_s(x))
```
 - Allowed variables are `BG_o`, `BG_s`, `FG_o`, `FG_f`, `FG_s`, and any required variables.


* `generateCorrectionFunction` - Generate the correction function that corresponds to the `Image_Correction_Expression` property.


* `clearCorrections` - Delete all current corrections
```Matlab
t.BG_offset = [];
t.BG_stripeX = [];
t.FG_offset = [];
t.FG_factor = [];
t.FG_stripeX = [];
```

* `plot` - Plot one of the computed corrections
```Matlab
t.plot(name)
```
 - Name can be one of the following strings: `threshold`, `background`, `foreground_factor`, `foreground_offset`, `sharpness`, `stripe_background`, `stripe_foreground`

### List of low level methods
* `isEmpty`
* `delete`
* `open`
* `close`
* `getColumn`
* `getRow`
* `getTile`
* `getBlock`
* `computeTileNumbers`
* `generateCorrectionFunction`
* `smoothSurf`
