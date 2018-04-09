# Feature extraction
These classes are used for extracting features from each channel for each object.

## `FeatureExtractor` class
This class holds a list of the different `featureGroups` that will be extracted from each object in an image and it contains the method for extracting all of the features.

The class has two properties.
* `channels` - An array of the channel names. This is automatically set to the `Channel_Names` when the `Compute_Object_Features` method of `CellExperiment` is used.

* `featureGroups` - An array of `FeatureGroup` classes. These are the features that will be extracted from the objects.

The class methods are

* `Compute` - Compute the features for objects in an image
  ```Matlab
  [x, names] = Compute(obj, mask, I, Image_Offset)
  % x : features (NxM) where N is number of objects in MASK and M is the number of features
  % names: array of strings giving the feature names
  % mask : object mask
  % I : image Ny x Nx x Nc, where Nc is the number of channels
  % Image_Offset : the location of the upper right pixel of the image I. Use this when extracting
  % features (locations) from a patch of a large image.
  ```

* `addFeatureGroup` - Add a new `FeatureGroup` to the `featureGroup` list

## `FeatureGroup` abstract class
An abstract class defining from which the different types of feature groups will be based on.

Properties:
* `Channel` - Name of the channel from which to extract the features from
* `Options` - (struct) Options for feature extraction
* `GroupName` - (protected) Name of the FeatureGroup
* `requriedOptions` - (cell, protected) list of options required for feature group
* `FeatureNames` - (array, protected, abstract) List of feature names computed by group

Methods:
* `Compute` - (abstract) Compute features. , .
  ```Matlab
  x = Compute(obj, I, L)
  % I is the image
  % L is the label matrix
  ```
  - If `Channel` is empty, then the mask is used for computing the features. This is only true for the `Location` and `Shape` `FeatureGroups`

### Concrete `FeatureGroup` classes
See each class for full documentation

#### `Location`
Compute the centroids and bounding boxes of each object.

#### `Shape`
Compute
* Area - the object area
* Major axis length - the major axis length of an ellipse fit the object
* Minor axis length - the minor axis length of an ellipse fit to the object
* Eccentricity - the eccentricity of the ellipse fit to the object
* Solidity - the area of the object divided by the area of the object's convex hull
* Form factor - the area of the object divided by the area of a circle with the same circumference as the object's perimeter
* Fourier descriptors - translation, rotation, and scale invariant Fourier descriptors that describe the object's shape

#### `Intensity`
Compute the sum, mean, std, skewness, and kurtosis of the object intensity

#### `RadialIntensity`
Compute the mean and std of the intensity from concentric rings around the object and the center region. For example, see the image below. Here we have three rings and a center region. From each of these four regions, we compute the mean and std of the intensity.


![Radial_Intensity](/docs/radialIntensity.PNG)
#### `HaralickTexture`
Compute Haralick texture features (averaged over the 0, 45, 90, 135 degree directions). This computes the same 13 Haralick texture features used by CellProfiler.

#### `Granularity`
Compute the granularity spectrum of each object.

#### `BasicProps`
A feature group that only extracts the centroid location, area, and object integrated (sum) intensity.

#### `EmptyFeature`
Empty (null) feature group.
