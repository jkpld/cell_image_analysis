# Partitioning objects
These classes are for partitioning nuclei that are touching each other.

## `ObjectPartitioner` abstract class
An abstract class from which different object partitioners can be created.

#### Properties
* `Minimum_Object_Size` - Remove all objects smaller than this (units in pixels)
* `Minimum_Hole_Size` - Fill in all holes small than this.
* `Use_Parallel` - Use parallel computing (multiple CPU's) to speed up

The remaining properties of this class are used to reduce the number of objects we try to partition. For example, using the object areas and the DAPI intensities, we can remove junk objects, blurry objects, or objects that likely contains more than some larger number of nuclei.

* `Area_Normalizer` - Function that takes in image locations and outputs the median object area at those positions. The areas are normalized before going through the `attemptPartitioning()` functions.

* `Intensity_Normalizer` - Function that takes in image locations and outputs the median object intensity at those positions. The intensities are normalized before going through the `attemptPartitioning()` method.

* `Restricted_Partitioning` - Logical flag. If true, then the object area and intensitiy will be normalized (using the `Area_Normalizer` and `Intensity_Normalizer`) and the results passed to the `attemptPartitioning()` method. If false, then all objects will be kept, and all partitioning will be attempted on all objects.

* `Acceptance_Contour` - Nx2 array describing a contour in the normalized DAPI-Area plane. Objects inside the contour will attempt to be partitioned. Objects outside of the contour will be ignored. Default value is (Normalized DAPI Intensity, Normalized Area) `[3, 1.1; 10, 4; 10, 7; 5, 7; 2, 3; 1.6, 2.5; 1.6, 1.1; 3, 1.1]`

#### Methods
* `attemptPartitioning` - Return true for all objects that we should attempt to partition. If the normalized object area and intensities are inside of the Acceptance_Contour, then we will try and partition the objects.
  ```Matlab
  tf = attemptPartitioning(obj, A, I)
  ```

* `preProcess` - Determine the objects to try and partition and remove small holes from the objects.
  ```Matlab
  [CC, attempt] = preProcess(obj, I, BW, Image_Offset)
  ```

* `postProcess` - Remove boundary objects and objects smaller than the minimum allowed size.
  ```Matlab
  BW = postProcess(obj, BW)
  ```

* `partition` - (abstract) Partition objects.
  ```Matlab
  mask_partioned = partition(obj, BW, I, Image_Offset)
  ```

## Concrete `ObjectPartitioner` classes

### `None`
Perform no partitioning.

### `DTWatershed`
Use the watershed method on the smoothed distance transform for partitioning.

### `SALRGeoPartition`
Use SALR clustering to locate the center of each nuclei, and then use a seed-point based geometric partitioning method to partition. This requires that both the SALR Clustering and the geometricPartitioing methods have been installed. (These are two other repositories on my github.)
