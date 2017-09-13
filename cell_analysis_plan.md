# Cell analysis plan

## DAPI Image
1. Compute the image threshold
1. Compute smoothed background and image blur
1. Compute background stripe artifact
1. Create total background and offset foreground function
 * Save this total background
 * Also save the smoothed background, mean background, and background stripe.
1. Re-compute the threshold after applying the background correction
 * Save this threshold
1. Compute smoothed foreground surface (used for computing FG stripe artifact)
1. Compute foreground stripe artifact
1. Extract initial nuclei info (centroid, area, integrated intensity)
 * Read in block
 * Create object mask
     - Apply saved background correction
     - Apply saved image threshold
 * Subtract away foreground stripe artifact
 * Extract the initial features for each object
 * Check to make sure the strip artifact does not show up when the nuclei data is plotted.
1. Create smooth surface giving the median area and integrated intensity.
 * This should represent the G1 nuclei band.
 * This is used to normalize the area and integrated intensity later on when determining which objects should be kept and which objects we should attempt to partition.
1. Partition the objects and save mask
 * Read in block
 * Create object mask
     - Apply saved background correction
     - Apply saved image threshold
 * Extract initial features for each object
     - First subtract away foreground stripe artifact
     - _Note. These features are the same as were previously extracted, but it could be messy to try to match the features to the object number from each block of the image; therefore, even though it takes more time, I am re-extracting the object initial features._
 * Determine of what objects to keep and what objects to try to partition
 * Partition objects in need.
 * Write block of image mask

### We now have a binary mask for all objects in the image!

## Additional image channels
1. Read in object mask from DAPI
2. Use the dilated inverse mask to compute the smoothed background surface for each channel
3. Compute background stripe artifact for each channel
4. Create total background function for each channel
 * Save these backgrounds

## Eextract features
1. Read in object mask from DAPI
2. Extract object features from each channel. _This could be done at a later time using the saved image mask and the saved background correction._
