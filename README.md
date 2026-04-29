## Agrawal Cluster Disassembly

This plugin takes 2 channel (green and red) .lif or .tif timelapse images and first segments the pupae ROIS in the green
channel. Per pupa the cluster are is found in the red channel. The cluster area is further segmented in the green channel
to give areas of positive and an area of negative signal. 

### Method

- Pupae segmentation is performed using cellpose (user specified environment and model we use version 3.10, and the cyto 
3 model) to segment a time-projected image of the green channel.
- Cluster segmentation is then performed by applying a Labkit classifier to the red channel full stack. Clusters are filtered
based on feret ratio and distance such that cluster with a high feret ratio (>2) that are less than 100 pixels from the pupa
edge are excluded.
- Areas within the clusters that have green marker signal are next segmented using imageJ thresholding "RenyiEntropy" performed
on the entire stack. Only ROIs that are greater than 10 pixels in area and within a detected cluster are selected.
- Statistics are then calculated per pupa and per timepoint. These are;
  - Cluster Area
  - Marker Area
  - Mean Cluster Intensity (red)
  - Mean Cluster Intensity (green)
  - Mean Marker Intensity (red)
  - Mean Marker Intensity (green)

### Output

The plugin indexes the pupae left to right(A-Z) and outputs a spreadsheet for the folder with each line detailing the 
image name,timepoint, pupal number and the statistics given above.

An Output image displaying both channels and an additional channel showing the detected, labelled and measured ROIs is output 
for each image stack.
