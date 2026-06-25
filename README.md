## Agrawal Cluster Disassembly

This plugin takes 2 channel (green and red) .lif or .tif timelapse images and first segments the pupae ROIS in the green
channel. Per pupa the cluster are is found in the red channel. The cluster area is further segmented in the green channel
to give areas of positive and an area of negative signal. 

### Method

- Pupae segmentation is performed manually on a maximum time projection of the image, this is done for all images in the
folder before analysis is initiated.
- A reference timepoint is requested for each file, this is the timepoint against which area and intensity are normalised.
- Cluster segmentation is then performed by applying a rolling ball (50 pixel) background subtraction to the red channel
and then a Moments threshold is applied. Clusters are filtered based on size with a user input Minimum area (we use 50 
pixels). 
- The green signal outline is thresholded based on the max-projection and a Default threshold on the green channel.
- Areas within the clusters that have green marker signal are next segmented using ImageJ thresholding 'RenyiEntropy' 
performed on the entire stack with the outline applied. The outline makes this thresholding more robust to variation in
background levels of green signal. Only ROIs that are greater than 10 pixels in area and within a detected cluster are 
selected.
- Statistics are then calculated per pupa and per timepoint along with the thresholds used for cluster and marker segmentation.
These are;
  - Cluster Area
  - Marker Area
  - Normalised Cluster Area
  - Normalised Marker Area
  - Mean Cluster Intensity (red)
  - Mean Cluster Intensity (green)
  - Total Cluster Intensity (red)
  - Total Cluster Intensity (green)  
  - Normalised Cluster Intensity (red)
  - Normalised Cluster Intensity (green)
  - Maximum Cluster Intensity (red)
  - Maximum Cluster Intensity (green)
  - Mean Marker Intensity (red)
  - Mean Marker Intensity (green)
  - Total Marker Intensity (red)
  - Total Marker Intensity (green)
  - Normalised Marker Intensity (red)
  - Normalised Marker Intensity (green)
  - Area Ratio (marker/cluster)
  - Change in Cluster Area since previous timepoint
  - Change in Marker Area since previous timepoint

### Output

The plugin indexes the pupae left to right(A-Z) and outputs a spreadsheet for the folder with each line detailing the 
image name,timepoint, pupal number and the statistics given above.

An Output image displaying both channels and an additional channel showing the detected, labeled and measured ROIs is output 
for each image stack.
