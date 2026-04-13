## Agrawal Cluster Disassembly

This plugin takes 2 channel (green and red) .lif or .tif timelapse images and first segments the pupae ROIS in the green
channel. Per pupa the cluster are is found in the red channel. The cluster area is  further segmented in the green channel
to give areas of positive and an area of negative signal. 

### Output

The plugin indexes the pupae left to right(A-Z) and outputs a spreadsheet for the folder with each line detailing the 
image name, pupae letter, total cluster area and green positive area over the course of the timelapse.

An Output image displaying both channels and an additional channel showing the detected, labelled and measured ROIs is output 
for each image stack.
