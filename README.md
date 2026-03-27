## Agrawal Pupae Quantification

This plugin takes a mixed folder of .lif and .itf files as input. For each file it extracts each image sequentially, finds the number 
of channels and determines which is green, white or red (if there are 3 or more channels).  

A segmentation is performed on the white (Transmitted Light) channel to find the Pupae outlines. This is done first using
a Huang automatic threshold, the resulting masks are then split using the 'Watershed' function. The ROIs created from this
splitting are then recombined vertically to create the final Pupae masks.

For each pupae cell clusters are then detected in the green channel using the Triangle automatic threshold, the resulting 
masks are filtered by size (500-Infinity pixels) and circularity (0.05-1.00).

If three channels are present the red channel is then segmented to detect the nuclei. First a 25 pixel rolling ball background
subtraction is performed, this is followed by thresholding with the Yen automatic thresholder. ROIs are then selected, again
filtered by size (10-500 pixels) and circularity (0.05-1.00).

### Edit

There is the option to edit detected pupae and cluster ROIs. If this option is selected then after each segmentation the 
green channel will be displayed with first the Pupae ROIs in the ROI manager and then the cluster ROIs. The ROis can then
be edited (added/deleted/redrawn) using the built in IJ functionality. Once the user is happy with the ROIs displayed then
the plugin will analyse and report the results of the ROis before moving on to the next image in the folder.

### Output

The plugin indexes the pupae left to right(A-Z), and the clusters top to bottom (1-10) and outputs a spreadsheet for the
folder with each line detailing the image name, pupae letter, pupae area, the cluster number and area for each detected cluster.
If there are 3 channels, the total area and number of nuceli found for each cluster will be output to a second results file. 

An Output image displaying all 2/3 channels and an additional channel showing the detected, labelled and measured ROIs is output 
for each image stack.
