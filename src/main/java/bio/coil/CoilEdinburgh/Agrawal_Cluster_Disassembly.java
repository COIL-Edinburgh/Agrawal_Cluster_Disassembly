/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package bio.coil.CoilEdinburgh;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.ChannelSplitter;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.ZProjector;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.LUT;
import io.scif.services.DatasetIOService;
import io.scif.services.FormatService;

import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imagej.roi.ROIService;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.IOException;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;


/**
 * This example illustrates how to create an ImageJ {@link Command} plugin.
 * <p>
 * </p>
 * <p>
 * You should replace the parameter fields with your own inputs and outputs,
 * and replace the {@link run} method implementation with your own logic.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins>Users Plugins>Agrawal Cluster Disassembly")
public class Agrawal_Cluster_Disassembly<T extends RealType<T>> implements Command {
    //
    // Feel free to add more parameters here...
    //
    @Parameter
    private FormatService formatService;

    @Parameter
    private DatasetIOService datasetIOService;

    @Parameter
    private UIService uiService;

    @Parameter
    private OpService ops;

    @Parameter
    private ROIService roiService;

    @Parameter(label = "Open Folder: ", style = "directory")
    public File filePath;

    @Parameter(label = "Min Cluster Area: ")
    public int minArea;

    @Parameter(label = "Red Channel: ")
    public int red = 1;

    @Parameter(label = "Green Channel: ")
    public int green = 0;

    RoiManager roiManager;
    double pixelSize;
    int nChannels;
    String lifFileName;
    String filename;
    ArrayList<ArrayList<Roi[]>> clusters = new ArrayList<>();
    ArrayList<ArrayList<Roi[]>> greenAreas = new ArrayList<>();
    ArrayList<Roi[]> pupae = new ArrayList<>();
    ArrayList<ArrayList<double[]>> stats = new ArrayList<>();
    ArrayList<String> letters;
    ArrayList<double[]> redThreshold = new ArrayList<>();
    ArrayList<double[]> greenThreshold = new ArrayList<>();
    int[] refFrame;
    File[] files;

    @Override
    public void run() {

        files = filePath.listFiles();
        refFrame = new int[files.length];

        if (RoiManager.getInstance() != null) {
            roiManager = RoiManager.getInstance();
        } else {
            roiManager = new RoiManager();
        }
        roiManager.reset();

        //open folder
        assert files != null;

        //Open each file in the folder and allow the user to draw the Pupae ROIS
        for (int i=0; i<files.length; i++) {
            File file = files[i];
            ImagePlus outputImp = null;
            if ((file.toString().contains(".tif")||file.toString().contains(".lif")) &&
                    !(file.toString().contains("Output")||file.toString().contains(".roi")||file.toString().contains(".zip"))) {
                IJ.run("Bio-Formats Importer", "open=[" + file.getAbsolutePath() + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
                ImagePlus imp = WindowManager.getCurrentImage();
                Roi[] tempPupae = drawPupae(imp);
                savePupae(tempPupae, file);
                setRefFrame(i);
                imp.close();
            }
        }

        //For each file in the folder use the saved Rois to run the analysis and create the output
        for (int i=0; i<files.length; i++) {
            File file = files[i];
            ImagePlus outputImp = null;
            if ((file.toString().contains(".tif")||file.toString().contains(".lif")) &&
                    !(file.toString().contains("Output")||file.toString().contains(".roi")||file.toString().contains(".zip"))) {
                IJ.run("Bio-Formats Importer", "open=[" + file.getAbsolutePath() + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
                ImagePlus imp = WindowManager.getCurrentImage();
                try {
                    runAll(imp, outputImp, i);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }

    //Allows the user to draw ROI on a max-projection of the timeseries and add them to the ROI manager
    private Roi[] drawPupae(ImagePlus imp){
        WaitForUserDialog drawPupae = new WaitForUserDialog("Draw ROIs around pupae", "Draw ROIs around pupae and press 'T' to add to ROI manager");
        ImagePlus impZ  = ZProjector.run(imp,"Max");
        impZ.show();
        drawPupae.show();
        Roi[] pupaeRois = roiManager.getRoisAsArray();
        roiManager.reset();
        leftToRight(pupaeRois);//Sorts pupae left to right
        roiManager.reset();
        impZ.close();
        return pupaeRois;
    }

    //Sorts Rois left to right
    private void leftToRight(Roi[] rois) {
        Roi temp = null;
        for (int i = 0; i < rois.length; i++) {
            // inner loop to compare and swap elements
            for (int j = i + 1; j < rois.length; j++) {
                if (rois[i].getBounds().x > rois[j].getBounds().x) {
                    temp = rois[i];
                    rois[i] = rois[j];
                    rois[j] = temp;
                }
            }
        }
    }

    //Saves the ROIs in the ROI manager so they can be used later
    private void savePupae(Roi[] pupae, File file){
        roiManager.reset();
        for (Roi points : pupae) {
            roiManager.addRoi(points);
        }
        String pupaeName = file.getAbsolutePath();
        String fileExt = ".tif";
        String newExt = "_roi.zip";
        if(pupaeName.contains(".lif")){fileExt=".lif";}
        if (pupae.length<2){newExt = ".roi";}
        pupaeName = pupaeName.replace(fileExt,newExt);
        roiManager.save(pupaeName);
        roiManager.reset();
    }

    private void setRefFrame(int i){
        NonBlockingGenericDialog output = new NonBlockingGenericDialog( "Select the reference frame:");
        output.addNumericField( "Reference Frame: ", 0);
        output.showDialog();
        double refframe = output.getNextNumber();
        refFrame[i]= (int) refframe;
    }

    //Initializes the ROI and Stats arraylists and then runs the segmentation, gets the stats and outputs the results
    // file and overlay image.
    private void runAll(ImagePlus imp, ImagePlus outputImp, int i) throws IOException {
        File file = files[i];
        clusters = new ArrayList<>();
        pupae = new ArrayList<>();
        greenAreas = new ArrayList<>();
        stats = new ArrayList<>();

        runSegmentation(imp, file);
        IJ.log("Segmentation Done. Starting getStats()");

        getStats(imp);
        IJ.log("Stats Done. Creating Results file");

        createResultsFile(String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Results.csv")), i);
        IJ.log("Results Saved. Creating output Imp");

        outputImp = drawOverlay(imp, outputImp);
        IJ.log(String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Output.tif")));
        IJ.save(outputImp, String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Output.tif")));
        IJ.run("Close All", " ");
    }

    //Segments the clusters (red channel) and markers (green channel) for each time point and fills the cluster. pupae
    // and marker ArrayLists with one ROI[] or Arraylist<Roi[]> for each timepoint.
    private void runSegmentation(ImagePlus imp, File file) {

        //Gets the filename without extensions
        IJ.log(imp.getTitle());
        filename = imp.getTitle();
        filename = filename.replace(".tif","");
        filename = filename.replace(".lif","");
        lifFileName = file.getName();

        //Gets the channel order
        ImagePlus[] split = ChannelSplitter.split(imp);
        if (!file.toString().contains(".tif")) {
            getChannelOrder(imp);//for .lif files
        } else if (file.toString().contains(".tif")) {
            getChannelOrderTif(imp);//for .tif files
        }

        //Opens the saved pupae ROI files
        Roi[] pupa_all = openPupae(file);
        imp.show();
        //Gets an outline of the pupae in the green channel (used later to make the marker thresholding less influenced
        // by number of pupae/ amount of background pixels)
        Roi outline = getOutline(split[green]);

        ImagePlus clusterMasks = split[red];

        //for each timepoint
        for(int i =0; i< imp.getNFrames();i++) {

            pupae.add(pupa_all);
            ArrayList<Roi[]> pupae_clusters = new ArrayList<>();
            ArrayList<Roi[]> greenAreaList = new ArrayList<>();
            redThreshold.add(new double[pupae.get(i).length]);
            greenThreshold.add(new double[pupae.get(i).length]);
            //for each pupae
            for (int j =0; j< pupae.get(i).length; j++) {
                //Segment Clusters (red channel)
                Roi pupa = pupae.get(i)[j];
                clusterMasks.show();
                Roi[] cluster = getClusters(clusterMasks, pupa, i, j);
                pupae_clusters.add(cluster);
                //for each cluster get the contained marker ROIs (green channel)
                Roi[][] greenArea = getMarkers(cluster, split[green], i, outline, j);
                greenAreaList.add(mergeArrays(greenArea));//merges all marker regions within one pupa to a single array
            }
            clusters.add(pupae_clusters); //update the Arraylists for this timepoint
            greenAreas.add(greenAreaList);
            IJ.log("Frame: "+ i);//report progress - useful for finding bugs
        }
    }

    //Sets the index of the red and green channels from image metadata (.lif files)
    private void getChannelOrder(ImagePlus imp) {
        LUT[] luts = imp.getLuts();
        for (int i = 0; i < luts.length; i++) {
            String string = luts[i].toString();
            if (string.contains("green")) {
                green = i;
            } else if (string.contains("red")) {
                red = i;
            }
        }
    }

    //Sets the index of the red and green channels from image metadata (.tif files)
    private void getChannelOrderTif(ImagePlus imp) {
        String info = imp.getInfoProperty();
        for (int i = 0; i < nChannels; i++) {
            if(info.contains("Image #"+i+"|ChannelDescription #"+i+"|LUTName = Green")){green=i;}
            if(info.contains("Image #"+i+"|ChannelDescription #"+i+"|LUTName = Red")){red=i;}
        }
    }

    //Opens previously saved user drawn pupae ROIs
    private Roi[] openPupae(File file){
        String roiFile = file.getAbsolutePath();
        String fileExt = ".tif";
        if(roiFile.contains(".lif")){fileExt=".lif";}
        roiFile = roiFile.replace(fileExt,"_roi.zip");
        if(new File(roiFile).exists()){
            roiManager.open(roiFile);}
        else{
            roiFile = roiFile.replace("_roi.zip", ".roi");
            if (new File(roiFile).exists()){
                roiManager.open(roiFile);}
        }
        return roiManager.getRoisAsArray();
    }

    //Get the outline ROI of the green channel from a Max-Projection.
    private Roi getOutline(ImagePlus imp){
        ImagePlus impZ  = ZProjector.run(imp,"Max");
        impZ.show();
        IJ.setAutoThreshold(impZ, "Default dark");
        IJ.run(impZ, "Analyze Particles...", "size=10-Infinity pixel circularity=0.0-1.00 add");
        Roi[] pupae = roiManager.getRoisAsArray();
        impZ.close();
        roiManager.reset();
        return pupae[0];
    }

    //For a single timepoint and pupa, subtracts background and then sets the pupa ROI before thresholding (Moments).
    //Returns an Roi[] of clusters for that pupa.
    private Roi[] getClusters(ImagePlus imp, Roi pupa, int i, int j) {
        imp.show();
        imp.setT(i+1);
        IJ.run(imp, "Subtract Background...", "rolling=50 slice");
        imp.setRoi(pupa);
        IJ.setAutoThreshold(imp, "Moments dark");
        redThreshold.get(i)[j] = IJ.getProcessor().getAutoThreshold();
        IJ.run(imp, "Analyze Particles...", "size=" + minArea + "-Infinity pixel circularity=0.0-1.00 add");
        Roi[] allClusters = roiManager.getRoisAsArray();
        if (allClusters.length == 0) {
            return allClusters;
        }
        roiManager.reset();
        return allClusters;
    }

    //For a single timepoint and clusters ROI[] (one of these per pupa), segments the region of maker signal within the
    // clusters. Thresholds within the outline Roi on the whole stack and then for each cluster Roi analyses particles to
    //return thresholded marker within that cluster.
    private Roi[][] getMarkers(Roi[] clusters, ImagePlus imp, int frame, Roi outline, int j) {

        imp.show();
        imp.setT(frame+1);
        Roi[][] output = new Roi[clusters.length][];
        imp.setRoi(outline);
        IJ.setAutoThreshold(imp, "RenyiEntropy dark stack");
        greenThreshold.get(frame)[j] = IJ.getProcessor().getAutoThreshold();
        for (int i = 0; i < clusters.length; i++) {
            imp.setRoi(clusters[i]);
            IJ.run(imp, "Analyze Particles...", "size=10-Infinity pixel circularity=0.0-1.00 add");
            output[i] = roiManager.getRoisAsArray();
            roiManager.reset();
        }
        return output;
    }

    //merges an array of arrays into one array
    private Roi[] mergeArrays(Roi[][] arrays){
        ArrayList<Roi> outputlist = new ArrayList<>();
        for(Roi[] array:arrays){
            outputlist.addAll(Arrays.asList(array));
        }
        Roi[] output = new Roi[outputlist.size()];
        for (int i = 0; i< output.length; i++){
            output[i] = outputlist.get(i);
        }
        return output;
    }

    //Uses the Cluster and Marker Roi to measure areas and intensities, populates the stats ArrayList<ArrayList<double[]>>
    // with an arraylist of timepoints, each containing a further arraylist of pupae. For each pupae is a double[] of
    //measured statistics.
    private void getStats(ImagePlus imp) {
        ImagePlus[] split = ChannelSplitter.split(imp);

        for (int i = 0; i < imp.getNFrames(); i++) { //for each timepoint
            imp.setT(i+1);
            ArrayList<double[]> frameStats = new ArrayList<>();//for each timepoint create a new ArrayList<>

            for (int j=0; j<pupae.get(i).length;j++){ //for each pupa
                double[] pupaeStats = new double[8]; //create a new double[]

                for (int k = 0; k<clusters.get(i).get(j).length; k++){ //for each cluster
                    if(clusters.get(i).get(j).length>0) {
                        //Add Cluster Area
                        pupaeStats[0] = pupaeStats[0] + clusters.get(i).get(j)[k].getStatistics().area;

                        //Add total Cluster intensity in green channel
                        ImageProcessor ipg = split[green].getChannelProcessor();
                        ipg.setRoi(clusters.get(i).get(j)[k]);
                        pupaeStats[1] = pupaeStats[1] + ipg.getStatistics().mean * ipg.getStatistics().area;

                        //Add total Cluster intensity in red channel
                        ImageProcessor ipr = split[red].getChannelProcessor();
                        ipr.setRoi(clusters.get(i).get(j)[k]);
                        pupaeStats[2] = pupaeStats[2] + ipr.getStatistics().mean * ipr.getStatistics().area;
                        if(ipg.getStatistics().max > pupaeStats[7]){
                            pupaeStats[7] = ipg.getStatistics().max;
                        }
                        if(ipr.getStatistics().max > pupaeStats[6]){
                            pupaeStats[6] = ipr.getStatistics().max;
                        }
                    }
                }
                for (int m = 0; m<greenAreas.get(i).get(j).length; m++){ //for each marker region
                    if(greenAreas.get(i).get(j).length>0) {
                        //Add Marker Area
                        pupaeStats[3] = pupaeStats[3] + greenAreas.get(i).get(j)[m].getStatistics().area;

                        //Add total marker intensity in the green channel
                        ImageProcessor ipg = split[green].getChannelProcessor();
                        ipg.setRoi(greenAreas.get(i).get(j)[m]);
                        pupaeStats[4] = pupaeStats[4] + ipg.getStatistics().mean * ipg.getStatistics().area;

                        //Add total marker intensity in the red channel
                        ImageProcessor ipr = imp.getChannelProcessor();
                        ipr.setRoi(greenAreas.get(i).get(j)[m]);
                        pupaeStats[5] = pupaeStats[5] + ipr.getStatistics().mean * ipr.getStatistics().area;
                    }
                }
                frameStats.add(pupaeStats);
            }
            IJ.log("Stats: "+ i);
            stats.add(frameStats);
        }
    }

    //Create the Results .csv
    public void createResultsFile(String name, int i) throws IOException {
        Date date = new Date(); // This object contains the current date value
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy, hh:mm:ss");
        boolean exists = new File(name).exists();
        try {
            FileWriter fileWriter = new FileWriter(name, true);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            //If this is a new file
            if(!exists) {
                bufferedWriter.newLine();
                bufferedWriter.write(formatter.format(date));
                bufferedWriter.newLine();
                bufferedWriter.write("Min Cluster Size: " + minArea); //User input parameter
                bufferedWriter.newLine();
                bufferedWriter.write("Reference Frame: " + refFrame[i]); //User input parameter
                bufferedWriter.newLine();
                String heading = "";

                //Create column headings
                for (int n = 0; n< pupae.get(0).length; n++) {
                    heading = heading + "Timepoint, Red Threshold, Green Threshold, Pupal No,Cluster Area,Marker Area," +
                            "Norm Cluster Area,Norm Marker Area," +
                            "Mean Cluster Intensity (red),Mean Cluster Intensity (green)," +
                            "Total Cluster Intensity (red),Total Cluster Intensity (green)," +
                            "Norm Cluster Intensity (red),Norm Cluster Intensity (green)," +
                            "Max Cluster Intensity (red),Max Cluster Intensity (green)," +
                            "Mean Marker Intensity (red),Mean Marker Intensity (green), " +
                            "Total Marker Intensity (red),Total Marker Intensity (green), " +
                            "Norm Marker Intensity (red),Norm Marker Intensity (green)," +
                            "Area Ratio, Change in area (cluster), Change in area (marker), ,";
                }
                bufferedWriter.write(heading);//write header 1
                bufferedWriter.newLine();
            }

            //For each timepoint
            for (int t = 0; t < stats.size(); t++) {
                StringBuilder data = new StringBuilder();
                //For each pupa
                for (int j = 0; j < stats.get(t).size(); j++) {
                    data.append(t); //timepoint
                    data.append(",").append(redThreshold.get(t)[j]);
                    data.append(",").append(greenThreshold.get(t)[j]);
                    data.append(",").append(j); //pupae number
                    double[] results = stats.get(t).get(j);
                    double[] results_t1 =  stats.get(t).get(j);
                    if(t>0){ results_t1 =  stats.get(t-1).get(j);}
                    data.append(",").append(results[0]); //cluster area
                    data.append(",").append(results[3]); //marker area
                    data.append(",").append(results[0]/stats.get(refFrame[i]).get(j)[0]); //norm cluster area
                    data.append(",").append(results[3]/stats.get(refFrame[i]).get(j)[3]); //norm marker area
                    data.append(",").append(results[1]/results[0]); //mean cluster intensity (red)
                    data.append(",").append(results[2]/results[0]); //mean cluster intensity (green)
                    data.append(",").append(results[1]); //total cluster intensity (red)
                    data.append(",").append(results[2]); //total cluster intensity (green)
                    data.append(",").append(results[1]/stats.get(refFrame[i]).get(j)[1]); //norm cluster intensity (red)
                    data.append(",").append(results[2]/stats.get(refFrame[i]).get(j)[2]); //norm cluster intensity (green)
                    data.append(",").append(results[6]); //max cluster intensity (red)
                    data.append(",").append(results[7]); //max cluster intensity (green)
                    data.append(",").append(results[4]/results[3]); //mean marker intensity (red)
                    data.append(",").append(results[5]/results[3]); //mean marker intensity (green)
                    data.append(",").append(results[4]); //total marker intensity (red)
                    data.append(",").append(results[5]); //total marker intensity (green)
                    data.append(",").append(results[4]/stats.get(refFrame[i]).get(j)[4]); //norm marker intensity (red)
                    data.append(",").append(results[5]/stats.get(refFrame[i]).get(j)[5]); //norm marker intensity (green)
                    data.append(",").append(results[3]/results[0]); //area ratio (marker/cluster)
                    data.append(",").append(results[0]-results_t1[0]); //Change in area (cluster)
                    data.append(",").append(results[3]-results_t1[3]); //Change in area (marker)
                    data.append(", ,");

                }
                bufferedWriter.write(String.valueOf(data));
                bufferedWriter.newLine();
            }
            bufferedWriter.close();
            IJ.log("Results file updated");
        } catch (IOException ex) {
            System.out.println(
                    "Error writing to file '" + name + "'");
        }
    }

    //Creates an ImagePlus with the drawn and thresholded Rois drawn and labelled and merged with the original
    // timeseries to give a 3 Channel Output image.
    private ImagePlus drawOverlay(ImagePlus imp, ImagePlus outputImp){

        //For each Frame
        for(int i=0; i< imp.getNFrames();i++) {
            //Create a single frame ImagePlus of the timepoint (i)
            ImagePlus impFrame = new Duplicator().run(imp, 1, 2, 1, 1, i+1, i+1);
            //Draw on the ROIs and labels and merge channels
            ImagePlus overlay = createOutputImage(impFrame, pupae.get(i), clusters.get(i), greenAreas.get(i));
            IJ.log("Overlay: "+ i);
            impFrame.close();

            //Add frame to output image
            if (outputImp == null) {
                outputImp = overlay;
            } else {
                outputImp = Concatenator.run(outputImp, overlay);
            }
        }
        return outputImp;
    }

    //Draws Rois and Labels on blank ImagePlus and then merges with the red and green channels of the original image to
    //create a 3 channel image.
    private ImagePlus createOutputImage(ImagePlus imp, Roi[] pupae, ArrayList<Roi[]> clusters, ArrayList<Roi[]> greenArea) {

        //Create blank image
        ImagePlus masks = IJ.createImage("Output_Masks", "16-bit", imp.getWidth(), imp.getHeight(), 1, 1, 1);
        masks.show();
        ImageProcessor ip = masks.getProcessor();
        Font font = new Font("SansSerif", Font.BOLD, 40);
        ip.setFont(font);
        ip.setColor(Color.getHSBColor(0, 0, 100));
        letters = new ArrayList<>();
        for (char letter = 'A'; letter <= 'Z'; letter++) {
            letters.add(Character.toString(letter));
        }

        ip.drawString(imp.getTitle(), 50, 50); //Write filename

        for (int i = 0; i < pupae.length; i++) { //Draw pupa and labels
            ip.drawString(letters.get(i), pupae[i].getBounds().x + 150, pupae[i].getBounds().y + 150);
            ip.draw(pupae[i]);
        }

        font = new Font("SansSerif", Font.BOLD, 30);
        ip.setFont(font);
        for (Roi[] roiArray : clusters) { //Draw clusters
            for (Roi points : roiArray) {
                ip.draw(points);
            }
        }

        ip.setLineWidth(4); //Change linewidth
        for (Roi[] roiArray : greenArea) { //Draw marker regions
            for (Roi points : roiArray) {
                ip.draw(points);
            }
        }

        masks.updateAndDraw();
        masks.show();

        //Split channels of original image and display so that the 'Marge Channels...' command works
        ImagePlus[] split = ChannelSplitter.split(imp);
        split[green].show();
        split[red].show();
        IJ.run("Merge Channels...", "c1=[" + split[red].getTitle() + "] c2=[" + split[green].getTitle() +
                    "] c5=[" + masks.getTitle() + "] create");
        return WindowManager.getCurrentImage();
    }




    /**
     * This main function serves for development purposes.
     * It allows you to run the plugin immediately out of
     * your integrated development environment (IDE).
     *
     * @param args whatever, it's ignored
     * @throws Exception
     */
    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        ij.command().run(Agrawal_Cluster_Disassembly.class, true);
    }

}
