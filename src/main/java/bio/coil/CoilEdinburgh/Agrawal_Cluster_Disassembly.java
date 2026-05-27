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
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.ChannelSplitter;
import ij.plugin.Concatenator;
import ij.plugin.Duplicator;
import ij.plugin.ZProjector;
import ij.plugin.frame.RoiManager;
import ij.process.FloatPolygon;
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

    //@Parameter(label = "Open Env Path: ", style = "directory")
    //public File envpath;

    //@Parameter(label = "Open Model: ", style = "file")
    //public File modelpath;

    //@Parameter(label = "Open Classifier: ", style = "file")
    //public File classifier;

    @Parameter(label = "Min Cluster Area: ")
    public int minArea;

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
    int green = 0;
    int red = 1;

    @Override
    public void run() {

        File[] files = filePath.listFiles();

        if (RoiManager.getInstance() != null) {
            roiManager = RoiManager.getInstance();
        } else {
            roiManager = new RoiManager();
        }
        roiManager.reset();

        //open folder
        assert files != null;

        for (File file : files) {
            ImagePlus outputImp = null;
            if ((file.toString().contains(".tif")||file.toString().contains(".lif")) &&
                    !(file.toString().contains("Output")||file.toString().contains(".roi")||file.toString().contains(".zip"))) {
                IJ.run("Bio-Formats Importer", "open=[" + file.getAbsolutePath() + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
                ImagePlus imp = WindowManager.getCurrentImage();
                Roi[] tempPupae = drawPupae(imp);
                savePupae(tempPupae, file);
                imp.close();
            }
        }

        for (File file : files) {
            ImagePlus outputImp = null;
            if ((file.toString().contains(".tif")||file.toString().contains(".lif")) &&
                    !(file.toString().contains("Output")||file.toString().contains(".roi")||file.toString().contains(".zip"))) {
                IJ.run("Bio-Formats Importer", "open=[" + file.getAbsolutePath() + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
                ImagePlus imp = WindowManager.getCurrentImage();
                try {
                    runAll(imp, outputImp, file);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }

    private void runAll(ImagePlus imp, ImagePlus outputImp, File file) throws IOException {
        clusters = new ArrayList<>();
        pupae = new ArrayList<>();
        greenAreas = new ArrayList<>();
        stats = new ArrayList<>();
        runAnalysis(imp, file, outputImp);
        IJ.log("Analysis Done. Starting getStats()");
        getStats(imp);
        IJ.log("Stats done. Creating output Imp");
        outputImp = drawOverlay(imp, outputImp);
        IJ.log("Output Imp Created. Creating Results file");
        createResultsFile(String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Results.csv")));
        IJ.log(String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Output.tif")));
        IJ.save(outputImp, String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Output.tif")));
        IJ.run("Close All", " ");
    }

    private void runAnalysis(ImagePlus imp, File file, ImagePlus outputImp) {

        IJ.log(imp.getTitle());
        filename = imp.getTitle();
        filename = filename.replace(".tif","");
        filename = filename.replace(".lif","");
        lifFileName = file.getName();

        ImagePlus[] split = ChannelSplitter.split(imp);
        if (!file.toString().contains(".tif")) {
            getChannelOrder(imp);
        } else if (file.toString().contains(".tif")) {
            getChannelOrderTif(imp);
        }

        Roi[] pupa_all = openPupae(file);
        imp.show();

        Roi outline = getOutline(split[green]);
        ImagePlus clusterMasks = split[red];
        //ImagePlus clusterMasks = WindowManager.getCurrentImage();
        for(int i =0; i< imp.getNFrames();i++) {
            //segmentPupae(split,i);
            pupae.add(pupa_all);
            ArrayList<Roi[]> pupae_clusters = new ArrayList<>();
            ArrayList<Roi[]> greenAreaList = new ArrayList<>();
            //for each pupae
            for (Roi pupa : pupae.get(i)) {
                //Segment Clusters
                //Roi[] cluster = getClusters(split[red], pupa, i);
                clusterMasks.show();
                Roi[] cluster = getClusters(clusterMasks, pupa, i);
                pupae_clusters.add(cluster);
                //for each cluster measure ares of nuclei
                Roi[][] greenArea = getMarkers(cluster,split[green], i,outline);
                greenAreaList.add(mergeArrays(greenArea));
            }
            clusters.add(pupae_clusters);
            greenAreas.add(greenAreaList);
            IJ.log("Frame: "+ i);
        }
    }

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

    private ImagePlus drawOverlay(ImagePlus imp, ImagePlus outputImp){

        //Draw overlay
        for(int i=0; i< imp.getNFrames();i++) {
            ImagePlus impFrame = new Duplicator().run(imp, 1, 2, 1, 1, i+1, i+1);
            ImagePlus overlay = createOutputImage(impFrame, pupae.get(i), clusters.get(i), greenAreas.get(i));
            IJ.log("Overlay: "+ i);
            impFrame.close();
            if (outputImp == null) {
                outputImp = overlay;
            } else {
                outputImp = Concatenator.run(outputImp, overlay);
            }
        }
        return outputImp;
    }

    private Roi[] drawPupae(ImagePlus imp){
        WaitForUserDialog drawPupae = new WaitForUserDialog("Draw ROIs around pupae", "Draw ROIs around pupae and press 'T' to add to ROI manager");
        ImagePlus impZ  = ZProjector.run(imp,"Max");
        impZ.show();
        drawPupae.show();
        Roi[] pupaeRois = roiManager.getRoisAsArray();
        roiManager.reset();
        leftToRight(pupaeRois);
        roiManager.reset();
        impZ.close();
        return pupaeRois;
    }

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

    private void getChannelOrderTif(ImagePlus imp) {
        String info = imp.getInfoProperty();
        for (int i = 0; i < nChannels; i++) {
            if(info.contains("Image #"+i+"|ChannelDescription #"+i+"|LUTName = Green")){green=i;}
            if(info.contains("Image #"+i+"|ChannelDescription #"+i+"|LUTName = Red")){red=i;}
        }
    }

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

    private Roi[][] getMarkers(Roi[] clusters, ImagePlus imp, int frame, Roi outline) {

        imp.show();
        imp.setT(frame+1);
        //IJ.run(imp, "Subtract Background...", "rolling=50 slice");
        Roi[][] output = new Roi[clusters.length][];
        imp.setRoi(outline);
        IJ.setAutoThreshold(imp, "RenyiEntropy dark stack");
        for (int i = 0; i < clusters.length; i++) {
            imp.setRoi(clusters[i]);
            IJ.run(imp, "Analyze Particles...", "size=10-Infinity pixel circularity=0.0-1.00 add");
            output[i] = roiManager.getRoisAsArray();
            roiManager.reset();
        }
        return output;
    }

    private Roi[] combineRois(Roi[] rois, ImagePlus imp) {
        // outer loop to iterate through the array
        Roi temp = null;
        for (int i = 0; i < rois.length; i++) {
            // inner loop to compare and swap elements
            for (int j = i + 1; j < rois.length; j++) {
                if (rois[i].getStatistics().area < rois[j].getStatistics().area) {
                    temp = rois[i];
                    rois[i] = rois[j];
                    rois[j] = temp;
                }
            }
        }
        ArrayList<ArrayList<Roi>> groupedRois = new ArrayList<>();
        for (int k = 0; k < rois.length; k++) {
            if (rois[k] != null) {
                int width = (int) rois[k].getStatistics().roiWidth;
                int[] limits = new int[]{(int) rois[k].getXBase(), (int) rois[k].getXBase() + width};
                ArrayList<Roi> group = new ArrayList<>();
                group.add(rois[k]);
                for (int x = k + 1; x < rois.length; x++) {
                    if (rois[x] != null) {
                        int xcoord = (int) rois[x].getContourCentroid()[0];
                        if (xcoord > limits[0] && xcoord < limits[1]) {
                            group.add(rois[x]);
                            rois[x] = null;
                        }
                    }
                }
                groupedRois.add(group);
            }
        }

        Roi[] output = new Roi[groupedRois.size()];
        roiManager.reset();

        for (ArrayList<Roi> roiArrayList : groupedRois) {

            ImagePlus masks = IJ.createImage("Masks_combiner", "8-bit", imp.getWidth(), imp.getHeight(),
                    imp.getNChannels(), imp.getNSlices(), imp.getNFrames());
            masks.show();

            for (Roi roi : roiArrayList) {
                roiManager.addRoi(roi);
                roiManager.runCommand(masks, "Fill");
                roiManager.reset();
            }

            IJ.run(masks, "Options...", "iterations=2 count=1 black do=Close");
            roiManager.reset();
            IJ.setAutoThreshold(masks, "Huang dark no-reset");
            IJ.run(masks, "Analyze Particles...", "size=50-Infinity pixel include add");
            if (roiManager.getRoisAsArray().length > 0) {
                output[groupedRois.indexOf(roiArrayList)] = roiManager.getRoisAsArray()[0];
            }
            roiManager.reset();
            masks.changes = false;
            masks.close();
        }
        return removeNull(output);
    }

    private Roi[] removeNull(Roi[] array) {
        ArrayList<Roi> list = new ArrayList<>();
        for (Roi roi : array) {
            if (roi != null) {
                if (roi.getStatistics().area > 100000) {
                    list.add(roi);
                }
            }
        }
        Roi[] output = new Roi[list.size()];
        for (int i = 0; i < list.size(); i++) {
            output[i] = list.get(i);
        }
        return output;
    }

    private Roi[] getClusters(ImagePlus imp, Roi pupa, int i) {

        imp.show();
        imp.setT(i+1);
        IJ.run(imp, "Subtract Background...", "rolling=50 slice");
        imp.setRoi(pupa);
        IJ.setAutoThreshold(imp, "Moments dark");
        IJ.run(imp, "Analyze Particles...", "size=" + minArea + "-Infinity pixel circularity=0.0-1.00 add");
        Roi[] allClusters = roiManager.getRoisAsArray();
        if (allClusters.length == 0) {
            return allClusters;
        }
        roiManager.reset();
        //allClusters = trimClusters(pupa, allClusters);
        return allClusters;
    }

    private Roi[] trimClusters(Roi outline, Roi[] clusters){

        for (int i = 0; i< clusters.length;i++){
            FloatPolygon outlineCluster = clusters[i].getInterpolatedPolygon();
            FloatPolygon outlinePoly = outline.getInterpolatedPolygon();
            double feretRatio = clusters[i].getFeretValues()[0]/clusters[i].getFeretValues()[2];
            double dist = 100;
                for (int n = 0; n < outlinePoly.npoints; n++) {
                    for (int j = 0; j < outlineCluster.npoints; j++) {
                        double disttemp = Math.pow(outlinePoly.xpoints[n] - outlineCluster.xpoints[j], 2) +
                                Math.pow(outlinePoly.ypoints[n] - outlineCluster.ypoints[j], 2);
                        disttemp = Math.sqrt(disttemp);
                        if (disttemp < dist) {
                            dist = disttemp;
                        }
                    }
                }
                if ( feretRatio>2) { //dist < 99 &&
                    clusters[i] = null;
                }

        }
        return removeNullClusters(clusters);
    }

    private Roi[] removeNullClusters(Roi[] array) {
        ArrayList<Roi> list = new ArrayList<>();
        for (Roi roi : array) {
            if (roi != null) {
                list.add(roi);
            }
        }
        Roi[] output = new Roi[list.size()];
        for (int i = 0; i < list.size(); i++) {
            output[i] = list.get(i);
        }
        return output;
    }

    private void getStats(ImagePlus imp) {
        ImagePlus[] split = ChannelSplitter.split(imp);
        //Area cluster, area green, intensity green cluster, intensity red cluster, intensity green green, intensity red green
        for (int i = 0; i < imp.getNFrames(); i++) {
            imp.setT(i+1);
            ArrayList<double[]> frameStats = new ArrayList<>();
            for (int j=0; j<pupae.get(i).length;j++){
                double[] pupaeStats = new double[6];
                for (int k = 0; k<clusters.get(i).get(j).length; k++){
                    if(clusters.get(i).get(j).length>0) {
                        pupaeStats[0] = pupaeStats[0] + clusters.get(i).get(j)[k].getStatistics().area;

                        ImageProcessor ipg = split[green].getChannelProcessor();
                        ipg.setRoi(clusters.get(i).get(j)[k]);
                        pupaeStats[1] = pupaeStats[1] + ipg.getStatistics().mean * ipg.getStatistics().area;

                        ImageProcessor ipr = split[red].getChannelProcessor();
                        ipr.setRoi(clusters.get(i).get(j)[k]);
                        pupaeStats[2] = pupaeStats[2] + ipr.getStatistics().mean * ipr.getStatistics().area;
                    }
                }
                for (int m = 0; m<greenAreas.get(i).get(j).length; m++){
                    if(greenAreas.get(i).get(j).length>0) {
                        pupaeStats[3] = pupaeStats[3] + greenAreas.get(i).get(j)[m].getStatistics().area;

                        ImageProcessor ipg = split[green].getChannelProcessor();
                        ipg.setRoi(greenAreas.get(i).get(j)[m]);
                        pupaeStats[4] = pupaeStats[4] + ipg.getStatistics().mean * ipg.getStatistics().area;

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

    private ImagePlus createOutputImage(ImagePlus imp, Roi[] pupae, ArrayList<Roi[]> clusters, ArrayList<Roi[]> greenArea) {

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

        ip.drawString(imp.getTitle(), 50, 50);

        for (int i = 0; i < pupae.length; i++) {
            ip.drawString(letters.get(i), pupae[i].getBounds().x + 150, pupae[i].getBounds().y + 150);
            ip.draw(pupae[i]);
        }

        font = new Font("SansSerif", Font.BOLD, 30);
        ip.setFont(font);
        for (Roi[] roiArray : clusters) {
            for (Roi points : roiArray) {
                ip.draw(points);
                //ip.drawString(String.valueOf(j + 1), roiArray[j].getBounds().x, roiArray[j].getBounds().y);
            }
        }

        ip.setLineWidth(4);
        for (Roi[] roiArray : greenArea) {
            for (Roi points : roiArray) {
                ip.draw(points);
                //ip.drawString(String.valueOf(j + 1), roiArray[j].getBounds().x, roiArray[j].getBounds().y);
            }
        }

        masks.updateAndDraw();
        masks.show();

        ImagePlus[] split = ChannelSplitter.split(imp);
        split[green].show();
        split[red].show();
        IJ.run("Merge Channels...", "c1=[" + split[red].getTitle() + "] c2=[" + split[green].getTitle() +
                    "] c5=[" + masks.getTitle() + "] create");
        return WindowManager.getCurrentImage();
    }

    public void createResultsFile(String name) throws IOException {
        Date date = new Date(); // This object contains the current date value
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy, hh:mm:ss");
        boolean exists = new File(name).exists();
        try {
            FileWriter fileWriter = new FileWriter(name, true);
                BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            if(!exists) {
                bufferedWriter.newLine();
                bufferedWriter.write(formatter.format(date));
                bufferedWriter.newLine();
                bufferedWriter.write("Min Cluster Size: " + minArea);
                bufferedWriter.newLine();
               // bufferedWriter.write("Cellpose model: "+ modelpath.getName());
               // bufferedWriter.newLine();
                //bufferedWriter.write("Labkit model: "+ classifier.getName());
               // bufferedWriter.newLine();

                String heading = "Image,Timepoint, Pupal No,Cluster Area,Marker Area,Mean Cluster Intensity (red),Mean Cluster Intensity (green)," +
                        "Mean Marker Intensity (red),Mean Marker Intensity (green)";
                bufferedWriter.write(String.valueOf(heading));//write header 1
                bufferedWriter.newLine();
            }

            for (int i = 0; i < stats.size(); i++) {
                for (int j = 0; j < stats.get(i).size(); j++) {
                    StringBuilder data = new StringBuilder();
                    data.append(filename).append(",").append(i);
                    data.append(",").append(j);
                    double[] results = stats.get(i).get(j);
                    data.append(",").append(results[0]);
                    data.append(",").append(results[3]);
                    data.append(",").append(results[1]/results[0]);
                    data.append(",").append(results[2]/results[0]);
                    data.append(",").append(results[4]/results[3]);
                    data.append(",").append(results[5]/results[3]);
                    bufferedWriter.write(String.valueOf(data));
                    bufferedWriter.newLine();
                }
            }
            bufferedWriter.close();
            IJ.log("Results file updated");
        } catch (IOException ex) {
            System.out.println(
                    "Error writing to file '" + name + "'");
        }
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
