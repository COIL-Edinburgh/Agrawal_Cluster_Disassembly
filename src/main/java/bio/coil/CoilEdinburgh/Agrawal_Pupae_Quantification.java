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
import ij.plugin.ChannelSplitter;
import ij.plugin.Concatenator;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.LUT;
import io.scif.services.DatasetIOService;
import io.scif.services.FormatService;
import loci.formats.ChannelSeparator;
import loci.formats.FormatException;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imagej.roi.ROIService;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.awt.color.ColorSpace;
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
@Plugin(type = Command.class, menuPath = "Plugins>Users Plugins>Agrawal Pupae Quantification")
public class Agrawal_Pupae_Quantification<T extends RealType<T>> implements Command {
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

    @Parameter(label = "Edit? ")
    public boolean edit;

    RoiManager roiManager;
    double pixelSize;
    int nChannels;
    String lifFileName;
    String filename;
    ArrayList<Roi[]> clusters;
    ArrayList<Roi[][]> nuclei;
    Roi[] pupae;
    ArrayList<String> letters;
    int green = 0;
    int tl = 1;
    int red = 2;

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
            if (file.toString().endsWith(".lif") && !file.toString().contains(".tif")) {
                for (int i = 1; i <= getNSeries(file); i++) {
                    IJ.run("Bio-Formats Importer", "open=[" + file.getAbsolutePath() + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_" + i);
                    ImagePlus imp = WindowManager.getCurrentImage();
                    runAll(imp, outputImp, file);
                }

            } else if (file.toString().contains(".tif") && !file.toString().contains("Output")) {
                IJ.run("Bio-Formats Importer", "open=[" + file.getAbsolutePath() + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
                ImagePlus imp = WindowManager.getCurrentImage();
                runAll(imp, outputImp, file);
            }
        }
    }

    private void runAll(ImagePlus imp, ImagePlus outputImp, File file){
        runAnalysis(imp, file, outputImp);
        if(edit) {
            EditSegmentation ES = new EditSegmentation(imp, roiManager, pupae, clusters);
            ES.run();
            pupae = ES.pupae;
            clusters = ES.clusters;
        }
        drawOverlay(imp, outputImp);
        makeResults();
        IJ.log(String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Output.tif")));
        IJ.save(outputImp, String.valueOf(Paths.get(String.valueOf(filePath), filename + "_Output.tif")));
        IJ.run("Close All", " ");
    }

//    private File saveRois(){
//
//        //create folder
//        String roiFolderName = String.valueOf(Paths.get(String.valueOf(filePath), filename+"_ROIs"));
//        File roiFolder = new File(roiFolderName);
//        if (!roiFolder.exists()) {
//            new File(roiFolderName).mkdir();
//        }
//            roiManager.reset();
//            for (int j=0; j< pupae.length; j++) {
//                roiManager.addRoi(pupae[j]);
//                String pupaeName = String.valueOf(Paths.get(roiFolderName,filename+"_"+j));
//                roiManager.save(pupaeName+".roi");
//                roiManager.reset();
//                for(int k=0; k<clusters.get(j).length; k++){
//                    if(clusters.get(j)[k]!=null) {
//                        roiManager.addRoi(clusters.get(j)[k]);
//                    }
//                }
//                if(roiManager.getRoisAsArray().length==1) {
//                    roiManager.save(pupaeName + "_clusters.roi");
//                }else if (roiManager.getRoisAsArray().length>1){
//                    roiManager.save(pupaeName + "_clusters.zip");
//                }
//                roiManager.reset();
//            }
//            return roiFolder;
//    }

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
        segmentPupae(split);

        clusters = new ArrayList<>();
        nuclei = new ArrayList<>();
        //for each pupae
        for (Roi pupa : pupae) {
            //Segment Clusters
            Roi[] cluster = getClusters(split[green], pupa);
            this.clusters.add(cluster);
            //for each cluster measure ares of nuclei
            if (nChannels == 3) {
                nuclei.add(getNuclei(cluster, split[red]));
            }
        }

    }

    private void drawOverlay(ImagePlus imp, ImagePlus outputImp){

        //Draw overlay
        ImagePlus overlay = createOutputImage(imp, pupae, clusters);
        if (outputImp == null) {
            outputImp = overlay;
        } else {
            outputImp = Concatenator.run(outputImp, overlay);
        }
    }

    private void makeResults(){
        //output results
        try {
            createResultsFile("Results");
            if (nChannels == 3) {
                createResultsFile("Nuclei");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void segmentPupae(ImagePlus[] split) {
        //segment pupae
        nChannels = split.length;

        split[tl].show();

        IJ.setAutoThreshold(split[tl], "Huang no-reset");
        IJ.run(split[tl], "Analyze Particles...", "size=500-Infinity pixel show=Masks include");
        ImagePlus masks = WindowManager.getCurrentImage();
        IJ.run("Options...", "iterations=1 count=1");
        IJ.run(masks, "Watershed", "");
        IJ.run(masks, "Analyze Particles...", "size=0-Infinity pixel add include");
        masks.changes = false;
        masks.close();
        Roi[] pupaeRois = roiManager.getRoisAsArray();
        pupae = combineRois(pupaeRois, split[tl]);
        leftToRight(pupae);
        roiManager.reset();
    }

    private void getChannelOrder(ImagePlus imp) {
        LUT[] luts = imp.getLuts();
        for (int i = 0; i < luts.length; i++) {
            String string = luts[i].toString();
            if (string.contains("green")) {
                green = i;
            } else if (string.contains("white")) {
                tl = i;
            } else if (string.contains("red")) {
                red = i;
            }
        }
    }

    private void getChannelOrderTif(ImagePlus imp) {
        String info = imp.getInfoProperty();
        for (int i = 0; i < nChannels; i++) {
            if(info.contains("Image #"+i+"|ChannelDescription #"+i+"|LUTName = Green")){green=i;}
            if(info.contains("Image #"+i+"|ChannelDescription #"+i+"|LUTName = Gray")){tl=i;}
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

    private Roi[][] getNuclei(Roi[] clusters, ImagePlus imp) {
        imp.show();
        IJ.run(imp, "Subtract Background...", "rolling=25 slice");
        Roi[][] output = new Roi[clusters.length][];
        for (int i = 0; i < clusters.length; i++) {
            imp.setRoi(clusters[i]);
            IJ.setAutoThreshold(imp, "Yen dark");
            IJ.run(imp, "Analyze Particles...", "size=10-500 pixel circularity=0.05-1.00 exclude add");
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

    private int getNSeries(File file) {
        ImageProcessorReader r = new ImageProcessorReader(
                new ChannelSeparator(LociPrefs.makeImageReader()));
        try {
            r.setId(file.getAbsolutePath());
        } catch (FormatException | IOException e) {
            throw new RuntimeException(e);
        }
        return r.getSeriesCount();
    }

    private Roi[] getClusters(ImagePlus imp, Roi pupa) {
        imp.show();
        imp.setRoi(pupa);
        IJ.setAutoThreshold(imp, "Triangle dark");
        IJ.run(imp, "Analyze Particles...", "size=" + minArea + "-Infinity pixel circularity=0.05-1.00 exclude add");
        IJ.run(imp, "Enhance Contrast", "saturated=0.35");
        Roi[] allClusters = roiManager.getRoisAsArray();
        if (allClusters.length == 0) {
            return allClusters;
        }
        roiManager.reset();
        return allClusters;
    }

    private double[][] getNucleiAreas(int i) {
        double[][] averages = new double[clusters.get(i).length][2];
        for (int j = 0; j < clusters.get(i).length; j++) {
            int count = 0;
            for (int k = 0; k < nuclei.get(i)[j].length; k++) {
                if (nuclei.get(i)[j][k] != null) {
                    averages[j][0] = averages[j][0] + nuclei.get(i)[j][k].getStatistics().area;
                    count++;
                }
            }
            averages[j][1] = count;
        }
        return averages;
    }

    private ImagePlus createOutputImage(ImagePlus imp, Roi[] pupae, ArrayList<Roi[]> clusters) {

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
            ip.drawString(letters.get(i), pupae[i].getBounds().x + 50, pupae[i].getBounds().y + 50);
            ip.draw(pupae[i]);
        }

        font = new Font("SansSerif", Font.BOLD, 30);
        ip.setFont(font);
        for (Roi[] roiArray : clusters) {
            for (int j = 0; j < roiArray.length; j++) {
                ip.draw(roiArray[j]);
                ip.drawString(String.valueOf(j + 1), roiArray[j].getBounds().x, roiArray[j].getBounds().y);
            }
        }

        for (Roi[][] roiArray : nuclei) {
            for (Roi[] rois : roiArray) {
                for (Roi points : rois) {
                    ip.draw(points);
                }
            }
        }


        masks.updateAndDraw();
        masks.show();

        ImagePlus[] split = ChannelSplitter.split(imp);
        split[green].show();
        split[tl].show();
        if (nChannels == 3) {
            split[red].show();
            IJ.run("Merge Channels...", "c1=[" + split[red].getTitle() + "] c2=[" + split[green].getTitle() + "] c4=[" + split[tl].getTitle() +
                    "] c5=[" + masks.getTitle() + "] create");
        } else {
            //c1=["+split[red.getTitle()+"]
            IJ.run("Merge Channels...", " c2=[" + split[green].getTitle() + "] c4=[" + split[tl].getTitle() +
                    "] c5=[" + masks.getTitle() + "] create");
        }
        return WindowManager.getCurrentImage();
    }

    public void createResultsFile(String name) throws IOException {
        Date date = new Date(); // This object contains the current date value
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy, hh:mm:ss");
        String foldername = filePath.getName();
        String CreateName = String.valueOf(Paths.get(String.valueOf(filePath), foldername + "_" + name + ".csv"));
        Boolean exists = new File(CreateName).exists();
        try {
            FileWriter fileWriter = new FileWriter(CreateName, true);
                BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            if(!exists) {
                bufferedWriter.newLine();
                bufferedWriter.write(formatter.format(date));
                bufferedWriter.newLine();
                bufferedWriter.write("Min Cluster Size: " + minArea);
                bufferedWriter.newLine();
                bufferedWriter.newLine();

                StringBuilder heading = new StringBuilder();
                if (Objects.equals(name, "Results")) {
                    heading.append("Image,Pupal No,Cluster 1,Cluster 2,Cluster 3,Cluster 4,Cluster 5,Cluster 6,Cluster 7,Cluster 8,Cluster 9,Cluster 10,Average Cluster," +
                            "Pupal Area,Ratio(Cluster/Pupae),");
                }
                if (Objects.equals(name, "Nuclei")) {
                    heading.append("Image,Pupal No,Cluster No,Total N,Total Area,Area 1,Area 2,Area 3,Area 4,Area 5,Area 6,Area 7,Area 9,Area 10,Area 11,Area 12");
                }
                bufferedWriter.write(String.valueOf(heading));//write header 1
                bufferedWriter.newLine();
            }
            int count = 0;
            double average = 0;

            if (Objects.equals(name, "Results")) {
                for (int i = 0; i < pupae.length; i++) {
                    StringBuilder data = new StringBuilder();
                    data.append(filename).append(",").append(letters.get(i));
                    for (int j = 0; j < 10; j++) {
                        if (j < clusters.get(i).length) {
                            data.append(",").append(clusters.get(i)[j].getStatistics().area);
                            average = average + clusters.get(i)[j].getStatistics().area;
                            count = count + 1;
                        } else {
                            data.append(",");
                        }
                    }

                    double ratio = (average / count) / pupae[i].getStatistics().area;
                    data.append(",").append(average / count).append(",").append(pupae[i].getStatistics().area).append(",").append(ratio).append(",");
                    bufferedWriter.write(String.valueOf(data));
                    bufferedWriter.newLine();
                }
            }

            if (Objects.equals(name, "Nuclei")) {
                for (int i = 0; i < pupae.length; i++) {
                    for (int j = 0; j < clusters.get(i).length; j++) {
                        StringBuilder data = new StringBuilder();
                        data.append(filename).append(",").append(letters.get(i));
                        data.append(",Cluster ").append(j);
                        double[][] nucleiStats = getNucleiAreas(i);
                        data.append(",").append(nucleiStats[j][1]).append(",").append(nucleiStats[j][0]);
                        for (int k = 0; k < nuclei.get(i)[j].length; k++) {
                            double area = nuclei.get(i)[j][k].getStatistics().area;
                            data.append(",").append(area);
                        }
                        bufferedWriter.write(String.valueOf(data));
                        bufferedWriter.newLine();
                    }
                }
            }
            bufferedWriter.close();
            IJ.log("Results file updated");
        } catch (IOException ex) {
            System.out.println(
                    "Error writing to file '" + CreateName + "'");
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
        ij.command().run(Agrawal_Pupae_Quantification.class, true);
    }

}
