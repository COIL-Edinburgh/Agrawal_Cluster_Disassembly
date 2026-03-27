package bio.coil.CoilEdinburgh;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.Concatenator;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

public class EditSegmentation {

    ImagePlus imp;
    RoiManager roiManager;
    Roi[] pupae;
    ArrayList<Roi[]> clusters;

    EditSegmentation(ImagePlus imp, RoiManager roiManager, Roi[] pupae, ArrayList<Roi[]> clusters) {
        this.imp = imp;
        this.roiManager = roiManager;
        this.clusters = clusters;
        this.pupae = pupae;
    }

    void run() {

        ImagePlus slideImp = imp;
        Prefs.showAllSliceOnly = false;

        slideImp.show();
        for (int j = 0; j < pupae.length; j++) {
            roiManager.addRoi(pupae[j]);
        }
        roiManager.runCommand(slideImp, "Show All");
        WaitForUserDialog next = new WaitForUserDialog("Edit Clusters?");
        next.show();
        Roi[] newPupae = roiManager.getRoisAsArray();
        roiManager.reset();
        leftToRight(newPupae);
        pupae = newPupae;

        for (Roi[] allCluster : clusters) {
            for (int j = 0; j < allCluster.length; j++) {
                roiManager.addRoi(allCluster[j]);
            }
        }
        roiManager.runCommand(slideImp, "Show All");
        WaitForUserDialog next2 = new WaitForUserDialog("Next Image?");
        next2.show();

        Roi[] newClusters = roiManager.getRoisAsArray();
        clusters = arrangeClusters(newClusters, newPupae);

        roiManager.reset();
    }

    ArrayList<Roi[]> arrangeClusters(Roi[] clusters, Roi[] pupae) {
        ArrayList<Roi[]> output = new ArrayList<>();

        for (Roi pupa : pupae) {
            ArrayList<Roi> tempArrayList = new ArrayList<>();
            for (Roi cluster : clusters) {
                if (pupa.containsPoint(cluster.getContourCentroid()[0], cluster.getContourCentroid()[1])) {
                    tempArrayList.add(cluster);
                }
            }
            Roi[] tempArray = toRoiArray(tempArrayList);
            output.add(tempArray);
        }

        return output;
    }

    Roi[] toRoiArray(ArrayList<Roi> input) {
        Roi[] output = new Roi[input.size()];
        for (int i = 0; i < input.size(); i++) {
            output[i] = input.get(i);
        }
        return output;
    }

    void leftToRight(Roi[] rois) {
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

}
