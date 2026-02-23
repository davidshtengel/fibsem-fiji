package com.shtengel.fib_sem.plugins;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.GradientMapAnalyzer;
import com.shtengel.fib_sem.data.GradientMapData;
import com.shtengel.fib_sem.util.ImageResolver;
import com.shtengel.fib_sem.util.FigBuilder;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.process.ImageProcessor;

/**
 * Fiji/ImageJ command plugin for computing gradient magnitude maps
 * from FIB-SEM images.
 *
 * <p>
 * Gradients are computed using central finite differences with optional
 * pre-smoothing and optional normalization by local intensity.
 * </p>
 *
 * <p>
 * <b>ROI Support:</b> When an ROI is active, only pixels within the ROI are
 * analyzed.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Map Gradient")
public class GradientMap implements Command {

	@Parameter(type = ItemIO.INPUT, required = false)
    private ImagePlus imp;

    @Parameter(label = "Perform smoothing")
    private boolean performSmoothing = true;

    @Parameter(label = "Normalize gradient")
    private boolean normalize = false;

    @Parameter(label = "Lower display threshold", min = "0.0", max = "1.0")
    private double thrMinDisp = 1e-3;

    @Parameter(label = "Upper display threshold", min = "0.0", max = "1.0")
    private double thrMaxDisp = 1e-3;

    @Parameter(label = "Export as figure")
    private boolean saveFig = false;
    
    /**
     * Executes the gradient map computation workflow.
     * 
     * <p>
     * Processing steps:
     * <ol>
     *   <li>Validates image and processor availability</li>
     *   <li>Identifies and logs active ROI if present</li>
     *   <li>Computes gradient map using {@link GradientMapAnalyzer}</li>
     *   <li>Displays result if requested</li>
     * </ol>
     * </p>
     * 
     * <p>
     * When an ROI is active, only pixels within the ROI are processed.
     * The output gradient map will have dimensions matching the ROI bounds.
     * </p>
     * 
     * @see GradientMapAnalyzer#computeGradientMap
     */
    @Override
    public void run() {
        // Validate image presence
    	imp = ImageResolver.resolveSourceImage(imp);
        if (imp == null) {  
            IJ.error("Gradient Mapping", "No image open.");
            return;
        }
        
        // Get the ImageProcessor
        ImageProcessor ip = imp.getProcessor();
        if (ip == null) {
            IJ.error("Gradient Mapping", "Could not access image processor.");
            return;
        }
        
        // Get the ROI, if one exists
        Roi roi = imp.getRoi();
        
        // Set ROI on the processor, if one exists
        if (roi != null) {
        	ip.setRoi(roi);
        	String roiName = roi.getName() != null ? roi.getName() : "unnamed ROI";
            IJ.log("Processing ROI: " + roiName + " (type: " + roi.getTypeAsString() + ")");
        } else {
        	IJ.log("Processing entire image (no ROI selected)");
        }

        // Compute gradient map
        GradientMapData result = GradientMapAnalyzer.computeGradientMap(ip, performSmoothing, normalize);
        ImageProcessor gradient = result.getGradient();

        // Display result
        if (gradient != null) {
            ImagePlus gradImp = new ImagePlus("Gradient Map: " + imp.getTitle(), gradient);
            
            if (saveFig) {
            	// Get directory from the source image
        	    FileInfo fi = imp.getOriginalFileInfo();
        	    String dir = "";
        	    
        	    if (fi != null && fi.directory != null) {
        	        dir = fi.directory;
        	    }
        	    
        	    String baseName = getBaseName(imp);
            	FigBuilder.createAndSave(gradImp,
            					"Gradient Map",
            					dir + baseName + "_gradient_map.png");
            }
            
            gradImp.show();
        }
        
        IJ.log("Gradient map computed with parameters.");
    }
    
    /**
	 * Helper to get base filename without extension
	 */
	private String getBaseName(ImagePlus imp) {
	    String title = imp.getTitle();
	    int dot = title.lastIndexOf('.');
	    return (dot > 0) ? title.substring(0, dot) : title;
	}
}
