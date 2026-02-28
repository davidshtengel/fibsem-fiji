package com.shtengel.fib_sem.plugins;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.ThresholdAnalyzer;
import com.shtengel.fib_sem.data.ThresholdData;
import com.shtengel.fib_sem.util.FigBuilder;
import com.shtengel.fib_sem.util.ImageResolver;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.measure.ResultsTable;
import ij.process.ImageProcessor;

/**
 * Fiji/ImageJ command plugin for estimating lower and upper intensity
 * thresholds using histogram-derived quantiles.
 *
 * <p>
 * Thresholds are computed from the cumulative distribution function (CDF)
 * of pixel intensities. The user specifies fractions of low- and high-intensity
 * pixels to discard.
 * </p>
 *
 * <p>
 * <b>ROI Support:</b> When an ROI is active, only pixels within the ROI are
 * analyzed.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Calculate Min/Max Thresholds")
public class Thresholds implements Command {
	
    @Parameter(type = ItemIO.INPUT, required = false)
    private ImagePlus imp;

    @Parameter(label = "Lower CDF Threshold (fraction)", min = "0.0", max = "1.0")
    private double thrMin = 0.001;
    
    @Parameter(label = "Upper CDF Threshold (fraction)", min = "0.0", max = "1.0")
    private double thrMax = 0.001;

    @Parameter(label = "Number of bins", min = "2", max = "65536")
    private int nbins = 256;
        
    @Parameter(label = "Logarithmic histogram")
    private boolean logHist = false;
    
    @Parameter(label = "Export as figure")
    private boolean saveFig = false;
    
    /**
     * Executes threshold computation and optional visualization.
     *
     * <p>
     * Processing steps:
     * <ol>
     *   <li>Validate image and processor</li>
     *   <li>Apply ROI if present</li>
     *   <li>Compute thresholds via {@link ThresholdAnalyzer}</li>
     *   <li>Optionally display PDF/CDF plot and export table</li>
     * </ol>
     * </p>
     */
    @Override
    public void run() {
        // Validate image presence
    	imp = ImageResolver.resolveSourceImage(imp);
        if (imp == null) {  
            IJ.error("Threshold Analysis", "No image open.");
            return;
        }
        
        // Get the ImageProcessor
        ImageProcessor ip = imp.getProcessor();
        if (ip == null) {
            IJ.error("Threshold Analysis", "Could not access image processor.");
            return;
        }
        
        // Get the ROI, if one exists
        Roi roi = imp.getRoi();
        
        // If there's an ROI, set it on the processor, so statistics calculations respect ROI bounds
        if (roi != null) {
        	ip.setRoi(roi);
        	String roiName = roi.getName() != null ? roi.getName() : "unnamed ROI";
            IJ.log("Processing ROI: " + roiName + " (type: " + roi.getTypeAsString() + ")");
        } else {
        	IJ.log("Processing entire image (no ROI selected)");
        }
        
        // Compute thresholds and get fields for plots, tables, and logs
        ThresholdData result = ThresholdAnalyzer.computeThresholds(ip, thrMin, thrMax, nbins);
        double minInt = result.getMinIntensity();
        double maxInt = result.getMaxIntensity();
        double minThr = result.getMinThreshold();
        double maxThr = result.getMaxThreshold();
        double[] pdf = result.getPDF();
        double[] cdf = result.getCDF();
        
        // Verify computation was successful
        if (Double.isNaN(minThr) || Double.isNaN(maxThr)) {
        	IJ.error("Threshold Analysis", "Failed to compute valid thresholds. Check image data.");
            return;
        }

        // Display results
        Plot threshold = displayThresholdPlot(minInt, maxInt, minThr, maxThr, pdf, cdf);
        
        if (saveFig) {
        	// Get directory from the source image
    	    FileInfo fi = imp.getOriginalFileInfo();
    	    String dir = "";
    	    
    	    if (fi != null && fi.directory != null) {
    	        dir = fi.directory;
    	    }
    	    
    	    String baseName = ImageResolver.getBaseName(imp);
    		
    		FigBuilder.createAndSave(threshold, 
    				"Threshold Plot",
    				dir + baseName + "_threshold_plot.png");
        }
        
        threshold.show();
        
        // Log threshold values
        IJ.log("Min: " + minThr);
        IJ.log("Max: " + maxThr);
        IJ.showStatus("Thresholding complete.");
    }
    
    /**
     * Displays a combined PDF/CDF plot with threshold markers and
     * exports the underlying histogram data to a {@link ResultsTable}.
     */
    private Plot displayThresholdPlot(double minInt, 
    		double maxInt, 
    		double minThr, 
    		double maxThr, 
    		double[] pdf, 
    		double[] cdf) {
    	
    	double binWidth = (maxInt - minInt) / nbins;
    	
    	// Calculate bin centers
    	double[] binCenters = new double[nbins]; 
        for (int i = 0; i < nbins; i++) { 
            binCenters[i] = minInt + (i + 0.5) * binWidth; 
        }
        
        // Normalize PDF for display only (scale to [0, 1] range)
        double[] pdfNormalized = normalizePDF(pdf);
        
        // Apply log transform if requested (to normalized version for display)
        double[] pdfToDisplay = logHist ? logArray(pdfNormalized) : pdfNormalized;
        
        String title = "Threshold Analysis: " + (logHist ? "Log PDF" : "PDF") + " & CDF";
        String roiInfo = (imp.getRoi() != null) ? "[ROI]" : "";
        Plot plot = new Plot(title + roiInfo, 
        				"Intensity", 
        				"Probability"
        ); 
                
        // Add PDF (as bars)
        plot.setColor("blue"); 
        plot.add("bar", binCenters, pdfToDisplay);
        
        // Add CDF (as line)
        plot.setColor("magenta");
        plot.add("line", binCenters, cdf); 
        
        // Add threshold lines
        plot.setColor("red");
        plot.drawDottedLine(minThr, 0.0, minThr, 1.0, 1);
        plot.setColor("cyan");
        plot.drawDottedLine(maxThr, 0.0, maxThr, 1.0, 1);
        
        // Add legend
        plot.setColor("black");
        plot.addLegend(String.format(
        		"PDF (normalized)\nCDF\nLower threshold = %.3g\nUpper threshold = %.3g",
        	    minThr, maxThr
        ));
        
        // Create a ResultsTable with original values for export
        ResultsTable rt = new ResultsTable();
        for (int i = 0; i < nbins; i++) {
            rt.incrementCounter();
            rt.addValue("Intensity", binCenters[i]);
            rt.addValue("PDF (original)", pdf[i]);
            rt.addValue("PDF (normalized)", pdfNormalized[i]);
            rt.addValue("CDF", cdf[i]);
        }
        
        rt.show("Threshold Analysis Data: " + imp.getTitle());
        
        return plot;
    }

    /**
     * Scales a PDF to the [0, 1] range for display purposes.
     */
    private static double[] normalizePDF(double[] pdf) {
        double[] normalized = new double[pdf.length];
        
        // Find max value
        double maxVal = 0;
        for (double v : pdf) {
            if (v > maxVal) maxVal = v;
        }
        
        // Scale to 0-1 range
        if (maxVal > 0) {
            for (int i = 0; i < pdf.length; i++) {
                normalized[i] = pdf[i] / maxVal;
            }
        }
        
        return normalized;
    }
    
    /**
     * Applies a base-10 logarithmic transform to an array.
     * Zero or negative values are mapped to zero.
     * @param arr	Input array
     * @return out 	Log-transformed array
     */
    private static double[] logArray(double[] arr) {
        double[] out = new double[arr.length];
        for (int i = 0; i < arr.length; i++) {
            out[i] = arr[i] > 0 ? Math.log10(arr[i]) : 0;
        }
        return out;
    }

}
