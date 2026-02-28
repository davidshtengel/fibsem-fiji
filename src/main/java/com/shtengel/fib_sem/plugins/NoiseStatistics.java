package com.shtengel.fib_sem.plugins;

import java.awt.Color;
import java.util.Arrays;
import java.util.LinkedHashMap;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.GradientMapAnalyzer;
import com.shtengel.fib_sem.core.NoiseStatisticsAnalyzer;
import com.shtengel.fib_sem.data.NoiseStatisticsData;
import com.shtengel.fib_sem.util.FigBuilder;
import com.shtengel.fib_sem.util.ImageResolver;
import com.shtengel.fib_sem.util.ParamPersister;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * ImageJ plugin for characterizing noise properties in FIB-SEM images.
 * 
 * <p>
 * Fiji/ImageJ command plugin for analyzing the relationship between FIB-SEM image 
 * intensity and noise variance. The analysis is based on the
 * variance-mean relationship, where variance is expected to be linearly related 
 * to mean intensity.
 * </p>
 * 
 * <p>
 * The plugin performs the following analysis:
 * <ul>
 *   <li>Smooths the image and computes local intensity means and variances</li>
 *   <li>Bins the data by intensity and plots variance vs. mean</li>
 *   <li>Performs linear regression with two models:
 *     <ul>
 *       <li><b>Free fit (SNR):</b> var = slope × mean + intercept</li>
 *       <li><b>Dark count fit (SNR1):</b> var = slope × (mean - darkCount)</li>
 *     </ul>
 *   </li>
 *   <li>Computes signal-to-noise ratios for both models</li>
 * </ul>
 * </p>
 * 
 * <p>
 * <b>ROI Support:</b> When an ROI is active, only pixels within the ROI are
 * analyzed.
 * </p>
 * 
 * @see NoiseStatisticsAnalyzer
 */
@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Noise Statistics Analysis (Single Image)")
public class NoiseStatistics implements Command {
	
	private double thrMinAnalysis;
    private double thrMaxAnalysis;
    private int nbinsAnalysis;
    private double gradientThreshold;
    private boolean displaySNR;
    private double darkCount;
    private boolean displaySNR1;
	// Fixed display-range parameters (not exposed in the dialog); used only to compute
	// the intensity axis limits for the variance-mean plot.
    private double thrMinDisp = 0.001;
    private double thrMaxDisp = 0.001;
    private int nbinsDisp = 256;
    private boolean showMaskVisualization;
    private boolean saveFigs;
    private ImagePlus imp;
    
    @Override
    public void run() {
    	// Get active image
        imp = ImageResolver.resolveSourceImage();
        if (imp == null) {
            IJ.error("Noise Statistics Analysis", "No image open.");
            return;
        }
        
        if (!showDialog()) {
            return;
        }
        
        ImageProcessor ip = imp.getProcessor();
        if (ip == null) {
            IJ.error("Noise Statistics Analysis", "Could not access image processor.");
            return;
        }

		Roi roi = imp.getRoi();
        
        // If there's an ROI, set it on the processor - statistics calculations will respect the ROI bounds
        if (roi != null) {
        	ip.setRoi(roi);
        	String roiName = roi.getName() != null ? roi.getName() : "unnamed ROI";
            IJ.log("Processing ROI: " + roiName + " (type: " + roi.getTypeAsString() + ")");
        } else {
        	IJ.log("Processing entire image (no ROI selected)");
        }
		
		double[] thresholdsDisp = new double[] {thrMinDisp, thrMaxDisp};
		double[] thresholdsAnalysis = new double[] {thrMinAnalysis, thrMaxAnalysis};
		
		// Compute noise statistics
		NoiseStatisticsData result = NoiseStatisticsAnalyzer.computeNoiseStatistics(
			ip,
			darkCount,
			null, // default smoothing kernel
			nbinsDisp,
			thresholdsDisp,
			nbinsAnalysis,
			thresholdsAnalysis,
			gradientThreshold,
			null // no filter mask
		);
        
        logResults(imp, result);

		// Persist computed I0 for contrast calculation
		double i0 = result.getI0();
		ParamPersister.set(imp, "C_ranSNR", true);
		ParamPersister.set(imp, "C_i0", i0);
        
        Plot noiseDist = plotNoiseDistribution(result, darkCount, displaySNR, displaySNR1, imp.getTitle());
        
        // Create mask visualization if requested
        ImagePlus mask = null;
        if (showMaskVisualization) {
            mask = drawMaskVisual(imp, result);            
        }
        
        if (saveFigs) {
    	    FileInfo fi = imp.getOriginalFileInfo();
    	    String dir = (fi != null && fi.directory != null) ? fi.directory : "";
    	    String baseName = ImageResolver.getBaseName(imp);
    	    
        	FigBuilder.createAndSave(noiseDist, "Noise Distribution", dir + baseName + "_noise_distribution.png");
        	
        	if (mask != null) {        		
        		FigBuilder.createAndSave(mask, "Pixel Exclusion Mask", dir + baseName + "_exclusion_mask.png");
        	}
        }
        
        noiseDist.show();
        if (showMaskVisualization) {        	
        	mask.show();
        }
        
        IJ.showStatus("Noise analysis complete.");
    }
    
    private boolean showDialog() {
		getAllPersistedParams();
		
    	GenericDialog gd = new GenericDialog("Noise Statistics Analysis (Single Image)");
    	
    	gd.addMessage("These parameters are set to exclude the pixels with highest and lowest intensities\nfrom the analysis. Set the thresholds to 0.0 to include all pixels.");
        gd.addNumericField("CDF threshold (lower)", thrMinAnalysis, 3, 6, "");
        gd.addNumericField("CDF threshold (upper)", thrMaxAnalysis, 3, 6, "");
        gd.addNumericField("Number of bins", nbinsAnalysis, 0, 6, "");
        
        gd.addMessage("This parameter is set to exclude the pixels with highest local gradient from the\nanalysis. Set to 1.0 to include all pixels. Set to 0.50 to only consider 50% of pixels\nwith the lowest local intensity gradient.");
        gd.addNumericField("Gradient threshold (upper)", gradientThreshold, 2, 6, "");
        gd.addCheckbox("Show SNR0 (free fit)", displaySNR);
        gd.addNumericField("User-defined Dark Count (intensity at zero variance)", darkCount, 1, 6, "");
        gd.addCheckbox("Show SNR1 (with user-defined dark count)", displaySNR1);
        
        gd.addMessage("Visualization options:");
        gd.addCheckbox("Show mask visualization", showMaskVisualization);
        
	    gd.addMessage("Export:");
	    gd.addCheckbox("Save image/plot as titled figures", saveFigs);
        
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        
        // Retrieve values
        thrMinAnalysis = gd.getNextNumber();
        thrMaxAnalysis = gd.getNextNumber();
        nbinsAnalysis = (int) gd.getNextNumber();
        gradientThreshold = gd.getNextNumber();
        displaySNR = gd.getNextBoolean();
        darkCount = gd.getNextNumber();
        displaySNR1 = gd.getNextBoolean();
        showMaskVisualization = gd.getNextBoolean();
        saveFigs = gd.getNextBoolean();
        
        // Validate inputs
        if (thrMinAnalysis < 0 || thrMinAnalysis > 1 || thrMaxAnalysis < 0 || thrMaxAnalysis > 1) {
            IJ.error("Invalid threshold values. Must be between 0 and 1.");
            return false;
        }
        
        if (nbinsAnalysis < 10 || nbinsAnalysis > 1024) {
            IJ.error("Number of bins must be between 10 and 1024.");
            return false;
        }

        if (gradientThreshold < 0 || gradientThreshold > 1) {
			IJ.error("Invalid gradient threshold. Must be between 0 and 1.");
            return false;
		}

		setAllPersistedParams();        
        return true;
    }
    
    /**
     * Creates a composite visualization showing which pixels were excluded
     * from analysis due to intensity thresholds and gradient filtering.
     */
    private ImagePlus drawMaskVisual(ImagePlus imp, NoiseStatisticsData result) {
    	String imageTitle = imp.getTitle();
    	ImageProcessor ip = imp.getProcessor();
		ImageProcessor ipToProcess = ImageResolver.cropToRoiIfPresent(ip);
		int width = ipToProcess.getWidth();
        int height = ipToProcess.getHeight();
        
        // Get analysis parameters from result
        double[] rangeAnalysis = result.getRangeAnalysis();
        double minThreshold = rangeAnalysis[0];
        double maxThreshold = rangeAnalysis[1];
        
        // Recompute smoothed image and gradients for visualization
        FloatProcessor fpSmoothed = ipToProcess.convertToFloatProcessor();
        fpSmoothed = (FloatProcessor) fpSmoothed.duplicate();
        GradientMapAnalyzer.applyDefaultSmoothing(fpSmoothed);
        float[] smoothed = (float[]) fpSmoothed.getPixels();
        float[] gradients = NoiseStatisticsAnalyzer.computeGradientMagnitudes(ipToProcess.convertToFloatProcessor());
        
        FloatProcessor fp = ipToProcess.convertToFloatProcessor();
        double min = fp.getMin();
        double max = fp.getMax();
        
        float gradientCutoff = NoiseStatisticsAnalyzer.computeGradientCutoff(gradients, gradientThreshold);
        
        ColorProcessor maskVis = new ColorProcessor(width, height);
        
        // Create color overlay
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int idx = y * width + x;
                float val = smoothed[idx];
                
                boolean isBorder = (y == 0 || y == height - 1 || x == 0 || x == width - 1);
                                
                // Check exclusion criteria
                boolean belowMin = (val < minThreshold);
                boolean aboveMax = (val > maxThreshold);
                boolean highGradient = !isBorder && (gradients[idx] > gradientCutoff);
                
                int color;
                
                if (isBorder) {
                	color = 0x404040;
                } else if (belowMin) {
                	color = 0xFF0000;
                } else if (aboveMax) {
                	color = 0x00FFFF;
                } else if (highGradient) {
                	color = 0x0000FF;
                } else {
                	float pixelValue = fp.getf(x, y);
                    int scaledValue = (int) (255.0 * (pixelValue - min) / (max - min));
                    scaledValue = Math.min(255, Math.max(0, scaledValue));    
                    color = (scaledValue << 16) | (scaledValue << 8) | scaledValue;                
				}
                maskVis.set(x, y, color);
            }
        }
        
        ImagePlus maskImp = new ImagePlus("Mask Visualization: " + imageTitle, maskVis);
        
        // Create and display results table
        displayCutoffValuesTable(minThreshold, maxThreshold, gradientCutoff, imageTitle);
        
        // Add overlay text explaining colors
        IJ.log("=== Mask Visualization Legend ===");
        IJ.log("Cyan: Pixels excluded (above intensity threshold)");
        IJ.log("Red: Pixels excluded (below intensity threshold)");
        IJ.log("Blue: Pixels excluded (high local gradient)");
        IJ.log("Grayscale: Pixels included in analysis");
        IJ.log("Dark Gray: Pixels excluded (border, not analyzed)");
        
        return maskImp;
    }
    
    /** Displays a results table showing the cutoff values for each exclusion category. */
    private void displayCutoffValuesTable(double minThreshold, double maxThreshold, float gradientCutoff, String imageTitle) {
        ij.measure.ResultsTable rt = new ij.measure.ResultsTable();
        
        rt.incrementCounter();
        rt.addLabel("Below Min Intensity (Red)");
        rt.addValue("Cutoff Value", minThreshold);

        rt.incrementCounter();
        rt.addLabel("Above Max Intensity (Cyan)");
        rt.addValue("Cutoff Value", maxThreshold);

        rt.incrementCounter();
        rt.addLabel("High Gradient (Blue)");
        rt.addValue("Cutoff Value", gradientCutoff);
                
        rt.show("Mask Cutoff Values: " + imageTitle);
    }
        
    /**
     * Logs analysis results to ImageJ log window
     */
    private void logResults(ImagePlus imp, NoiseStatisticsData result) {
    	IJ.log("=== Noise Statistics Analysis ===");
        IJ.log("Image: " + imp.getTitle());
        IJ.log("");
        IJ.log("Gradient Threshold: " + gradientThreshold);
        IJ.log("Display range: " + String.format("%.2f to %.2f", 
            result.getRangeDisplay()[0], result.getRangeDisplay()[1]));
        IJ.log("Analysis range: " + String.format("%.2f to %.2f", 
            result.getRangeAnalysis()[0], result.getRangeAnalysis()[1]));
        IJ.log("Peak intensity (I_peak): " + String.format("%.2f", result.getIPeak()));
        IJ.log("");
        IJ.log("--- Free Fit Results ---");
        IJ.log("Zero intercept (I0): " + String.format("%.2f", result.getI0()));
        IJ.log("Slope: " + String.format("%.2f", result.getSlope()));
        IJ.log("SNR <S^2>/<N^2>: " + String.format("%.2f", result.getSNR()));
        IJ.log("");
        IJ.log("--- Fit with Dark Count ---");
        IJ.log("Dark Count: " + String.format("%.2f", darkCount));
        IJ.log("Slope: " + String.format("%.2f", result.getSlopeHeader()));
        IJ.log("SNR1 <S^2>/<N^2>: " + String.format("%.2f", result.getSNR1()));
    }
       
	/** Creates a plot showing variance vs mean intensity with fitted lines */
    private Plot plotNoiseDistribution(
    		NoiseStatisticsData result,
            double darkCount,
            boolean showSNR,
            boolean showSNR1,
            String imageTitle) {
    	
    	double[] meanVals = result.getMeanValues();
    	double[] varVals = result.getVarianceValues();
    	double[] rangeAnalysis = result.getRangeAnalysis();
    	double slope = result.getSlope();
    	double i0 = result.getI0();
    	double snr = result.getSNR();
    	double snr1 = result.getSNR1();
    	double slopeHeader = result.getSlopeHeader();
        
        Plot plot = new Plot(
            "Noise Distribution: " + imageTitle,
            "Image Intensity Mean",
            "Image Intensity Variance"
        );
        StringBuilder legend = new StringBuilder();
        
        // Plot data points
        plot.setColor("BLUE");
        plot.addPoints(meanVals, varVals, Plot.CIRCLE);
        legend.append(String.format("Binned Data (Grad. threshold = %.3f)\n", gradientThreshold));
        
        double[] fitX, fitY;        
        double[] limits = new double[] {Arrays.stream(meanVals).min().getAsDouble(),
								Arrays.stream(meanVals).max().getAsDouble(),
								Arrays.stream(varVals).min().getAsDouble(),
								Arrays.stream(varVals).max().getAsDouble()};

        // Add threshold lines
		plot.setLineWidth(1);
        plot.setColor("red");
        plot.drawDottedLine(rangeAnalysis[0], limits[2], rangeAnalysis[0], limits[3], 1);
        plot.setColor("cyan");
        plot.drawDottedLine(rangeAnalysis[1], limits[2], rangeAnalysis[1], limits[3], 1);
        
        // Add free fit line if requested
		plot.setLineWidth(1.5f);
        if (showSNR) {
            fitX = new double[] {rangeAnalysis[0], rangeAnalysis[1]};
            double intercept = -slope * i0;
            fitY = new double[] {
                slope * fitX[0] + intercept,
                slope * fitX[1] + intercept
            };
			plot.setColor(new Color(170, 0, 255));
			plot.addPoints(fitX, fitY, Plot.LINE);
            legend.append(String.format(
                "Free fit: SNR = %.3f, I0 = %.3f\n",
                snr, i0
            ));            
        }
        
        // Add dark count fit line if requested
        if (showSNR1) {
            fitX = new double[] {rangeAnalysis[0], rangeAnalysis[1]};
            fitY = new double[] {
                slopeHeader * (fitX[0] - darkCount),
                slopeHeader * (fitX[1] - darkCount)
            };
            plot.setColor(new Color(0, 221, 0));
            plot.addPoints(fitX, fitY, Plot.LINE);
            legend.append(String.format(
                "Constrained fit: SNR1 = %.3f, I1 = %.3f\n",
                snr1, darkCount
            ));
        }

        plot.setColor("BLACK");
		plot.addLegend(legend.toString());
        plot.setLimits(limits);
        plot.setFrameSize(900,600);
        return plot;
    }

	private void getAllPersistedParams() {
		thrMinAnalysis = ParamPersister.get(imp, "N_thrMinAnalysis", 0.01);
		thrMaxAnalysis = ParamPersister.get(imp, "N_thrMaxAnalysis", 0.01);
		nbinsAnalysis = ParamPersister.get(imp, "N_nbinsAnalysis", 256);
        gradientThreshold = ParamPersister.get(imp, "N_gradientThreshold", 0.50);
		displaySNR = ParamPersister.get(imp, "N_displaySNR", true);
		darkCount = ParamPersister.get(imp, "N_darkCount", 0.0);
        displaySNR1 = ParamPersister.get(imp, "N_displaySNR1", false);
        showMaskVisualization = ParamPersister.get(imp, "N_showMaskVisualization", true);
        saveFigs = ParamPersister.get(imp, "N_saveFigs", false);
	}

	private void setAllPersistedParams() {
		ParamPersister.set(imp, "N_thrMinAnalysis", thrMinAnalysis);
		ParamPersister.set(imp, "N_thrMaxAnalysis", thrMaxAnalysis);
		ParamPersister.set(imp, "N_nbinsAnalysis", nbinsAnalysis);
        ParamPersister.set(imp, "N_gradientThreshold", gradientThreshold);
		ParamPersister.set(imp, "N_displaySNR", displaySNR);
		ParamPersister.set(imp, "N_darkCount", darkCount);
        ParamPersister.set(imp, "N_displaySNR1", displaySNR1);
        ParamPersister.set(imp, "N_showMaskVisualization", showMaskVisualization);
        ParamPersister.set(imp, "N_saveFigs", saveFigs);
		logParams();
	}

	private void logParams() {
		LinkedHashMap<String, Object> params = new LinkedHashMap<>();
		params.put("Analysis CDF threshold (lower)", thrMinAnalysis);
		params.put("Analysis CDF threshold (upper)", thrMaxAnalysis);
		params.put("Number of bins (analysis)", nbinsAnalysis);
		params.put("Gradient threshold", gradientThreshold);
		params.put("Display SNR0", displaySNR);
		params.put("Dark count", darkCount);
		params.put("Display SNR1", displaySNR1);
		params.put("Show mask visualization", showMaskVisualization);
		params.put("Save figures", saveFigs);
		ParamPersister.logParams("Parameters - Noise Statistics", params);
	}
}