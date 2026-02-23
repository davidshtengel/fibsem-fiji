package com.shtengel.fib_sem.plugins;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import ij.process.ImageProcessor;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.Roi;
import java.awt.Color;
import java.util.Arrays;

import com.shtengel.fib_sem.core.EdgeTransitionAnalyzer;
import com.shtengel.fib_sem.data.EdgeTransitionData;
import com.shtengel.fib_sem.util.EllipseFitter;
import com.shtengel.fib_sem.util.FigBuilder;
import com.shtengel.fib_sem.util.ImageResolver;
import com.shtengel.fib_sem.util.ParamPersister;

@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Resolution; Edge Transition Analysis")
public class EdgeTransitions implements Command {
	
    private double lowerBound = 0.37;
	private double upperBound = 0.63;
	private double pixelSize = 1.0;
	private int subsetSize = 25;
	private double sectionLength = 25.0;
    private int minMaxAperture = 5;
	private float transitionLowLimit = 0.5f;
	private float transitionHighLimit = 10.0f;
	private int neighborExclusionRadius = 10;
    private boolean excludeCenter = false;
    private double centerExclusionRadius = 20;
	private float thrMinCriterion = 0.15f;
    private float thrMaxCriterion = 0.15f;
    private float gradientThreshold = 0.005f;
    private boolean saveFigs = false;
    private ImagePlus imp;
    	
	@Override
	public void run() {
		// Get active image
        imp = ImageResolver.resolveSourceImage();
        if (imp == null) {
            IJ.error("Edge Transition Analysis", "No image open.");
            return;
        }
        
        // Read pixel size from .dat header, if applicable
        double datPixelSize = getPixelSize(imp);
        if (datPixelSize > 0) {
        	pixelSize = datPixelSize;
        }

		// Show dialog
		if (!showDialog()) { return; }
		
		// Get the ImageProcessor
		ImageProcessor ip = imp.getProcessor();
		if (ip == null) {
            IJ.error("Edge Transition Analysis", "Could not access image processor.");
            return;
		}

		// Get the ROI if one exists
		Roi roi = imp.getRoi();
		int xOffset = 0;
		int yOffset = 0;

		if (roi != null) {
        	ip.setRoi(roi);
        	xOffset = roi.getBounds().x;
        	yOffset = roi.getBounds().y;
        	String roiName = roi.getName() != null ? roi.getName() : "unnamed ROI";
            IJ.log("Processing ROI: " + roiName + " (type: " + roi.getTypeAsString() + ")");
            IJ.log("ROI offset: (" + xOffset + ", " + yOffset + ")");
        } else {
        	IJ.log("Processing entire image (no ROI selected)");
		}

		IJ.showStatus("Analyzing edge transitions");

		// Analyze edge transitions
		EdgeTransitionData result = EdgeTransitionAnalyzer.analyzeEdgeTransitions(
			ip,
            lowerBound,
            upperBound,
            pixelSize,
            subsetSize,
            sectionLength,
            minMaxAperture,
            transitionLowLimit,
            transitionHighLimit,
            neighborExclusionRadius,
            thrMinCriterion,
            thrMaxCriterion,
            excludeCenter,
            centerExclusionRadius,
            gradientThreshold
        );

        logResults(imp, result);        
        createResultsTable(result, imp.getTitle());        
        
        ImagePlus edges = drawEdgePoints(imp, result, xOffset, yOffset);
    	Plot transition = plotDirDistribution(result, imp.getTitle());
    	Plot directions = plotTransitionHist(result, imp.getTitle());
    	Plot intensities = plotIntensityProfiles(result, imp.getTitle());

    	if (edges != null) edges.show();
    	if (transition != null) transition.show();
    	if (directions != null) directions.show();
    	if (intensities != null) intensities.show();
    	
    	if (saveFigs) {
    		// Get directory from the source image
    	    FileInfo fi = imp.getOriginalFileInfo();
    	    String dir = (fi != null && fi.directory != null) ? fi.directory : "";
    	    String baseName = ImageResolver.getBaseName(imp);
    		
    	    if (edges != null)	FigBuilder.createAndSave(edges.flatten(), 
    				"Image with Edge Points (colored by normalized transition length)",
    				dir + baseName + "_edge_points.png");
    		if (transition != null) FigBuilder.createAndSave(transition, 
    				String.format("%.2f to %.2f Transition Distribution", lowerBound, upperBound), 
    				dir + baseName + "_transition_distribution.png");
    		if (directions != null) FigBuilder.createAndSave(directions,
    				"Transition Distribution over Directions",
    				dir + baseName + "_transition_directional_distribution.png");
    		if (intensities != null) FigBuilder.createAndSave(intensities,
    				"Analyzed Transitions Intensity Profiles",
    				dir + baseName + "_intensity_profiles.png");
    	}
    	
        IJ.showStatus("Edge transition analysis complete.");
	}

	private boolean showDialog() {
		getAllPersistedParams();

		GenericDialog gd = new GenericDialog("Edge Transition Analysis");

		gd.addMessage("Transition bounds:");
	    gd.addNumericField("Lower bound", lowerBound, 2, 6, "");
	    gd.addNumericField("Upper bound", upperBound, 2, 6, "");

	    gd.addMessage("Image parameters:");	    
	    gd.addNumericField("Pixel size", pixelSize, 2, 6, "nm");
 	    gd.addNumericField("Subset size", subsetSize, 0, 6, "pixels");
	    gd.addNumericField("Section length", sectionLength, 0, 6, "pixels");
	    gd.addNumericField("Min/max aperture", minMaxAperture, 0, 6, "pixels");

	    gd.addMessage("Transition limits:");
	    gd.addNumericField("Transition low limit", transitionLowLimit, 1, 6, "pixels");
	    gd.addNumericField("Transition high limit", transitionHighLimit, 1, 6, "pixels");
	    
	    gd.addMessage("Exclusion parameters:");
	    gd.addNumericField("Neighbor exclusion radius", neighborExclusionRadius, 0, 6, "pixels");
	    gd.addCheckbox("Exclude center", excludeCenter);
	    gd.addNumericField("Center exclusion radius", centerExclusionRadius, 0, 6, "pixels");
	    
	    gd.addMessage("Threshold parameters:");
	    gd.addNumericField("Minimum threshold criterion", thrMinCriterion, 2, 6, "");
	    gd.addNumericField("Maximum threshold criterion", thrMaxCriterion, 2, 6, "");
	    gd.addNumericField("Gradient threshold", gradientThreshold, 4, 6, ""); // Might need to add clarifying comment

	    gd.addMessage("Export:");
	    gd.addCheckbox("Save image/plot as titled figures", saveFigs);

	    gd.showDialog();
		if (gd.wasCanceled()) { return false; }

		// Retrieve values
		lowerBound = gd.getNextNumber();
	    upperBound = gd.getNextNumber();
	    pixelSize = gd.getNextNumber();
	    subsetSize = (int) gd.getNextNumber();
	    sectionLength = gd.getNextNumber();
	    minMaxAperture = (int) gd.getNextNumber();
	    transitionLowLimit = (float) gd.getNextNumber();
	    transitionHighLimit = (float) gd.getNextNumber();
	    neighborExclusionRadius = (int) gd.getNextNumber();
	    excludeCenter = gd.getNextBoolean();
	    centerExclusionRadius = gd.getNextNumber();
	    thrMinCriterion = (float) gd.getNextNumber();
	    thrMaxCriterion = (float) gd.getNextNumber();
	    gradientThreshold = (float) gd.getNextNumber();
	    saveFigs = gd.getNextBoolean();

		// Validation
	    if (!validateBounds(lowerBound, upperBound)) {
	    	return false;
	    }
		if (!validateThresholds(thrMinCriterion, thrMaxCriterion)) {
			return false;
		}
		if (subsetSize < 10) {
			IJ.error("Invalid subset size. Must be at least 10.");
			return false;
		}
		if (minMaxAperture < 2) {
			IJ.error("Invalid min/max aperture value. Must be at least 2.");
			return false;
		}
		if (gradientThreshold < 0 || gradientThreshold > 1) {
			IJ.error("Invalid gradient threshold. Must be between 0 and 1.");
            return false;
		}

		setAllPersistedParams();
		
		return true;
	}
	
	public void logResults(ImagePlus imp, EdgeTransitionData result) {
	    IJ.log("=== Edge Transition Analysis Results ===");
	    IJ.log("Image: " + imp.getTitle());
	    IJ.log("Pixel size: "+ pixelSize);
	    IJ.log("Initial edge points detected: " + result.getInitialEdgeCount());
	    IJ.log("After neighbor exclusion: " + result.getFilteredEdgeCount());
	    IJ.log("Valid transitions analyzed: " + result.getValidCount());
	    IJ.log("");
	    
	    if (result.getValidCount() > 0) {
	        IJ.log("Transition Statistics:");
	        IJ.log(String.format("  Mean: %.3f pixels (%.3f nm)", 
	            result.getMeanTransition(), 
	            result.getMeanTransition() * pixelSize));
	        IJ.log(String.format("  Std Dev: %.3f pixels (%.3f nm)", 
	            result.getStdTransition(),
	            result.getStdTransition() * pixelSize));
	        IJ.log(String.format("  Min: %.3f pixels (%.3f nm)", 
	            result.getMinTransition(),
	            result.getMinTransition() * pixelSize));
	        IJ.log(String.format("  Max: %.3f pixels (%.3f nm)", 
	            result.getMaxTransition(),
	            result.getMaxTransition() * pixelSize));
	    } else {
	        IJ.log("No valid transitions found.");
	    }
	    IJ.log("========================================");
	}

	public void createResultsTable(EdgeTransitionData result, String imageTitle) {
	    if (result.getValidCount() == 0) {
	        IJ.log("No valid transitions to create results table.");
	        return;
	    }
	    
	    ij.measure.ResultsTable rt = new ij.measure.ResultsTable();
	    
	    int[] xSelected = result.getXSelected();
	    int[] ySelected = result.getYSelected();
	    double[] transitions = result.getTransitionDistances();
	    double[] cosX = result.getCosXSelected();
	    double[] cosY = result.getCosYSelected();
	    
	    for (int i = 0; i < xSelected.length; i++) {
	        rt.incrementCounter();
	        rt.addValue("X", xSelected[i]);
	        rt.addValue("Y", ySelected[i]);
	        rt.addValue("Transition (px)", transitions[i]);
	        rt.addValue("Transition (nm)", transitions[i] * pixelSize);
	        rt.addValue("Gradient_X", cosX[i]);
	        rt.addValue("Gradient_Y", cosY[i]);
	    }
	    
	    rt.show("Edge Transitions: " + imageTitle);
	}
	
	/**
     * Creates a visualization showing detected edge points on the original image
     */
	private ImagePlus drawEdgePoints(ImagePlus imp, EdgeTransitionData result, int xOffset, int yOffset) {
		double[] transitions = result.getTransitionDistances();
		int[] xSelected = result.getXSelected();
		int[] ySelected = result.getYSelected();

		if (xSelected.length == 0) {
			IJ.log("No valid edge points to display");
			return null;
		}
				
		ImagePlus impCopy = imp.duplicate();
		impCopy.setTitle("Visualized Edge Points: " + imp.getTitle());
        Overlay overlay = new Overlay();

        for (int i = 0; i < xSelected.length; i++) {
        	Color color = getTransitionColor(transitions[i], result);

        	PointRoi point = new PointRoi(xSelected[i] + xOffset, ySelected[i] + yOffset);
        	point.setPointType(2); 
        	point.setSize(2);
        	point.setStrokeColor(color);
        	
        	overlay.add(point);
        }
        
        impCopy.setOverlay(overlay);
        
        return impCopy;
	}

	/**
     * Creates a histogram plot of transition distances
     */
	private Plot plotTransitionHist(EdgeTransitionData result, String imageTitle) {
		double[] transitions = result.getTransitionDistances();
        
        if (transitions.length == 0) {
            IJ.log("No transitions to plot.");
            return null;
        }

        int nBins = (int) (4 * Math.cbrt(transitions.length));	// rice rule for bin size
        double minVal = result.getMinTransition();
        double maxVal = result.getMaxTransition();
        double binWidth = (maxVal - minVal) / nBins;
        
        int[] hist = new int[nBins];
        double[] binCenters = new double[nBins];
        
        for (int i = 0; i < nBins; i++) {
            binCenters[i] = minVal + (i + 0.5) * binWidth;
        }
        
        for (double trans : transitions) {
            int binIdx = (int) ((trans - minVal) / binWidth);
            if (binIdx >= nBins) binIdx = nBins - 1;
            if (binIdx < 0) binIdx = 0;
            hist[binIdx]++;
        }

        double[] histDouble = new double[nBins];
        for (int i = 0; i < nBins; i++) {
            histDouble[i] = hist[i];
        }
        
        String plotTitle = String.format(
        	"Transition Distribution [%.3f, %.3f]: %s",
        	lowerBound,
        	upperBound,
        	imageTitle
        );
        Plot plot = new Plot(plotTitle, "Transition Distance (pixels)", "Count");
        
		// Add mean line
        plot.setColor(Color.RED);
        plot.drawLine(result.getMeanTransition(), 0, result.getMeanTransition(), Arrays.stream(histDouble).max().getAsDouble());
        
        // Add histogram
        plot.setLineWidth(1.5f);
		plot.setColor(Color.BLUE);
        plot.addPoints(binCenters, histDouble, Plot.BAR);        
        
        // Add statistics labels
        plot.setColor(Color.BLACK);
        plot.addLabel(0.025, 0.05, String.format("N = %d", transitions.length));
        plot.addLabel(0.025, 0.075, String.format("Mean = %.3f px", result.getMeanTransition()));
        plot.addLabel(0.025, 0.1, String.format("Std = %.3f px", result.getStdTransition()));
        plot.addLabel(0.025, 0.125, String.format("Pixel size = %.3f nm", pixelSize));
        plot.addLabel(0.025, 0.15, String.format("Mean = %.3f nm", result.getMeanTransition() * pixelSize));
        plot.addLabel(0.025, 0.175, String.format("Std = %.3f nm", result.getStdTransition() * pixelSize));
        plot.setFrameSize(800,800);
        return plot;
	}
	
	/**
	 * Creates a scatter plot showing transition distribution over directions
	 * (X and Y components of transition vectors) with ellipse fit
	 */
	private Plot plotDirDistribution(EdgeTransitionData result, String imageTitle) {
	    double[] transitions = result.getTransitionDistances();
	    double[] cosX = result.getCosXSelected();
	    double[] cosY = result.getCosYSelected();
	    
	    if (transitions.length == 0) {
	        IJ.log("No transitions to plot.");
	        return null;
	    }
	    
	    // Calculate X and Y components of transition vectors
	    double[] transX = new double[transitions.length];
	    double[] transY = new double[transitions.length];
	    for (int i = 0; i < transitions.length; i++) {
	        transX[i] = cosX[i] * transitions[i];
	        transY[i] = cosY[i] * transitions[i];
	    }
	    
	    // Find the range for symmetric axes
	    double maxAbsX = Math.max(Math.abs(Arrays.stream(transX).min().getAsDouble()),
	                               Math.abs(Arrays.stream(transX).max().getAsDouble()));
	    double maxAbsY = Math.max(Math.abs(Arrays.stream(transY).min().getAsDouble()),
	                               Math.abs(Arrays.stream(transY).max().getAsDouble()));
	    double maxAbs = Math.max(maxAbsX, maxAbsY);
	    
	    // Create the plot
	    String plotTitle = String.format("Transition Distribution over Directions: %s", imageTitle);
	    Plot plot = new Plot(plotTitle, "Transition X Component (pixels)", "Transition Y Component (pixels)");
	    plot.setLimits(-maxAbs * 1.1, maxAbs * 1.1, -maxAbs * 1.1, maxAbs * 1.1);
	    plot.setFrameSize(800,800);
	    plot.setColor(Color.LIGHT_GRAY);
	    plot.drawLine(-maxAbs * 1.1, 0, maxAbs * 1.1, 0);  // X axis
	    plot.drawLine(0, -maxAbs * 1.1, 0, maxAbs * 1.1);  // Y axis
	    
	    // Fit ellipse
	    try {
	        double[] ellipseFit = EllipseFitter.fitCenteredEllipse(cosX, cosY, transitions);
	        
	        if (ellipseFit != null) {
	            double a = ellipseFit[0];
	            double b = ellipseFit[1];
	            double phi = ellipseFit[2];
	            
	            // Generate ellipse curve
	            int nPoints = 360;
	            double[] ellipseX = new double[nPoints];
	            double[] ellipseY = new double[nPoints];
	            
	            for (int i = 0; i < nPoints; i++) {
	            	double theta = -Math.PI + (2.0 * Math.PI * i) / nPoints;
	                double thetaRot = theta - phi;
	                double cosT = Math.cos(thetaRot);
	                double sinT = Math.sin(thetaRot);
	                
	                double term1 = cosT / a;
	                double term2 = sinT / b;
	                double r = 1.0 / Math.sqrt(term1 * term1 + term2 * term2);
	                
	                ellipseX[i] = r * Math.cos(theta);
	                ellipseY[i] = r * Math.sin(theta);
	            }
	            
	            // Draw fitted ellipse
	            plot.setColor(Color.MAGENTA);
	            plot.setLineWidth(2);
	            plot.addPoints(ellipseX, ellipseY, Plot.LINE);
	            double ellipticity = EllipseFitter.computeEllipticity(a, b);
	            plot.addLegend(String.format(
	            		"Ellipticity: %.3f",
	            		ellipticity)
            	);
	        } else {
	            IJ.log("Could not fit ellipse to transition data");
	        }
	    } catch (Exception e) {
	        IJ.log("Error fitting ellipse: " + e.getMessage());
	        e.printStackTrace();
	    }
	    
 		// Add scatter points
 		plot.setLineWidth(1);
	    for(int i = 0; i < transitions.length; i++) {
        	Color color = getTransitionColor(transitions[i], result);
        	plot.setColor(color);
        	plot.addPoints(new double[] {transX[i]}, new double[] {transY[i]}, Plot.CIRCLE);
	    }
	    plot.setLineWidth(1.5f);
	    return plot;
	}
	
	/**
	 * Creates a plot showing intensity profiles for all analyzed edge transitions.
	 */
	private Plot plotIntensityProfiles(EdgeTransitionData result, String imageTitle) {
		double[][] imgValsAll = result.getImageValuesAll();
		
		if(imgValsAll == null || imgValsAll.length == 0) {
			IJ.log("No intensity profiles to plot");
			return null;
		}
		
		int nProfiles = imgValsAll.length;
		int profileLength = imgValsAll[0].length;
		
		// Create distance array (x-axis), centered around 0 and converted to nm
		double[] distNm = new double[profileLength];
		int halfLength = profileLength / 2;
		for (int i = 0; i < profileLength; i++) {
			distNm[i] = (i - halfLength + 1) * pixelSize;
		}
		
		// Find global min/max for y-axis limits
		double minIntensity = Double.MAX_VALUE;
		double maxIntensity = -Double.MAX_VALUE;
		for (double[] profile : imgValsAll) {
			for (double value : profile) {
				if (value < minIntensity) minIntensity = value;
				if (value > maxIntensity) maxIntensity = value;
			}
		}

		// Create the plot
		String plotTitle = String.format("Intensity Profiles: %s", imageTitle);
		Plot plot = new Plot(plotTitle, "Distance (nm)", "Interpolated Image Intensity");
		
		// Set plot limits w/ padding
		double xMin = distNm[0] - pixelSize;
		double xMax = distNm[profileLength - 1] + pixelSize;
		double yMax = maxIntensity + (maxIntensity - minIntensity) * 0.05;
		double yMin = minIntensity - (maxIntensity - minIntensity) * 0.05;
		plot.setLimits(xMin, xMax, yMin, yMax);
		plot.setFrameSize(600, 600);
		
		// Plot each profile w/ coloring
		for (int j = 0; j < nProfiles; j++) {
			Color color = getColor((double)(nProfiles - j) / nProfiles);
			plot.setColor(color);
			plot.setLineWidth(0.5f);
			plot.addPoints(distNm, imgValsAll[j], Plot.LINE);
		}
		
		return plot;
	}

	// Helper Functions

	/** Extracts pixel size from image info for .dat files */
	private double getPixelSize(ImagePlus imp) {
		String info = imp.getInfoProperty();
		if (info == null || info.isEmpty()) {		
			return -1.0;
		}
		
		// Look for pattern ""Pixel Size: X.XX nm"
		String[] lines = info.split("\n");
		for (String line : lines) {
			line = line.trim();
			
			if (line.startsWith("Pixel Size:")) {
	            try {
	                String valueStr = line.substring("Pixel Size:".length()).trim();
	                
	                if (valueStr.endsWith(" nm")) {
	                    valueStr = valueStr.substring(0, valueStr.length() - 3).trim();
	                }
	                double value = Double.parseDouble(valueStr);
	                
	                if (value > 0 && value < 1000) { // Sanity check
	                    return value;
	                }
	            } catch (NumberFormatException | IndexOutOfBoundsException e) {
	                IJ.log("Warning: Found 'Pixel Size:' line but could not parse value: " + line);
	            }
	        }
		}
		return -1.0;
	}

	private void getAllPersistedParams() {
		lowerBound = ParamPersister.get(imp, "R-lowerBound", 0.37);
		upperBound = ParamPersister.get(imp, "R-upperBound", 0.63);
		pixelSize = ParamPersister.get(imp, "R-pixelSize", pixelSize);
        subsetSize = ParamPersister.get(imp, "R-subsetSize", 25);
		sectionLength = ParamPersister.get(imp, "R-sectionLength", 25.0);
		minMaxAperture = ParamPersister.get(imp, "R-minMaxAperture", 5);
        transitionLowLimit = ParamPersister.get(imp, "R-transitionLowLimit", 0.5f);
        transitionHighLimit = ParamPersister.get(imp, "R-transitionHighLimit", 10.0f);
        neighborExclusionRadius = ParamPersister.get(imp, "R-neighborExclusionRadius", 20);
        excludeCenter = ParamPersister.get(imp, "R-excludeCenter", false);
        centerExclusionRadius = ParamPersister.get(imp, "R-centerExclusionRadius", 10);
		thrMinCriterion = ParamPersister.get(imp, "R-thrMinCriterion", 0.15f);
        thrMaxCriterion = ParamPersister.get(imp, "R-thrMaxCriterion", 0.15f);
        gradientThreshold = ParamPersister.get(imp, "R-gradientThreshold", 0.005f);
        saveFigs = ParamPersister.get(imp, "R-saveFigs", false);
	}

	private void setAllPersistedParams() {
		ParamPersister.set(imp, "R-lowerBound", lowerBound);
		ParamPersister.set(imp, "R-upperBound", upperBound);
		ParamPersister.set(imp, "R-pixelSize", pixelSize);
        ParamPersister.set(imp, "R-subsetSize", subsetSize);
		ParamPersister.set(imp, "R-sectionLength", sectionLength);
		ParamPersister.set(imp, "R-minMaxAperture", minMaxAperture);
        ParamPersister.set(imp, "R-transitionLowLimit", transitionLowLimit);
        ParamPersister.set(imp, "R-transitionHighLimit", transitionHighLimit);
        ParamPersister.set(imp, "R-neighborExclusionRadius", neighborExclusionRadius);
        ParamPersister.set(imp, "R-excludeCenter", excludeCenter);
        ParamPersister.set(imp, "R-centerExclusionRadius", centerExclusionRadius);
		ParamPersister.set(imp, "R-thrMinCriterion", thrMinCriterion);
        ParamPersister.set(imp, "R-thrMaxCriterion", thrMaxCriterion);
        ParamPersister.set(imp, "R-gradientThreshold", gradientThreshold);
        ParamPersister.set(imp, "R-saveFigs", saveFigs);
		logParams();
	}

	private void logParams() {
    IJ.log("--- Parameters - Edge Transitions ---");
    IJ.log("Lower bound: " + lowerBound);
    IJ.log("Upper bound: " + upperBound);
    IJ.log("Pixel size: " + pixelSize + " nm");
    IJ.log("Subset size: " + subsetSize + " px");
    IJ.log("Section length: " + sectionLength + " px");
    IJ.log("Min/max aperture: " + minMaxAperture + " px");
    IJ.log("Transition low limit: " + transitionLowLimit + " px");
    IJ.log("Transition high limit: " + transitionHighLimit + " px");
    IJ.log("Neighbor exclusion radius: " + neighborExclusionRadius + " px");
    IJ.log("Exclude center: " + excludeCenter);
    IJ.log("Center exclusion radius: " + centerExclusionRadius + " px");
    IJ.log("Min threshold criterion: " + thrMinCriterion);
    IJ.log("Max threshold criterion: " + thrMaxCriterion);
    IJ.log("Gradient threshold: " + gradientThreshold);
    IJ.log("Save figures: " + saveFigs);
    IJ.log("-------------------------------------");
}

	private boolean validateBounds(double lower, double upper) {
		if (lower < 0 || lower > 1 || upper < 0 || upper > 1) {
	        IJ.error("Invalid bound values. Must be between 0 and 1.");
	        return false;
	    }
	    if (lower >= upper) {
	        IJ.error("Lower bound must be less than upper bound.");
	        return false;
	    }
	    return true;
	}
	
	private boolean validateThresholds(float thrMin, float thrMax) {
	    if (thrMin < 0 || thrMin > 0.5 || thrMax < 0 || thrMax > 0.5) {
	        IJ.error("Invalid threshold values. Must be between 0 and 0.5.");
	        return false;
	    }
	    return true;
	}

	private Color getTransitionColor(double transitionValue, EdgeTransitionData result) {
	    double range = result.getMaxTransition() - result.getMinTransition();
	    double normalized = (transitionValue - result.getMinTransition()) / range;
	    return getColor(normalized);
	}
	
	/** Generate a rainbow color (similar to matplotlib's gist_rainbow_r colormap) */
	private Color getColor(double t) {
		float hue = (float)(0.8 - t * 0.8);
		if (hue < 0) hue += 1.0f;
		return Color.getHSBColor(hue, 0.9f, 0.9f);
	}
}