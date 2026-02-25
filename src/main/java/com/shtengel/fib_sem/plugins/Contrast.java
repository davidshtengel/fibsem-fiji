package com.shtengel.fib_sem.plugins;

import java.awt.Checkbox;
import java.awt.Color;
import java.awt.TextField;
import java.util.LinkedHashMap;
import java.util.Vector;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.ContrastAnalyzer;
import com.shtengel.fib_sem.data.ContrastData;
import com.shtengel.fib_sem.data.ThresholdData;
import com.shtengel.fib_sem.util.FigBuilder;
import com.shtengel.fib_sem.util.ImageResolver;
import com.shtengel.fib_sem.util.ParamPersister;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.process.ImageProcessor;

/**
 * ImageJ plugin for computing image contrast from FIB-SEM intensity distributions.
 *
 * <p>Uses gradient filtering to select low-texture pixels, identifies
 * high- and low-intensity peaks via CDF thresholds, and reports the
 * contrast metric {@code (I_high - I_low) / (I_mean - I0)}.</p>
 *
 * <p>The dark count (I0) can be pre-computed by running the Noise Statistics plugin first,
 * or entered manually.</p>
 *
 * @see ContrastAnalyzer
 * @see ContrastData
 */
@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Contrast Analysis")
public class Contrast implements Command {
	
	private boolean ranSNR;
	private double i0;
	private double gradientThreshold;
	private double thrMinContrast;
	private double thrMaxContrast;
	private int nbins;
	private boolean saveFigs;
	private ImagePlus imp;

	@Override
	public void run() {
		imp = ImageResolver.resolveSourceImage();
		if (imp == null) {
			IJ.error("Contrast", "No image open.");
			return;
		}

		if(!showDialog()) {
			return;
		}

		ImageProcessor ip = imp.getProcessor();
		if (ip == null) {
			IJ.error("Contrast", "Could not access image processor");
			return;
		}

		// Apply the ROI to processor if present
		Roi roi = imp.getRoi();
		if (roi != null) {
            ip.setRoi(roi);
            String roiName = roi.getName() != null ? roi.getName() : "unnamed ROI";
            IJ.log("Processing ROI: " + roiName + " (type: " + roi.getTypeAsString() + ")");
        } else {
            IJ.log("Processing entire image (no ROI selected)");
        }

		// Compute contrast
        ContrastData result = ContrastAnalyzer.computeContrast(
            ip, i0, gradientThreshold,
            thrMinContrast, thrMaxContrast, nbins
        );

		// Log results
        logResults(result);

        // Create PDF plot
        Plot pdfPlot = plotContrastPDF(result, imp.getTitle());

        if (saveFigs) {
            FileInfo fi = imp.getOriginalFileInfo();
            String dir = (fi != null && fi.directory != null) ? fi.directory : "";
            String baseName = ImageResolver.getBaseName(imp);

            FigBuilder.createAndSave(pdfPlot, "Contrast PDF",
                                     dir + baseName + "_contrast_pdf.png");
        }

        pdfPlot.show();
        IJ.showStatus("Contrast analysis complete.");
	}

	private boolean showDialog() {
		getAllPersistedParams();

		GenericDialog gd = new GenericDialog("Contrast Analysis");

		if (ranSNR) {
			gd.addMessage("I0 was obtained from a prior SNR (Noise Statistics) analysis.\n"
						 + "It represents the intensity at zero variance (zero variance intercept).");
			gd.addNumericField("I0 (Pre-computed)", i0, 2, 10, "");
		} else {
			gd.addMessage("I0 (dark count) was NOT obtained from a prior SNR (Noise Statistics) analysis.\n"
						 + "It represents the intensity at zero variance (zero variance intercept).");
			gd.addCheckbox("Check to run 'Noise Statistics' first, then calculate contrast with loaded I0.", false);
			gd.addNumericField("I0 (Not pre-computed)", i0, 2, 10, "");
		
			@SuppressWarnings("unchecked")
			Vector<Checkbox> checkboxes = gd.getCheckboxes();
			@SuppressWarnings("unchecked")
			Vector<TextField> numericFields = gd.getNumericFields();

			Checkbox runSNRCheckbox = checkboxes.lastElement();
			TextField i0Field = numericFields.lastElement();

			gd.addDialogListener((dialog, e) -> {
				i0Field.setEnabled(!runSNRCheckbox.getState());
				return true;
			});
		}

		gd.addMessage("Gradient threshold controls what fraction of pixels (by gradient magnitude)\n"
                     + "are retained for analysis. E.g. 0.75 keeps 75% of pixels with lowest local gradients.");
        gd.addNumericField("Gradient threshold", gradientThreshold, 2, 6, "");

        gd.addMessage("CDF thresholds for identifying the two intensity peaks.\n"
                     + "Adjust these to coincide with the peaks of low- and high- intensity modes in the PDF.");
        gd.addNumericField("Lower threshold (I_low)", thrMinContrast, 3, 6, "");
        gd.addNumericField("Upper threshold (I_high)", thrMaxContrast, 3, 6, "");
        gd.addNumericField("Number of bins", nbins, 0, 6, "");

        gd.addMessage("Export:");
        gd.addCheckbox("Save plot as titled figure", saveFigs);

        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }

        // Retrieve values
		if (ranSNR) {
			i0 = gd.getNextNumber();
		} else {
			boolean runSNR = gd.getNextBoolean();
			i0 = gd.getNextNumber();

			if (runSNR) {
				IJ.run("Noise Statistics (Single Image)");  // blocks until plugin's dialog completes
				getAllPersistedParams();      // re-read i0 from image properties
				if (!ranSNR) {
					IJ.error("Contrast", "SNR analysis did not complete. Please enter I0 manually.");
					return false;
				}
			}
		}
        gradientThreshold = gd.getNextNumber();
        thrMinContrast = gd.getNextNumber();
        thrMaxContrast = gd.getNextNumber();
        nbins = (int) gd.getNextNumber();
        saveFigs = gd.getNextBoolean();

        // Validate
        if (gradientThreshold < 0 || gradientThreshold > 1) {
            IJ.error("Invalid gradient threshold. Must be between 0 and 1.");
            return false;
        }
        if (thrMinContrast < 0 || thrMinContrast > 1 || thrMaxContrast < 0 || thrMaxContrast > 1) {
            IJ.error("Invalid CDF thresholds. Must be between 0 and 1.");
            return false;
        }
        if (nbins < 10 || nbins > 10000) {
            IJ.error("Number of bins must be between 10 and 10000.");
            return false;
        }

        setAllPersistedParams();
        return true;
	}

	/**
     * Creates a plot showing the histogram (PDF) of the gradient-filtered
     * smoothed subset, with vertical threshold markers and contrast annotation.
     *
     * @param result     the contrast analysis results
     * @param imageTitle title of the source image
     * @return {@code Plot} ready for display
     */
    private Plot plotContrastPDF(ContrastData result, String imageTitle) {
        ThresholdData td = result.getThresholdData();
        double[] pdf = td.getPDF();
        int nBins = pdf.length;
        double iLow = result.getILow();
        double iHigh = result.getIHigh();
        double minInt = td.getMinIntensity();
        double maxInt = td.getMaxIntensity();
		double contrast = result.getContrast();

        // Compute bin centers
        double binWidth = (maxInt - minInt) / nBins;
        double[] binCenters = new double[nBins];
        for (int i = 0; i < nBins; i++) {
            binCenters[i] = minInt + (i + 0.5) * binWidth;
        }

        // Find max PDF for y-axis scaling
        double maxPdf = 0;
        for (double v : pdf) {
            if (v > maxPdf) maxPdf = v;
        }

        Plot plot = new Plot(
            "Contrast PDF: " + imageTitle,
            "Intensity",
            "Probability Density"
        );

        StringBuilder legend = new StringBuilder();

        // PDF curve
        plot.setColor(Color.BLUE);
        plot.addPoints(binCenters, pdf, Plot.LINE);
        legend.append("PDF\n");

        // Vertical threshold lines
        plot.setColor(Color.RED);
        plot.drawDottedLine(iLow, 0, iLow, maxPdf, 1);
		plot.addPoints(new double[] {iLow, iLow}, new double[] {0, maxPdf}, Plot.DOT);
		legend.append(String.format("I_low = %.2f\n", iLow));
        plot.setColor(Color.CYAN);
        plot.drawDottedLine(iHigh, 0, iHigh, maxPdf, 1);
		plot.addPoints(new double[] {iHigh, iHigh}, new double[] {0, maxPdf}, Plot.DOT);
        legend.append(String.format("I_high = %.2f\n", iHigh));

        // Annotations
        plot.setColor(Color.BLACK);
		plot.addLabel(0.05, 0.25, String.format("I0 = %.2f\nContrast = %.3f", i0, contrast));

        plot.addLegend(legend.toString());
        plot.setLimits(minInt, maxInt, 0, maxPdf * 1.05);

        return plot;
    }

	/**
     * Logs analysis results to the ImageJ log window.
     */
    private void logResults(ContrastData result) {
        IJ.log("=== Contrast Analysis ===");
        IJ.log("Image: " + imp.getTitle());
        IJ.log("");
        IJ.log("Gradient threshold: " + String.format("%.2f", result.getGradientThreshold()));
        IJ.log("Pixels analyzed: " + result.getSubsetPixels()
             + " / " + result.getTotalPixels()
             + String.format(" (%.1f%%)", 100.0 * result.getSubsetPixels() / result.getTotalPixels()));
        IJ.log("");
        IJ.log("I_low (CDF = " + String.format("%.3f", result.getThrMinContrast()) + "): "
             + String.format("%.2f", result.getILow()));
        IJ.log("I_high (CDF = " + String.format("%.3f", 1.0 - result.getThrMaxContrast()) + "): "
             + String.format("%.2f", result.getIHigh()));
        IJ.log("I_mean: " + String.format("%.2f", result.getIMean()));
        IJ.log("I0 (dark count): " + String.format("%.2f", result.getI0()));
        IJ.log("");
        IJ.log("Contrast = " + String.format("%.3f", result.getContrast()));
    }

	private void getAllPersistedParams() {
		ranSNR = ParamPersister.get(imp, "C_ranSNR", false);
		i0 = ParamPersister.get(imp, "C_i0", 0.0);
        gradientThreshold = ParamPersister.get(imp, "C_gradientThreshold", 0.50);
        thrMinContrast = ParamPersister.get(imp, "C_thrMinContrast", 0.25);
        thrMaxContrast = ParamPersister.get(imp, "C_thrMaxContrast", 0.25);
        nbins = ParamPersister.get(imp, "C_nbins", 256);
        saveFigs = ParamPersister.get(imp, 	"C_saveFigs", false);
	}

	private void setAllPersistedParams() {
		ParamPersister.set(imp, "C_ranSNR", ranSNR);
		ParamPersister.set(imp, "C_i0", i0);
        ParamPersister.set(imp, "C_gradientThreshold", gradientThreshold);
        ParamPersister.set(imp, "C_thrMinContrast", thrMinContrast);
        ParamPersister.set(imp, "C_thrMaxContrast", thrMaxContrast);
        ParamPersister.set(imp, "C_nbins", nbins);
        ParamPersister.set(imp, "C_saveFigs", saveFigs);
        logParams();
	}

	private void logParams() {
        LinkedHashMap<String, Object> params = new LinkedHashMap<>();
        params.put("Dark count (I0)", i0);
        params.put("Gradient threshold", gradientThreshold);
        params.put("Lower CDF threshold", thrMinContrast);
        params.put("Upper CDF threshold", thrMaxContrast);
        params.put("Number of bins", nbins);
        params.put("Save figures", saveFigs);
        ParamPersister.logParams("Parameters - Contrast Analysis", params);
    }
}