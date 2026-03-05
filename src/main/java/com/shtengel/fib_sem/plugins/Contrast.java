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
import com.shtengel.fib_sem.util.DoubleGaussianFitter;
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
	private boolean autoMode;
	private double iLow;
	private double iHigh;
	private int nbins;
	private boolean showComponents;
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

		// Compute contrast (auto or manual)
		ContrastData result;
		if (autoMode) {
			result = ContrastAnalyzer.computeContrastAuto(
				ip, i0, gradientThreshold, nbins
			);
		} else {
			result = ContrastAnalyzer.computeContrastManual(
				ip, i0, gradientThreshold, iLow, iHigh, nbins
			);
		}

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

        gd.addMessage("Peak identification for I_low and I_high:");
		gd.addCheckbox("Determine I_low and I_high automatically (double-Gaussian fit)", autoMode);
		gd.addCheckbox("Show individual Gaussian components on plot", showComponents);

		gd.addNumericField("I_low (lower intensity peak)", iLow, 2, 10, "");
		gd.addNumericField("I_high (upper intensity peak)", iHigh, 2, 10, "");
		gd.addNumericField("Number of bins", nbins, 0, 6, "");

		@SuppressWarnings("unchecked")
		Vector<Checkbox> allCheckboxes = gd.getCheckboxes();
		@SuppressWarnings("unchecked")
		Vector<TextField> allNumericFields = gd.getNumericFields();
		int nCheckboxes = allCheckboxes.size();
		Checkbox autoCheckbox = allCheckboxes.get(nCheckboxes - 2);
		Checkbox componentsCheckbox = allCheckboxes.get(nCheckboxes - 1);

		int nFields = allNumericFields.size();
		TextField iLowField  = allNumericFields.get(nFields - 3);
		TextField iHighField = allNumericFields.get(nFields - 2);
		iLowField.setEnabled(!autoMode);
		iHighField.setEnabled(!autoMode);
		componentsCheckbox.setEnabled(autoMode);

		gd.addDialogListener((dialog, e) -> {
			boolean auto = autoCheckbox.getState();
			iLowField.setEnabled(!auto);
			iHighField.setEnabled(!auto);
			componentsCheckbox.setEnabled(auto);
			return true;
		});

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
				IJ.run("Noise Statistics Analysis (Single Image)");  // blocks until plugin's dialog completes
				getAllPersistedParams();      // re-read i0 from image properties
				if (!ranSNR) {
					IJ.error("Contrast", "SNR analysis did not complete. Please enter I0 manually.");
					return false;
				}
			}
		}
        gradientThreshold = gd.getNextNumber();
        autoMode = gd.getNextBoolean();
        showComponents = gd.getNextBoolean();
		iLow  = gd.getNextNumber();
		iHigh = gd.getNextNumber();
        nbins = (int) gd.getNextNumber();
        saveFigs = gd.getNextBoolean();

        // Validate
        if (gradientThreshold < 0 || gradientThreshold > 1) {
            IJ.error("Invalid gradient threshold. Must be between 0 and 1.");
            return false;
        }
        if (!autoMode) {
			if (Double.isNaN(iLow) || Double.isNaN(iHigh)) {
				IJ.error("I_low and I_high must be valid numbers.");
				return false;
			}
			if (iLow >= iHigh) {
				IJ.error("I_low must be less than I_high.");
				return false;
			}
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
		plot.setLineWidth(1.5f);
        plot.addPoints(binCenters, pdf, Plot.LINE);
        legend.append("PDF\n");

        // Fitted double-gaussian curve (auto mode only)
        if (result.isAutoMode() && result.getGaussianFit() != null) {
			plot.setLineWidth(0.75f);
            DoubleGaussianFitter.FitResult fit = result.getGaussianFit();
            double[] fitY = new double[nBins];
            for (int i = 0; i < nBins; i++) {
                fitY[i] = fit.evaluate(binCenters[i]);
            }
            plot.setColor(Color.MAGENTA);
            plot.addPoints(binCenters, fitY, Plot.LINE);
            legend.append(String.format("Double-Gaussian Fit (R\u00B2=%.4f)\n", fit.rSquared()));

            // Individual Gaussian components
            if (showComponents) {
                double a1 = fit.getA1(), mu1 = fit.getMu1(), s1 = fit.getSigma1();
                double a2 = fit.getA2(), mu2 = fit.getMu2(), s2 = fit.getSigma2();
                double[] g1 = new double[nBins];
                double[] g2 = new double[nBins];
                for (int i = 0; i < nBins; i++) {
                    double dx1 = binCenters[i] - mu1;
                    double dx2 = binCenters[i] - mu2;
                    g1[i] = a1 * Math.exp(-0.5 * dx1 * dx1 / (s1 * s1));
                    g2[i] = a2 * Math.exp(-0.5 * dx2 * dx2 / (s2 * s2));
                }
                plot.setColor(Color.RED);
                plot.addPoints(binCenters, g1, Plot.LINE);
                legend.append(String.format("Gaussian 1 (\u03BC=%.2f)\n", mu1));
                plot.setColor(Color.CYAN);
                plot.addPoints(binCenters, g2, Plot.LINE);
                legend.append(String.format("Gaussian 2 (\u03BC=%.2f)\n", mu2));
            }
        }

        // Vertical threshold lines
		plot.setLineWidth(1);
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
		plot.addLabel(0.7, 0.25, String.format("I0 = %.2f\nContrast = %.3f", i0, contrast));

        plot.addLegend(legend.toString());
		plot.setFrameSize(800,600);
        plot.setLimits(minInt, maxInt, 0, maxPdf * 1.05);
		plot.setLineWidth(1.5f);
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
        IJ.log("I_low:  " + String.format("%.2f", result.getILow()));
		IJ.log("I_high: " + String.format("%.2f", result.getIHigh()));
		IJ.log("I_mean: " + String.format("%.2f", result.getIMean()));
        IJ.log("I0 (dark count): " + String.format("%.2f", result.getI0()));
        IJ.log("");
        IJ.log("Contrast = " + String.format("%.3f", result.getContrast()));
    }

	private void getAllPersistedParams() {
		ranSNR = ParamPersister.get(imp, "C_ranSNR", false);
		i0 = ParamPersister.get(imp, "C_i0", 0.0);
        gradientThreshold = ParamPersister.get(imp, "C_gradientThreshold", 0.25);
        autoMode = ParamPersister.get(imp, "C_autoMode", true);
		iLow  = ParamPersister.get(imp, "C_iLow", 0.0);
		iHigh = ParamPersister.get(imp, "C_iHigh", 0.0);
		nbins = ParamPersister.get(imp, "C_nbins", 256);
        saveFigs = ParamPersister.get(imp, 	"C_saveFigs", false);
	}

	private void setAllPersistedParams() {
		ParamPersister.set(imp, "C_ranSNR", ranSNR);
		ParamPersister.set(imp, "C_i0", i0);
        ParamPersister.set(imp, "C_gradientThreshold", gradientThreshold);
        ParamPersister.set(imp, "C_autoMode", autoMode);
		ParamPersister.set(imp, "C_iLow", iLow);
		ParamPersister.set(imp, "C_iHigh", iHigh);
		ParamPersister.set(imp, "C_nbins", nbins);
        ParamPersister.set(imp, "C_saveFigs", saveFigs);
        logParams();
	}

	private void logParams() {
        LinkedHashMap<String, Object> params = new LinkedHashMap<>();
        params.put("Dark count (I0)", i0);
        params.put("Gradient threshold", gradientThreshold);
        params.put("Auto mode", autoMode);
		if (!autoMode) {
			params.put("I_low (manual)", iLow);
			params.put("I_high (manual)", iHigh);
		}
		params.put("Number of bins", nbins);
        params.put("Save figures", saveFigs);
        ParamPersister.logParams("Parameters - Contrast Analysis", params);
    }
}