package com.shtengel.fib_sem.plugins;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.LinkedHashMap;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.Timer;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.ContrastAnalyzer;
import com.shtengel.fib_sem.util.Col;
import com.shtengel.fib_sem.data.ContrastData;
import com.shtengel.fib_sem.data.ThresholdData;
import com.shtengel.fib_sem.util.DoubleGaussianFitter;
import com.shtengel.fib_sem.util.FigBuilder;
import com.shtengel.fib_sem.util.ImageResolver;
import com.shtengel.fib_sem.util.LinkedSliderField;
import com.shtengel.fib_sem.util.OverlayTint;
import com.shtengel.fib_sem.util.ParamPersister;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

/**
 * ImageJ plugin for computing image contrast from FIB-SEM intensity distributions.
 *
 * <p>Uses gradient filtering to select low-texture pixels, identifies
 * high- and low-intensity peaks via CDF thresholds, and reports the
 * contrast metric {@code (I_high - I_low) / (I_mean - I0)}.</p>
 * 
 * <p>Provides an interactive preview window with sliders for I-Low, I-High,
 * and gradient threshold. Excluded pixels are visualized with color overlays.
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
	private boolean showGradientPreview;
	private int iLowWindow;
	private int iHighWindow;
	private ImagePlus imp;
	private ImageProcessor croppedIp;
	private int cropWidth, cropHeight;

	// Cached intermediate data (computed once on open)
	private float[] cachedGradMagnitudes;
	private float[] cachedSmoothedPixels;
	private float minIntensity, maxIntensity;

	// Preview (allocated once, reused across updatePreview() calls)
	private ImagePlus previewImp;
	private int[] previewRgb;
	private ColorProcessor previewCp;

	// UI components
	private JFrame frame;
	private LinkedSliderField gradientSliderField, iLowSliderField, iHighSliderField;
	private JCheckBox autoModeCheckbox, showComponentsCheckbox, saveFigsCheckbox, showGradientCheckbox;
	private JTextField i0Field, nbinsField, iLowWindowField, iHighWindowField;
	private JButton applyButton;
	private Timer updateTimer;

	// Pre-sorted interior gradients for timely lookup
	private float[] sortedInteriorGradients;

	// Track whether auto-fit values have been populated
	private boolean autoFitPopulated = false;

	@Override
	public void run() {
		imp = ImageResolver.resolveSourceImage();
		if (imp == null) {
			IJ.error("Contrast", "No image open.");
			return;
		}

		getAllPersistedParams();

		ImageProcessor ip = imp.getProcessor();
		if (ip == null) {
			IJ.error("Contrast", "Could not access image processor");
			return;
		}

		Roi roi = imp.getRoi();
		if (roi != null) {
			ip.setRoi(roi);
			croppedIp = ip.crop();
			String roiName = roi.getName() != null ? roi.getName() : "unnamed ROI";
			IJ.log("Preview using ROI: " + roiName + " (type: " + roi.getTypeAsString() + ")");
		} else {
			croppedIp = ip;
			IJ.log("Preview using entire image (no ROI selected)");
		}
		cropWidth = croppedIp.getWidth();
		cropHeight = croppedIp.getHeight();

		// Cache gradient magnitudes and smoothed pixels (one-time computation)
		IJ.showStatus("Computing gradient map for preview...");
		cachedGradMagnitudes = ContrastAnalyzer.computeGradientMagnitudes(croppedIp);
		cachedSmoothedPixels = ContrastAnalyzer.computeSmoothedPixels(croppedIp);

		// Determine intensity range from smoothed pixels
		minIntensity = Float.MAX_VALUE;
		maxIntensity = -Float.MAX_VALUE;
		for (float v : cachedSmoothedPixels) {
			if (v < minIntensity) minIntensity = v;
			if (v > maxIntensity) maxIntensity = v;
		}

		int interiorCount = (cropWidth - 2) * (cropHeight - 2);
		sortedInteriorGradients = new float[interiorCount];
		int k = 0;
		for (int y = 1; y < cropHeight - 1; y++) {
			for (int x = 1; x < cropWidth - 1; x++) {
				sortedInteriorGradients[k++] = cachedGradMagnitudes[y * cropWidth + x];
			}
		}
		java.util.Arrays.sort(sortedInteriorGradients);

		// Allocate reusable preview buffers
		previewRgb = new int[cropWidth * cropHeight];
		previewCp = new ColorProcessor(cropWidth, cropHeight, previewRgb);

		// Pre-fill border pixels as grayscale (never changes)
		float range = maxIntensity - minIntensity;
		if (range == 0) range = 1;
		for (int x = 0; x < cropWidth; x++) {
			int topIdx = x;
			int botIdx = (cropHeight - 1) * cropWidth + x;
			previewRgb[topIdx] = OverlayTint.grayscale(
				OverlayTint.toGray(cachedSmoothedPixels[topIdx], minIntensity, range));
			previewRgb[botIdx] = OverlayTint.grayscale(
				OverlayTint.toGray(cachedSmoothedPixels[botIdx], minIntensity, range));
		}
		for (int y = 1; y < cropHeight - 1; y++) {
			int leftIdx = y * cropWidth;
			int rightIdx = y * cropWidth + cropWidth - 1;
			previewRgb[leftIdx] = OverlayTint.grayscale(
				OverlayTint.toGray(cachedSmoothedPixels[leftIdx], minIntensity, range));
			previewRgb[rightIdx] = OverlayTint.grayscale(
				OverlayTint.toGray(cachedSmoothedPixels[rightIdx], minIntensity, range));
		}

		// Create preview image with the reusable processor
		previewImp = new ImagePlus("Contrast Preview: " + imp.getTitle(), previewCp);
		previewImp.show();

		// Build and show the interactive UI
		SwingUtilities.invokeLater(() -> {
			buildFrame();
			updatePreview();
		});
	}

	/** Constructs and displays the non-modal control JFrame with all panels. */
	private void buildFrame() {
		frame = new JFrame("Contrast Analysis");
		frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		frame.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				onClose();
			}
		});

		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
		mainPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

		mainPanel.add(buildI0Panel());
		mainPanel.add(Box.createVerticalStrut(8));
		mainPanel.add(buildSlidersPanel());
		mainPanel.add(Box.createVerticalStrut(8));
		mainPanel.add(buildOptionsPanel());
		mainPanel.add(Box.createVerticalStrut(8));

		// Apply/Close button panel
		JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
		applyButton = new JButton("Apply");
		applyButton.addActionListener(e -> onApply());
		buttonPanel.add(applyButton);
		JButton close = new JButton("Close");
		close.addActionListener(e -> onClose());
		buttonPanel.add(close);
		mainPanel.add(buttonPanel);

		frame.getContentPane().add(mainPanel, BorderLayout.CENTER);
		frame.pack();
		frame.setMinimumSize(new Dimension(380, 0));
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);

		// Debounce timer for preview updates (50ms)
		updateTimer = new Timer(50, e -> updatePreview());
		updateTimer.setRepeats(false);
	}

	/** Builds the I0 (dark count) input panel with optional SNR button. */
	private JPanel buildI0Panel() {
		JPanel panel = new JPanel(new GridBagLayout());
		panel.setBorder(BorderFactory.createTitledBorder("Dark Count (I0)"));
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(2, 4, 2, 4);
		gbc.anchor = GridBagConstraints.WEST;

		if (ranSNR) {
			gbc.gridx = 0; gbc.gridy = 0;
			panel.add(new JLabel("I0 (from Noise Statistics):"), gbc);
			gbc.gridx = 1;
			i0Field = new JTextField(String.format("%.2f", i0), 8);
			panel.add(i0Field, gbc);
		} else {
			gbc.gridx = 0; gbc.gridy = 0; gbc.gridwidth = 2;
			panel.add(new JLabel("<html><i>I0 not pre-computed. Enter manually or run Noise Statistics.</i></html>"), gbc);

			gbc.gridy = 1; gbc.gridwidth = 1;
			gbc.gridx = 0;
			panel.add(new JLabel("I0:"), gbc);
			gbc.gridx = 1;
			i0Field = new JTextField(String.format("%.2f", i0), 8);
			panel.add(i0Field, gbc);

			gbc.gridy = 2; gbc.gridx = 0; gbc.gridwidth = 2;
			JButton runSNRButton = new JButton("Run Noise Statistics...");
			runSNRButton.addActionListener(e -> onRunSNR());
			panel.add(runSNRButton, gbc);
		}

		return panel;
	}

	/** Builds the threshold sliders panel using LinkedSliderField components. */
	private JPanel buildSlidersPanel() {
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		panel.setBorder(BorderFactory.createTitledBorder("Preview Thresholds"));

		gradientSliderField = new LinkedSliderField(
			"Gradient threshold:", 0.0, 1.0, gradientThreshold, "%.3f");
		gradientSliderField.addChangeCallback(this::schedulePreviewUpdate);
		panel.add(gradientSliderField);

		showGradientCheckbox = new JCheckBox("Show gradient threshold in preview", showGradientPreview);
		showGradientCheckbox.addActionListener(e -> schedulePreviewUpdate());
		panel.add(showGradientCheckbox);

		iLowSliderField = new LinkedSliderField(
			"I-Low:", minIntensity, maxIntensity, iLow, "%.1f");
		iLowSliderField.addChangeCallback(this::schedulePreviewUpdate);
		panel.add(iLowSliderField);

		JPanel iLowWindowPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
		iLowWindowPanel.add(new JLabel("    Window:"));
		iLowWindowField = new JTextField(String.valueOf(iLowWindow), 4);
		iLowWindowField.addActionListener(e -> schedulePreviewUpdate());
		iLowWindowPanel.add(iLowWindowField);
		panel.add(iLowWindowPanel);

		iHighSliderField = new LinkedSliderField(
			"I-High:", minIntensity, maxIntensity, iHigh, "%.1f");
		iHighSliderField.addChangeCallback(this::schedulePreviewUpdate);
		panel.add(iHighSliderField);

		JPanel iHighWindowPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 4, 0));
		iHighWindowPanel.add(new JLabel("    Window:"));
		iHighWindowField = new JTextField(String.valueOf(iHighWindow), 4);
		iHighWindowField.addActionListener(e -> schedulePreviewUpdate());
		iHighWindowPanel.add(iHighWindowField);
		panel.add(iHighWindowPanel);

		// Set initial enabled state based on auto mode
		iLowSliderField.setEnabled(!autoMode);
		iHighSliderField.setEnabled(!autoMode);
		iLowWindowField.setEnabled(!autoMode);
		iHighWindowField.setEnabled(!autoMode);

		return panel;
	}

	/** Builds the options panel (auto mode, show components, bins, save). */
	private JPanel buildOptionsPanel() {
		JPanel panel = new JPanel(new GridBagLayout());
		panel.setBorder(BorderFactory.createTitledBorder("Options"));
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.insets = new Insets(2, 4, 2, 4);
		gbc.anchor = GridBagConstraints.WEST;

		gbc.gridx = 0; gbc.gridy = 0; gbc.gridwidth = 2;
		autoModeCheckbox = new JCheckBox("Auto detect I-Low / I-High (double-Gaussian fit)", autoMode);
		autoModeCheckbox.addActionListener(e -> onAutoModeChanged());
		panel.add(autoModeCheckbox, gbc);

		gbc.gridy = 1;
		showComponentsCheckbox = new JCheckBox("Show individual Gaussian components on plot", showComponents);
		showComponentsCheckbox.setEnabled(autoMode);
		panel.add(showComponentsCheckbox, gbc);

		gbc.gridy = 2; gbc.gridwidth = 1;
		gbc.gridx = 0;
		panel.add(new JLabel("Number of bins:"), gbc);
		gbc.gridx = 1;
		nbinsField = new JTextField(String.valueOf(nbins), 6);
		panel.add(nbinsField, gbc);

		gbc.gridy = 3; gbc.gridx = 0; gbc.gridwidth = 2;
		saveFigsCheckbox = new JCheckBox("Save plot as titled figure", saveFigs);
		panel.add(saveFigsCheckbox, gbc);

		return panel;
	}

	/** Restarts the debounce timer for preview updates. */
	private void schedulePreviewUpdate() {
		if (updateTimer != null) {
			updateTimer.restart();
		}
	}

	/** Recomputes the RGB preview overlay from cached data and current slider values. */
	private void updatePreview() {
		if (previewImp == null || previewImp.getWindow() == null) return;

		// Gradient cutoff from pre-sorted array
		double gradThresh = gradientSliderField.getValue();
		int cutoffIdx = Math.min(
			(int) (gradThresh * sortedInteriorGradients.length),
			sortedInteriorGradients.length - 1);
		float gradCutoff = sortedInteriorGradients[cutoffIdx];

		boolean showIntensityOverlay = !autoModeCheckbox.isSelected() || autoFitPopulated;
		boolean showGradient = showGradientCheckbox.isSelected();
		float iLowVal = (float) iLowSliderField.getValue();
		float iHighVal = (float) iHighSliderField.getValue();

		float iLowDelta = 50, iHighDelta = 50;
		try { iLowDelta = Integer.parseInt(iLowWindowField.getText().trim()); } catch (NumberFormatException ignored) {}
		try { iHighDelta = Integer.parseInt(iHighWindowField.getText().trim()); } catch (NumberFormatException ignored) {}

		float range = maxIntensity - minIntensity;
		if (range == 0) range = 1;

		// Interior pixels only (borders pre-filled in run())
		for (int y = 1; y < cropHeight - 1; y++) {
			for (int x = 1; x < cropWidth - 1; x++) {
				int idx = y * cropWidth + x;
				float val = cachedSmoothedPixels[idx];
				int gray = OverlayTint.toGray(val, minIntensity, range);

				if (showIntensityOverlay && val >= iLowVal - iLowDelta && val <= iLowVal + iLowDelta) {
					previewRgb[idx] = OverlayTint.tint(gray, Col.THR_MIN);
				} else if (showIntensityOverlay && val >= iHighVal - iHighDelta && val <= iHighVal + iHighDelta) {
					previewRgb[idx] = OverlayTint.tint(gray, Col.THR_MAX);
				} else if (showGradient && cachedGradMagnitudes[idx] >= gradCutoff) {
					previewRgb[idx] = OverlayTint.tint(gray, Col.GRADIENT);
				} else {
					previewRgb[idx] = OverlayTint.grayscale(gray);
				}
			}
		}

		previewImp.updateAndDraw();
	}

	/** Toggles I-Low/I-High slider and component checkbox enabled states. */
	private void onAutoModeChanged() {
		boolean auto = autoModeCheckbox.isSelected();
		iLowSliderField.setEnabled(!auto);
		iHighSliderField.setEnabled(!auto);
		iLowWindowField.setEnabled(!auto);
		iHighWindowField.setEnabled(!auto);
		showComponentsCheckbox.setEnabled(auto);
		schedulePreviewUpdate();
	}

	/** Launches the Noise Statistics plugin in a background thread, then reloads I0. */
	private void onRunSNR() {
		new SwingWorker<Void, Void>() {
			@Override
			protected Void doInBackground() {
				IJ.run("Noise Statistics Analysis (Single Image)");
				return null;
			}

			@Override
			protected void done() {
				getAllPersistedParams();
				if (ranSNR) {
					i0Field.setText(String.format("%.2f", i0));
					IJ.log("I0 loaded from Noise Statistics: " + String.format("%.2f", i0));
				} else {
					IJ.error("Contrast", "SNR analysis did not complete. Please enter I0 manually.");
				}
			}
		}.execute();
	}

	/**
	 * Reads all current parameter values from UI controls into instance fields.
	 *
	 * @param strict if {@code true}, shows error dialogs and returns {@code false}
	 *               on parse failures; if {@code false}, silently ignores them
	 * @return {@code true} if all values were read successfully
	 */
	private boolean readUIValues(boolean strict) {
		try {
			i0 = Double.parseDouble(i0Field.getText().trim());
		} catch (NumberFormatException ex) {
			if (strict) { IJ.error("Contrast", "I0 must be a valid number."); return false; }
		}
		gradientThreshold = gradientSliderField.getValue();
		autoMode = autoModeCheckbox.isSelected();
		showComponents = showComponentsCheckbox.isSelected();
		saveFigs = saveFigsCheckbox.isSelected();
		showGradientPreview = showGradientCheckbox.isSelected();
		try { iLowWindow = Integer.parseInt(iLowWindowField.getText().trim()); } catch (NumberFormatException ignored) {}
		try { iHighWindow = Integer.parseInt(iHighWindowField.getText().trim()); } catch (NumberFormatException ignored) {}
		try {
			nbins = Integer.parseInt(nbinsField.getText().trim());
		} catch (NumberFormatException ex) {
			if (strict) { IJ.error("Contrast", "Number of bins must be a valid integer."); return false; }
		}
		if (!autoMode) {
			iLow = iLowSliderField.getValue();
			iHigh = iHighSliderField.getValue();
		}
		return true;
	}

	/** Validates parameters, runs contrast analysis in a background thread, and shows results. */
	private void onApply() {
		if (!readUIValues(true)) return;

		// Validate
		if (gradientThreshold < 0 || gradientThreshold > 1) {
			IJ.error("Contrast", "Gradient threshold must be between 0 and 1.");
			return;
		}
		if (!autoMode && iLow >= iHigh) {
			IJ.error("Contrast", "I-Low must be less than I-High.");
			return;
		}
		if (nbins < 10 || nbins > 10000) {
			IJ.error("Contrast", "Number of bins must be between 10 and 10000.");
			return;
		}

		applyButton.setEnabled(false);
		applyButton.setText("Computing...");

		new SwingWorker<ContrastData, Void>() {
			@Override
			protected ContrastData doInBackground() {
				if (autoMode) {
					return ContrastAnalyzer.computeContrastAuto(croppedIp, i0, gradientThreshold, nbins);
				} else {
					return ContrastAnalyzer.computeContrastManual(croppedIp, i0, gradientThreshold, iLow, iHigh, nbins);
				}
			}

			@Override
			protected void done() {
				try {
					ContrastData result = get();

					// If auto mode, update sliders to fitted values
					if (autoMode) {
						iLow = result.getILow();
						iHigh = result.getIHigh();
						autoFitPopulated = true;
						iLowSliderField.setValue(iLow);
						iHighSliderField.setValue(iHigh);
					}

					logResults(result);
					Plot pdfPlot = plotContrastPDF(result, imp.getTitle());

        if (saveFigs) {
            FileInfo fi = imp.getOriginalFileInfo();
            String dir = (fi != null && fi.directory != null) ? fi.directory : "";
            String baseName = ImageResolver.getBaseName(imp);

            FigBuilder.createAndSave(pdfPlot, "Contrast PDF",
                                     dir + baseName + "_contrast_pdf.png");
        }

					pdfPlot.show();
					setAllPersistedParams(true);
					updatePreview();
					IJ.showStatus("Contrast analysis complete.");

				} catch (Exception ex) {
					IJ.error("Contrast", "Analysis failed: " + ex.getMessage());
				} finally {
					applyButton.setEnabled(true);
					applyButton.setText("Apply");
				}
			}
		}.execute();
	}

	/** Persists parameters, closes preview window, and disposes the frame. */
	private void onClose() {
		readUIValues(false);
		setAllPersistedParams(false);

		if (previewImp != null) {
			previewImp.close();
			previewImp = null;
		}
		if (updateTimer != null) {
			updateTimer.stop();
		}
		if (frame != null) {
			frame.dispose();
			frame = null;
		}
	}

	/**
	 * Creates a plot showing the histogram (PDF) of the gradient-filtered
	 * smoothed subset, with vertical threshold markers and contrast annotation.
	 */
	private Plot plotContrastPDF(ContrastData result, String imageTitle) {
		ThresholdData td = result.getThresholdData();
		double[] pdf = td.getPDF();
		int nBins = pdf.length;
		double iLowR = result.getILow();
		double iHighR = result.getIHigh();
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

		Plot plot = new Plot("Contrast PDF: " + imageTitle, "Intensity", "Probability Density");
		StringBuilder legend = new StringBuilder();

        // PDF curve
        plot.setColor(Col.PDF);
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
            plot.setColor(Col.GAUSS_FIT);
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
                plot.setColor(Col.GAUSS_1);
                plot.addPoints(binCenters, g1, Plot.LINE);
                legend.append(String.format("Gaussian 1 (\u03BC=%.2f)\n", mu1));
                plot.setColor(Col.GAUSS_2);
                plot.addPoints(binCenters, g2, Plot.LINE);
                legend.append(String.format("Gaussian 2 (\u03BC=%.2f)\n", mu2));
            }
        }

        // Vertical threshold lines
		plot.setLineWidth(1);
		plot.setColor(Col.THR_MIN);
		plot.drawDottedLine(iLowR, 0, iLowR, maxPdf, 1);
		plot.addPoints(new double[] {iLowR, iLowR}, new double[] {0, maxPdf}, Plot.DOT);
		legend.append(String.format("I_low = %.2f\n", iLowR));
		plot.setColor(Col.THR_MAX);
		plot.drawDottedLine(iHighR, 0, iHighR, maxPdf, 1);
		plot.addPoints(new double[] {iHighR, iHighR}, new double[] {0, maxPdf}, Plot.DOT);
		legend.append(String.format("I_high = %.2f\n", iHighR));

		// Annotations as legend entries (invisible dummy series)
		plot.setColor(Col.KEY);
		plot.addPoints(new double[] {Double.NaN}, new double[] {Double.NaN}, Plot.DOT);
		legend.append(String.format("I0 = %.2f\n", i0));
		plot.addPoints(new double[] {Double.NaN}, new double[] {Double.NaN}, Plot.DOT);
		legend.append(String.format("Contrast = %.3f\n", contrast));

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
		saveFigs = ParamPersister.get(imp, "C_saveFigs", false);
		showGradientPreview = ParamPersister.get(imp, "C_showGradientPreview", true);
		iLowWindow = ParamPersister.get(imp, "C_iLowWindow", 50);
		iHighWindow = ParamPersister.get(imp, "C_iHighWindow", 50);
	}

	private void setAllPersistedParams(boolean showLog) {
		ParamPersister.set(imp, "C_ranSNR", ranSNR);
		ParamPersister.set(imp, "C_i0", i0);
        ParamPersister.set(imp, "C_gradientThreshold", gradientThreshold);
        ParamPersister.set(imp, "C_autoMode", autoMode);
		ParamPersister.set(imp, "C_iLow", iLow);
		ParamPersister.set(imp, "C_iHigh", iHigh);
		ParamPersister.set(imp, "C_nbins", nbins);
		ParamPersister.set(imp, "C_saveFigs", saveFigs);
		ParamPersister.set(imp, "C_showGradientPreview", showGradientPreview);
		ParamPersister.set(imp, "C_iLowWindow", iLowWindow);
		ParamPersister.set(imp, "C_iHighWindow", iHighWindow);
		if (showLog) {
			logParams();
		}
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