package com.shtengel.fib_sem.data;

import com.shtengel.fib_sem.util.DoubleGaussianFitter;

/**
 * Data container for image contrast analysis results.
 *
 * <p>
 * Stores the outputs of a contrast calculation based on the bimodal intensity
 * distribution of a FIB-SEM image. The contrast metric quantifies separation
 * between high- and low-intensity material phases relative to the dark count.
 * </p>
 *
 * <p>
 * The contrast is defined as:
 * <pre>
 *   contrast = (I_high - I_low) / (I_mean - I0)
 * </pre>
 * where {@code I_mean = (I_high + I_low) / 2} and {@code I0} is the dark count
 * (intensity at zero signal).
 * </p>
 * 
 * <p>
 * In automatic mode, I_low and I_high are the fitted means (μ1, μ2) of a
 * double-Gaussian model; the full fit result is available via
 * {@link #getGaussianFit()}.
 * In manual mode, I_low and I_high are user-supplied values.
 * </p>
 *
 * @see com.shtengel.fib_sem.core.ContrastAnalyzer
 */
public class ContrastData {

	private final double iLow;
    private final double iHigh;
    private final double iMean;
    private final double i0;
    private final double contrast;
    private final double gradientThreshold;
    private final int totalPixels;
    private final int subsetPixels;
    private final double[] smoothedSubset;
	private final boolean autoMode;
	private final int nbins;
	private final ThresholdData thresholdData;
	private final DoubleGaussianFitter.FitResult gaussianFit;
	/**
     * Constructs a new {@code ContrastData} instance.
     *
     * @param iLow              lower intensity peak (CDF = thrMinContrast)
     * @param iHigh             upper intensity peak (CDF = 1 - thrMaxContrast)
     * @param i0                dark count (intensity at zero variance)
     * @param gradientThreshold gradient threshold used for pixel filtering
     * @param totalPixels       total number of pixels in the (ROI-cropped) image
     * @param subsetPixels      number of pixels retained after gradient filtering
     * @param smoothedSubset    the filtered subset of smoothed pixel values
	 * @param autoMode          {@code true} if peaks were detected via double-Gaussian fit
	 * @param nbins             number of histogram bins used
	 * @param thresholdData     the underlying threshold analysis (histogram, PDF, CDF)
	 * @param gaussianFit       the double-Gaussian fit result, or {@code null} for manual mode
	 */
	public ContrastData(double iLow, double iHigh, double i0,
						double gradientThreshold,
						int totalPixels, int subsetPixels,
						double[] smoothedSubset,
						boolean autoMode, int nbins,
						ThresholdData thresholdData,
						DoubleGaussianFitter.FitResult gaussianFit) {
		this.iLow = iLow;
		this.iHigh = iHigh;
		this.iMean = (iHigh + iLow) / 2.0;
		this.i0 = i0;
		this.contrast = (iHigh - iLow) / (this.iMean - i0);
		this.gradientThreshold = gradientThreshold;
		this.totalPixels = totalPixels;
		this.subsetPixels = subsetPixels;
		this.smoothedSubset = smoothedSubset;
		this.autoMode = autoMode;
		this.nbins = nbins;
		this.thresholdData = thresholdData;
		this.gaussianFit = gaussianFit;
	}

    public double getILow() { return iLow; }
    public double getIHigh() { return iHigh; }
    public double getIMean() { return iMean; }
    public double getI0() { return i0; }

    public double getContrast() { return contrast; }
    public double getGradientThreshold() { return gradientThreshold; }

	public int getTotalPixels() { return totalPixels; }
    public int getSubsetPixels() { return subsetPixels; }
    public double[] getSmoothedSubset() { return smoothedSubset; }

	/** Returns {@code true} if I_low and I_high were detected via double-Gaussian fit. */
	public boolean isAutoMode() { return autoMode; }

	public int getNbins() { return nbins; }
	
	public ThresholdData getThresholdData() { return thresholdData; }
	
	/**
	 * Returns the double-Gaussian fit result, or {@code null} if manual mode was used.
	 *
	 * <p>When non-null, the fit parameters are {@code [A1, μ1, σ1, A2, μ2, σ2]}
	 * with μ1 = I_low and μ2 = I_high.</p>
	 */
	public DoubleGaussianFitter.FitResult getGaussianFit() { return gaussianFit; }
}
