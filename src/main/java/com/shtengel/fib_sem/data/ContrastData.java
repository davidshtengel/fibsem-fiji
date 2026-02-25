package com.shtengel.fib_sem.data;

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
    private final double thrMinContrast;
    private final double thrMaxContrast;
	private final ThresholdData thresholdData;

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
     * @param thrMinContrast    lower CDF threshold used for peak identification
     * @param thrMaxContrast    upper CDF threshold used for peak identification
     * @param thresholdData     the underlying threshold analysis (histogram, PDF, CDF)
     */
	public ContrastData(double iLow, double iHigh, double i0,
                        double gradientThreshold,
                        int totalPixels, int subsetPixels,
                        double[] smoothedSubset,
                        double thrMinContrast, double thrMaxContrast,
					ThresholdData thresholdData) {
        this.iLow = iLow;
        this.iHigh = iHigh;
        this.iMean = (iHigh + iLow) / 2.0;
        this.i0 = i0;
        this.contrast = (iHigh - iLow) / (this.iMean - i0);
        this.gradientThreshold = gradientThreshold;
        this.totalPixels = totalPixels;
        this.subsetPixels = subsetPixels;
        this.smoothedSubset = smoothedSubset;
        this.thrMinContrast = thrMinContrast;
        this.thrMaxContrast = thrMaxContrast;
		this.thresholdData = thresholdData;
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

    public double getThrMinContrast() { return thrMinContrast; }
    public double getThrMaxContrast() { return thrMaxContrast; }

	/** Returns the underlying threshold analysis (histogram, PDF, CDF). */
    public ThresholdData getThresholdData() { return thresholdData; }
}
