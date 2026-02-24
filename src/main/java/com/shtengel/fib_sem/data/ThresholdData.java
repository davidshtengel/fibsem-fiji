package com.shtengel.fib_sem.data;

/**
 * Container for histogram-based threshold analysis results.
 *
 * <p>
 * Holds absolute intensity bounds, computed quantile thresholds,
 * and the underlying PDF/CDF used to derive them.
 * </p>
 */
public class ThresholdData {
    private double minIntensity;
    private double maxIntensity;
    private double minThreshold;
    private double maxThreshold;
    private double[] pdf;
    private double[] cdf;
    
    /**
     * Constructor for degenerate cases (e.g. uniform images).
     */
    public ThresholdData(double minInt, double maxInt) {
        this.minIntensity = minInt;
        this.maxIntensity = maxInt;
    }
    
    /**
     * Full result constructor.
     */
    public ThresholdData(double minInt, double maxInt, double minThr, double maxThr, double[] pdf, double[] cdf) {
        this.minIntensity = minInt;
        this.maxIntensity = maxInt;
        this.minThreshold = minThr;
        this.maxThreshold = maxThr;
        this.pdf = pdf;
        this.cdf = cdf;
    }

    public double getMinIntensity() { return minIntensity; }
    public double getMaxIntensity() { return maxIntensity; }
    public double getMinThreshold() { return minThreshold; }
    public double getMaxThreshold() { return maxThreshold; }
    public double[] getPDF() { return pdf; }
    public double[] getCDF() { return cdf; }


    @Override
    public String toString() {
        return String.format(" Computed Thresholds[min=%.2f, max=%.2f]", minThreshold, maxThreshold);
    }
}
