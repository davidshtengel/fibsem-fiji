package com.shtengel.fib_sem.data;

/**
 * Container for variance–mean noise analysis results.
 *
 * <p>
 * This class stores both raw binned statistics (mean, variance) and
 * derived quantities such as fitted slopes, intercepts, peak intensity,
 * and signal-to-noise ratios.
 * </p>
 *
 * <p>
 * All arrays correspond to the same binning used during analysis.
 * </p>
 */
public class NoiseStatisticsData {
	/** Mean intensity per analysis bin. */
    private final double[] meanValues;

    /** Variance of residual (noise) per analysis bin. */
    private final double[] varianceValues;

    /** Zero-intercept intensity (i_0) from free linear fit. */
    private final double i0;

    /** Signal-to-noise ratio from free fit. */
    private final double snr;

    /** Signal-to-noise ratio using fixed dark-count offset. */
    private final double snr1;

    /** Slope of variance–mean relationship from free fit. */
    private final double slope;

    /** Slope of variance–mean relationship with dark-count correction. */
    private final double slopeHeader;

    /** Peak intensity estimated from display-range histogram. */
    private final double iPeak;

    /** Variance evaluated at peak intensity. */
    private final double varPeak;

    /** Intensity range used for variance–mean analysis. */
    private final double[] rangeAnalysis;

    /** Intensity range used for display and histogram visualization. */
    private final double[] rangeDisplay;
    
    public NoiseStatisticsData(double[] meanValues, 
							double[] varianceValues, 
							double i0, 
							double snr, 
							double snr1, 
							double slope, 
							double slopeHeader, 
							double iPeak, 
							double varPeak, 
							double[] rangeAnalysis, 
							double[] rangeDisplay) {
        this.meanValues = meanValues;
        this.varianceValues = varianceValues;
        this.i0 = i0;
        this.snr = snr;
        this.snr1 = snr1;
        this.slope = slope;
        this.slopeHeader = slopeHeader;
        this.iPeak = iPeak;
        this.varPeak = varPeak;
        this.rangeAnalysis = rangeAnalysis;
        this.rangeDisplay = rangeDisplay;
    }
    
    public double[] getMeanValues() { return meanValues; }
    public double[] getVarianceValues() { return varianceValues; }
    public double getI0() { return i0; }
    public double getSNR() { return snr; }
    public double getSNR1() { return snr1; }
    public double getSlope() { return slope; }
    public double getSlopeHeader() { return slopeHeader; }
    public double getIPeak() { return iPeak; }
    public double getVarPeak() { return varPeak; }
    public double[] getRangeAnalysis() { return rangeAnalysis; }
    public double[] getRangeDisplay() { return rangeDisplay; }
    
    @Override
    public String toString() {
        return String.format(
            "Noise Statistics[SNR=%.2f, SNR1=%.2f, I0=%.2f, iPeak=%.2f]",
            snr, snr1, i0, iPeak
        );
    }
}