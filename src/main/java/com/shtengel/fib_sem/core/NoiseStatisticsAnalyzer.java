package com.shtengel.fib_sem.core;

import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Arrays;

import com.shtengel.fib_sem.data.NoiseStatisticsData;
import com.shtengel.fib_sem.data.ThresholdData;
import com.shtengel.fib_sem.util.ImageResolver;

/**
 * Core implementation of variance–mean noise analysis for FIB-SEM images.
 *
 * <p>
 * The algorithm estimates noise by subtracting a locally smoothed version
 * of the image from the original, treating the residual as noise. Pixels
 * are binned by smoothed intensity, and noise variance is computed per bin.
 * </p>
 *
 * <p>
 * Linear regression of variance vs. mean is used to estimate system gain
 * and signal-to-noise ratios. Optional dark-count correction allows
 * enforcing a known zero-signal offset.
 * </p>
 *
 * <p>
 * <b>ROI Support:</b> When an ROI is active, only pixels within the ROI are
 * analyzed.
 * </p>
 */
public class NoiseStatisticsAnalyzer {
	
	/**
     * Computes noise statistics from an image by analyzing variance vs mean intensity.
     * If an ROI is set on the ImageProcessor, only pixels within the ROI are analyzed.
     *
     * @param ip                   source image processor
     * @param darkCount            expected zero-signal intensity offset
     * @param kernel               smoothing kernel (null → default 3×3)
     * @param nbinsDisp            histogram bins for display-range estimation
     * @param thresholdsDisp       [lower, upper] CDF cutoffs for display range
     * @param nbinsAnalysis        bins for variance–mean analysis
     * @param thresholdsAnalysis   [lower, upper] CDF cutoffs for analysis range
     * @param gradientThreshold    upper CDF threshold for gradient filtering (0.0-1.0)
     * @param filterArray          optional pixel mask (null → use all pixels)
     * @return populated {@link NoiseStatisticsData}
     */
	public static NoiseStatisticsData computeNoiseStatistics(
			ImageProcessor ip, 
			double darkCount, 
			float[] kernel, 
			int nbinsDisp, 
			double[] thresholdsDisp, 
			int nbinsAnalysis, 
			double[] thresholdsAnalysis, 
			double gradientThreshold,
			boolean[] filterArray) {
		
    	// Crop to ROI if present; output dimensions match ROI bounds
    	ImageProcessor ipToProcess = ImageResolver.cropToRoiIfPresent(ip);
   
	    // Convert to float
	    FloatProcessor fp = ipToProcess.convertToFloatProcessor();
		int width = fp.getWidth();
		int height = fp.getHeight();
				
		// Use default kernel if not provided
		if (kernel == null) {
			kernel = getDefaultSmoothingKernel();
		}
		
		// Smooth the image (convolution respects ROI internally)
		FloatProcessor fpSmoothed = smoothImage(fp, kernel);
        
        // Compute gradient magnitudes for filtering
        float[] gradientMagnitudes = computeGradientMagnitudes(fp);
        
        // Collect pixels w/in ROI bounds (exclude a 1-pixel border due to convolution)
        ArrayList<Float> pixelsSmoothedList = new ArrayList<>();
        ArrayList<Float> pixelsDiffList = new ArrayList<>();
        ArrayList<Float> gradientList = new ArrayList<>();
        ArrayList<Boolean> maskList = new ArrayList<>();
        
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                float origVal = fp.getf(x, y);
                float smoothedVal = fpSmoothed.getf(x, y);
                float diffVal = origVal - smoothedVal;
                float gradMag = gradientMagnitudes[y * width + x];
                
                pixelsSmoothedList.add(smoothedVal);
                pixelsDiffList.add(diffVal);
                gradientList.add(gradMag);
                
                // Apply filter if provided (filterArray would need to be ROI-sized or mapped differently)
                if (filterArray != null && filterArray.length == width * height) {
                    maskList.add(filterArray[y * width + x]);
                } else {
                    maskList.add(true);
                }
            }
        }
        
        // Convert to arrays
        float[] pixelsSmoothed = toFloatArray(pixelsSmoothedList);
        float[] pixelsDiff = toFloatArray(pixelsDiffList);
        float[] gradients = toFloatArray(gradientList);
        boolean[] mask = toBooleanArray(maskList);
        
        // Apply gradient threshold filtering
        boolean[] gradientMask = applyGradientThreshold(gradients, gradientThreshold);
        
        // Combine masks
        for (int i = 0; i < mask.length; i++) {
            mask[i] = mask[i] && gradientMask[i];
        }
        
        // Apply mask to create filtered arrays
        float[] smoothedFilt = applyMask(pixelsSmoothed, mask);
        float[] diffFilt = applyMask(pixelsDiff, mask);
        
        // Check if there's enough data
        if (smoothedFilt.length < 10) {
            throw new IllegalArgumentException(
                "Insufficient data points in ROI for analysis (found " + smoothedFilt.length + ", need at least 10)"
            );
        }
        
        // Compute display and analysis ranges using ThresholdAnalyzer
        FloatProcessor fpSmoothedFilt = new FloatProcessor(smoothedFilt.length, 1, smoothedFilt);
        ThresholdData threshDisp = ThresholdAnalyzer.computeThresholds(
            fpSmoothedFilt, thresholdsDisp[0], thresholdsDisp[1], nbinsDisp
        );
        ThresholdData threshAnalysis = ThresholdAnalyzer.computeThresholds(
            fpSmoothedFilt, thresholdsAnalysis[0], thresholdsAnalysis[1], nbinsAnalysis
        );
        
        double[] rangeDisplay = new double[] {
            threshDisp.getMinThreshold(),
            threshDisp.getMaxThreshold()
        };
        double[] rangeAnalysis = new double[] {
            threshAnalysis.getMinThreshold(),
            threshAnalysis.getMaxThreshold()
        };
        
        // Build histogram for analysis range
        double[] bins = linspace(rangeAnalysis[0], rangeAnalysis[1], nbinsAnalysis);
        
        // Digitize smoothed data into bins
        int[] binIndices = new int[smoothedFilt.length];
        for (int i = 0; i < smoothedFilt.length; i++) {
            binIndices[i] = findBin(smoothedFilt[i], bins);
        }
        
        // For each bin, compute mean intensity and variance
        ArrayList<Double> meanVals = new ArrayList<>();
        ArrayList<Double> varVals = new ArrayList<>();
        
        for (int binIdx = 0; binIdx < nbinsAnalysis - 1; binIdx++) {
            ArrayList<Float> smoothedInBin = new ArrayList<>();
            ArrayList<Float> diffInBin = new ArrayList<>();
            
            for (int i = 0; i < smoothedFilt.length; i++) {
                if (binIndices[i] == binIdx) {
                    smoothedInBin.add(smoothedFilt[i]);
                    diffInBin.add(diffFilt[i]);
                }
            }
            
            if (smoothedInBin.size() >= 5) {
                double meanSmoothed = mean(smoothedInBin);
                double varDiff = variance(diffInBin);
                
                if (!Double.isNaN(meanSmoothed) && !Double.isNaN(varDiff)) {
                    meanVals.add(meanSmoothed);
                    varVals.add(varDiff);
                }
            }
        }
        
        double[] meanArray = toDoubleArray(meanVals);
        double[] varArray = toDoubleArray(varVals);
        
        // Find peak intensity from histogram
        double iPeak = findPeakIntensity(smoothedFilt, rangeDisplay, nbinsDisp);
        
        // Perform linear fit: var = slope * mean + intercept
        double[] fitParams = linearFit(meanArray, varArray);
        double slope = fitParams[0];
        double intercept = fitParams[1];
        double i0 = (Math.abs(slope) > 1e-10) ? -intercept / slope : 0.0;
        
        // Variance at peak
        double varPeak = slope * iPeak + intercept;
        
        // Fit with dark count (forced intercept)
        double slopeHeader = meanVarianceRatio(meanArray, varArray, darkCount);
        
        // Compute SNR values
        double snr = computeSNR(smoothedFilt, diffFilt, i0);
        double snr1 = computeSNR(smoothedFilt, diffFilt, darkCount);
        
        return new NoiseStatisticsData(
            meanArray,
            varArray,
            i0,
            snr,
            snr1,
            slope,
            slopeHeader,
            iPeak,
            varPeak,
            rangeAnalysis,
            rangeDisplay
        );
    }
	
	/**
	 * Computes local gradient magnitude for each pixel using central differences.
	 * Border pixels are set to zero.
	 */
	private static float[] computeGradientMagnitudes(FloatProcessor fp) {
	    int width = fp.getWidth();
	    int height = fp.getHeight();
	    float[] gradients = new float[width * height];

	    for (int y = 1; y < height - 1; y++) {
	        for (int x = 1; x < width - 1; x++) {
	            float gx = fp.getf(x + 1, y) - fp.getf(x - 1, y);
	            float gy = fp.getf(x, y + 1) - fp.getf(x, y - 1);
	            gradients[y * width + x] = (float) Math.sqrt(gx * gx + gy * gy);
	        }
	    }

	    return gradients;
	}
	
	/**
	 * Public wrapper for gradient magnitude computation.
	 *
	 * @param ip source image processor
	 * @return flattened array of gradient magnitudes
	 */
	public static float[] computeGradientMagnitudes(ImageProcessor ip) {
	    FloatProcessor fp = ip.convertToFloatProcessor();
	    return computeGradientMagnitudes(fp);
	}
	
	/**
	 * Returns a mask selecting pixels whose gradient magnitude lies below
	 * the specified CDF threshold.
	 *
	 * @param gradients         per-pixel gradient magnitudes
	 * @param gradientThreshold fraction (0–1) of lowest-gradient pixels to keep
	 */
	private static boolean[] applyGradientThreshold(
	        float[] gradients,
	        double gradientThreshold) {

	    int n = gradients.length;
	    boolean[] mask = new boolean[n];
	    double t = Math.max(0.0, Math.min(1.0, gradientThreshold));

	    if (t >= 1.0) { Arrays.fill(mask, true); return mask; }
	    if (t <= 0.0) { Arrays.fill(mask, false); return mask; }

	    float cutoff = computeGradientCutoff(gradients, gradientThreshold);

	    // Build mask: keep lowest-gradient pixels
	    for (int i = 0; i < n; i++) {
	        mask[i] = (gradients[i] <= cutoff);
	    }

	    return mask;
	}
	
	/**
	 * Returns elements of {@code data} where the corresponding {@code mask} entry is true.
	 *
	 * @param data  source array
	 * @param mask  boolean mask (same length as data)
	 * @return filtered array containing only unmasked elements
	 */
	public static float[] applyMask(float[] data, boolean[] mask) {
		ArrayList<Float> filtered = new ArrayList<>();
		for(int i = 0; i < data.length; i++) {
			if (mask[i]) {
				filtered.add(data[i]);
			}
		}
		
		return toFloatArray(filtered);
	}
	
	/**
	 * Computes the gradient threshold cutoff value based on CDF.
	 *
	 * @param gradients array of gradient magnitudes
	 * @param gradientThreshold CDF fraction (0-1)
	 * @return cutoff value for gradient filtering
	 */
	public static float computeGradientCutoff(float[] gradients, double gradientThreshold) {
	    int n = gradients.length;
	    
	    // Clamp threshold
	    double t = Math.max(0.0, Math.min(1.0, gradientThreshold));
	    
	    if (t >= 1.0) return Float.MAX_VALUE;
	    if (t <= 0.0) return Float.MIN_VALUE;
	    
	    // Copy & sort gradients
	    float[] sorted = gradients.clone();
	    Arrays.sort(sorted);
	    
	    // Index corresponding to CDF cutoff
	    int cutoffIndex = (int) Math.floor(t * (n - 1));
	    if (cutoffIndex < 0) cutoffIndex = 0;
	    if (cutoffIndex >= n) cutoffIndex = n - 1;
	    
	    return sorted[cutoffIndex];
	}
	
	/**
     * Default 3x3 smoothing kernel used to estimate local signal.
     */
    private static float[] getDefaultSmoothingKernel() {
        double st = 1.0 / Math.sqrt(2.0);
        float[] kernel = new float[] {
            (float)st, 1f, (float)st,
            1f, 1f, 1f,
            (float)st, 1f, (float)st
        };
        
        // Normalize
        float sum = 0;
        for (float v : kernel) sum += v;
        for (int i = 0; i < kernel.length; i++) {
            kernel[i] /= sum;
        }
        
        return kernel;
    }
    
    private static FloatProcessor smoothImage(FloatProcessor fp, float[] kernel) {
    	FloatProcessor fpSmoothed = (FloatProcessor) fp.duplicate();
    	Convolver conv = new Convolver();
    	conv.convolve(fpSmoothed, kernel, 3, 3);
    	return fpSmoothed;
    }
    
    /**
     * Computes smoothed version of image for visualization purposes.
     * Public wrapper around internal smoothing logic.
     *
     * @param ip source image processor
     * @return flattened array of smoothed pixel values
     */
    public static float[] computeSmoothedImage(ImageProcessor ip) {
        FloatProcessor fp = ip.convertToFloatProcessor();
        float[] kernel = getDefaultSmoothingKernel();
        
        return (float[]) smoothImage(fp, kernel).getPixels();
    }
    
    /**
     * Creates linearly spaced array
     */
    private static double[] linspace(double start, double end, int n) {
        double[] result = new double[n];
        double step = (end - start) / (n - 1);
        for (int i = 0; i < n; i++) {
            result[i] = start + i * step;
        }
        return result;
    }
    
    /**
     * Finds which bin a value belongs to
     */
    private static int findBin(float value, double[] bins) {
    	if (value < bins[0] || value >= bins[bins.length - 1]) return -1;
        
        for (int i = 0; i < bins.length - 1; i++) {
            if (value >= bins[i] && value < bins[i + 1]) {
                return i;
            }
        }
        return -1;
    }
    
    /**
     * Computes mean of a list
     */
    private static double mean(ArrayList<Float> values) {
        if (values.isEmpty()) return Double.NaN;
        double sum = 0;
        for (float v : values) sum += v;
        return sum / values.size();
    }
    
    /**
     * Computes variance of a list (uses N-1 denominator)
     */
    private static double variance(ArrayList<Float> values) {
        if (values.isEmpty()) return Double.NaN;
        double m = mean(values);
        double sumSq = 0;
        for (float v : values) {
            double diff = v - m;
            sumSq += diff * diff;
        }
        return sumSq / (values.size() - 1);
    }
    
    /**
     * Finds peak intensity from histogram
     */
    private static double findPeakIntensity(float[] data, double[] range, int nbins) {
        // Build histogram
        int[] hist = new int[nbins];
        double binWidth = (range[1] - range[0]) / nbins;
        
        for (float v : data) {
            if (v >= range[0] && v <= range[1]) {
                int binIdx = (int)((v - range[0]) / binWidth);
                if (binIdx >= nbins) binIdx = nbins - 1;
                if (binIdx < 0) binIdx = 0;
                hist[binIdx]++;
            }
        }
        
        // Find max bin
        int maxIdx = 0;
        int maxCount = hist[0];
        for (int i = 1; i < nbins; i++) {
            if (hist[i] > maxCount) {
                maxCount = hist[i];
                maxIdx = i;
            }
        }
        
        // Return bin center
        return range[0] + (maxIdx + 0.5) * binWidth;
    }
    
    /**
     * Ordinary least-squares linear fit: y = slope * x + intercept.
	 *
	 * @return {slope, intercept}, or {1.0, 0.0} if input is empty or mismatched
     */
    private static double[] linearFit(double[] x, double[] y) {
        if (x.length != y.length || x.length == 0) {
            return new double[] {1.0, 0.0};
        }
        
        int n = x.length;
        double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
        
        for (int i = 0; i < n; i++) {
            sumX += x[i];
            sumY += y[i];
            sumXY += x[i] * y[i];
            sumXX += x[i] * x[i];
        }
        
        double slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
        double intercept = (sumY - slope * sumX) / n;
        
        return new double[] {slope, intercept};
    }
    
    /**
     * Estimates the variance–mean slope assuming a fixed offset by
     * averaging per-bin ratios.
     */
    private static double meanVarianceRatio(double[] means, double[] vars, double offset) {
        double sumRatio = 0;
        int count = 0;
        
        for (int i = 0; i < means.length; i++) {
            double denominator = means[i] - offset;
            if (Math.abs(denominator) > 1e-6) {
                sumRatio += vars[i] / denominator;
                count++;
            }
        }
        
        return count > 0 ? sumRatio / count : 1.0;
    }
    
    /**
     * Computes signal-to-noise ratio as ⟨(s_d-B_d)²⟩ / ⟨n_d²⟩.
     *
     * <p>
     * Signal is taken as the smoothed intensity minus offset, and noise
     * as the residual between original and smoothed images.
     * </p>
     */
    private static double computeSNR(float[] signal, float[] noise, double offset) {
        double sumSigSq = 0;
        double sumNoiseSq = 0;
        
        for (int i = 0; i < signal.length; i++) {
            double s = signal[i] - offset;
            sumSigSq += s * s;
            sumNoiseSq += noise[i] * (double) noise[i];
        }
           
        return sumNoiseSq > 0 ? (sumSigSq / signal.length) / (sumNoiseSq / noise.length) : 0;
    }
    
    /**
     * Helper to convert ArrayList<Float> to float[]
     */
    private static float[] toFloatArray(ArrayList<Float> list) {
        float[] array = new float[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }
    
    /**
     * Helper to convert ArrayList<Double> to double[]
     */
    private static double[] toDoubleArray(ArrayList<Double> list) {
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }
    
    /**
     * Helper to convert ArrayList<Boolean> to boolean[]
     */
    private static boolean[] toBooleanArray(ArrayList<Boolean> list) {
        boolean[] array = new boolean[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }
}