package com.shtengel.fib_sem.core;

import com.shtengel.fib_sem.data.ThresholdData;

import ij.process.ImageProcessor;

/**
 * Core utility for computing intensity thresholds from image histograms.
 *
 * <p>
 * Thresholds are derived from histogram-based quantiles using the cumulative
 * distribution function (CDF). This implementation is lightweight and
 * independent of ImageJ's built-in auto-thresholding methods.
 * </p>
 */
public class ThresholdAnalyzer {
	
	/**
     * Computes lower and upper intensity thresholds using CDF cutoffs.
     *
     * <p>
     * If an ROI is set on the input {@link ImageProcessor}, only pixels within
     * the ROI are analyzed.
     * </p>
     *
     * @param ip      source image processor
     * @param thrMin  fraction of lowest-intensity pixels to discard [0, 1]
     * @param thrMax  fraction of highest-intensity pixels to discard [0, 1]
     * @param nbins   number of histogram bins
     * @return container with absolute intensity range, thresholds, PDF, and CDF
     */
    public static ThresholdData computeThresholds(ImageProcessor ip, double thrMin, double thrMax, int nbins) {
    	// Crop to ROI if present; output dimensions match ROI bounds
    	ImageProcessor ipToProcess = (ip.getRoi() != null) ? ip.crop() : ip;
        
        // Convert to float
        float[] pixels = (float[]) ipToProcess.convertToFloatProcessor().getPixels();
        int nPixels = pixels.length;
        
    	// Compute absolute min and max intensity
        float minInt = Float.MAX_VALUE;
        float maxInt = -Float.MAX_VALUE;
        for (float v : pixels) {
            if (v < minInt) minInt = v;
            if (v > maxInt) maxInt = v;
        }

        // Early exit for uniform images
        if (minInt == maxInt) {
            return new ThresholdData(minInt, maxInt);
        }

        // Compute histogram of pixel intensity
        double binWidth = (maxInt - minInt) / nbins;
        double[] hist = new double[nbins];
        for (float v : pixels) {
            int idx = (int) (((v - minInt) / (maxInt - minInt)) * nbins);
            if (idx >= nbins) idx = nbins - 1;
            if (idx < 0) idx = 0;
            hist[idx]++;
        }
        
        // Compute probability density function (PDF) and cumulative distribution function (CDF)
        // PDF normalizes histogram counts by total pixels, and CDF determines probability a pixel is <= a given intensity
        double[] pdf = new double[nbins];
        double[] cdf = new double[nbins];
        double cumsum = 0.0;
        for (int i = 0; i < nbins; i++) {
            pdf[i] = hist[i] / nPixels;
            cumsum += pdf[i];
            cdf[i] = cumsum;
        }

        // Find the threshold bin indices, the lowest bin where CDF >= thrMin and the highest bin where CDF >= 1 - thrMax
        int idxDataMin = 0;
        int idxDataMax = nbins - 1;

        // Lower threshold
        for (int i = 0; i < nbins; i++) {
            if (cdf[i] >= thrMin) {
                idxDataMin = i;
                break;
            }
        }
        // Upper threshold
        for (int i = 0; i < nbins; i++) {
            if (cdf[i] >= (1.0 - thrMax)) { 
                idxDataMax = i; 
                break; 
            }
        }
        // Convert bin indices to intensity values
        double dataMin = minInt + ((idxDataMin + 0.5) * binWidth);
        double dataMax = minInt + ((idxDataMax + 0.5) * binWidth);

        return new ThresholdData(minInt, maxInt, dataMin, dataMax, pdf, cdf);
    }
    
}
