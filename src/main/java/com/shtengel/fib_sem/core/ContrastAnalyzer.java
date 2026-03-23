package com.shtengel.fib_sem.core;

import java.util.ArrayList;

import com.shtengel.fib_sem.data.*;
import com.shtengel.fib_sem.util.DoubleGaussianFitter;

import ij.process.*;

/**
 * Computes the contrast metric for bimodal FIB-SEM intensity distributions.
 *
 * @see ContrastData
 */
public class ContrastAnalyzer {

	private ContrastAnalyzer() {}

	/**
	 * Computes contrast with automatically detected bimodal peaks.
	 *
	 * <p>Applies gradient filtering to retain low-texture pixels, smooths,
	 * then fits a double-Gaussian to the resulting PDF to derive I_low and I_high
	 * from the two mode centres.</p>
	 *
	 * @param ip                source image processor (ROI is respected if set)
	 * @param darkCount         dark count I0 (intensity at zero variance)
	 * @param gradientThreshold fraction of pixels to retain by gradient magnitude (0–1)
	 * @param nbins             number of histogram bins
	 * @return contrast analysis results (including the fitted double-Gaussian)
	 * @see ContrastData
	 */
	public static ContrastData computeContrastAuto(ImageProcessor ip,
		double darkCount,
		double gradientThreshold,
		int nbins) {

		double[][] subset = extractSubset(ip, gradientThreshold);
		double[] values = subset[0];
		int totalPixels = (int) subset[1][0];
		int subsetPixels = (int) subset[2][0];

		// Build histogram, PDF, and CDF
		ThresholdData thresholds = ThresholdAnalyzer.computeThresholds(
			values, 0.001, 0.001, nbins
		);
		
		// Compute bin centres for the fitter
		double[] pdf = thresholds.getPDF();
		double minInt = thresholds.getMinIntensity();
		double maxInt = thresholds.getMaxIntensity();
		int nBins = pdf.length;
		double binWidth = (maxInt - minInt) / nBins;
		double[] binCenters = new double[nBins];
		for (int i = 0; i < nBins; i++) {
			binCenters[i] = minInt + (i + 0.5) * binWidth;
		}

		// Fit double-Gaussian → I_low = μ1, I_high = μ2
		DoubleGaussianFitter.FitResult fit = DoubleGaussianFitter.fit(binCenters, pdf);

		double iLow  = fit.getMu1();
		double iHigh = fit.getMu2();

		return new ContrastData(iLow, iHigh, darkCount,
			gradientThreshold,
			totalPixels, subsetPixels,
			values,
			true,   // autoMode
			nbins,
			thresholds,
			fit
		);
	}

	/**
	 * Computes contrast using user-supplied I_low and I_high values.
	 *
	 * <p>Applies gradient filtering and smoothing as usual, but uses the provided
	 * intensity values directly rather than fitting peaks.</p>
	 *
	 * @param ip                source image processor (ROI is respected if set)
	 * @param darkCount         dark count I0 (intensity at zero variance)
	 * @param gradientThreshold fraction of pixels to retain by gradient magnitude (0–1)
	 * @param iLow              user-specified lower intensity peak value
	 * @param iHigh             user-specified upper intensity peak value
	 * @param nbins             number of histogram bins (used for visualization)
	 * @return contrast analysis results
	 * @see ContrastData
	 */
	public static ContrastData computeContrastManual(ImageProcessor ip,
		double darkCount,
		double gradientThreshold,
		double iLow,
		double iHigh,
		int nbins) {

		double[][] subset = extractSubset(ip, gradientThreshold);
		double[] values = subset[0];
		int totalPixels = (int) subset[1][0];
		int subsetPixels = (int) subset[2][0];

		ThresholdData thresholds = ThresholdAnalyzer.computeThresholds(
			values, 0.001, 0.001, nbins
		);

		return new ContrastData(iLow, iHigh, darkCount,
			gradientThreshold,
			totalPixels, subsetPixels,
			values,
			false,  // manual mode
			nbins,
			thresholds,
			null    // no fit in manual mode
		);
	}

	/**
	 * Computes gradient magnitudes for the given image processor.
	 *
	 * <p>Crops to ROI if set, converts to float, then computes gradient
	 * magnitude via central finite differences (no smoothing, no normalization).</p>
	 *
	 * @param ip source image processor (ROI is respected if set)
	 * @return flattened gradient magnitude array (length = croppedWidth * croppedHeight)
	 */
	public static float[] computeGradientMagnitudes(ImageProcessor ip) {
		ImageProcessor cropped = (ip.getRoi() != null) ? ip.crop() : ip;
		FloatProcessor fp = cropped.convertToFloatProcessor();
		float[][] components = GradientMapAnalyzer.computeGradientComponents(
			fp.duplicate().convertToFloatProcessor(), false);
		return GradientMapAnalyzer.computeGradientMagnitude(
			components[0], components[1], fp, false);
	}

	/**
	 * Computes smoothed pixel values for the given image processor.
	 *
	 * <p>Crops to ROI if set, converts to float, duplicates, applies
	 * default smoothing, and returns the smoothed pixel array.</p>
	 *
	 * @param ip source image processor (ROI is respected if set)
	 * @return flattened smoothed pixel array (length = croppedWidth * croppedHeight)
	 */
	public static float[] computeSmoothedPixels(ImageProcessor ip) {
		ImageProcessor cropped = (ip.getRoi() != null) ? ip.crop() : ip;
		FloatProcessor fp = cropped.convertToFloatProcessor();
		FloatProcessor smoothedFp = fp.duplicate().convertToFloatProcessor();
		GradientMapAnalyzer.applyDefaultSmoothing(smoothedFp);
		return (float[]) smoothedFp.getPixels();
	}

	/**
	 * Extracts the gradient-filtered, smoothed pixel subset from the image.
	 */
	private static double[][] extractSubset(ImageProcessor ip, double gradientThreshold) {
		ImageProcessor cropped = (ip.getRoi() != null) ? ip.crop() : ip;
		int width = cropped.getWidth();
		int height = cropped.getHeight();
		int totalPixels = width * height;

		float[] gradMag = computeGradientMagnitudes(ip);
		float gradientCutoff = computeGradientCutoff(gradMag, width, height, gradientThreshold);
		float[] smoothed = computeSmoothedPixels(ip);

		// Extract subset where gradient < cutoff (excluding 1-pixel border)
		ArrayList<Double> subsetList = new ArrayList<>();
		for (int y = 1; y < height - 1; y++) {
			for (int x = 1; x < width - 1; x++) {
				int idx = y * width + x;
				if (gradMag[idx] < gradientCutoff) {
					subsetList.add((double) smoothed[idx]);
				}
			}
		}
		double[] values = toArray(subsetList);
		return new double[][] {values, {totalPixels}, {subsetList.size()}};
	}

	/**
     * Determines the gradient cutoff value such that a given fraction of interior
     * pixels have gradient magnitude below this cutoff.
     *
     * @param gradients          array of gradient magnitudes (full image, flattened)
     * @param width              image width
     * @param height             image height
     * @param gradientThreshold  fraction of pixels to retain (0 to 1)
     * @return the gradient cutoff value
     */
    public static float computeGradientCutoff(float[] gradients, int width, int height, double gradientThreshold) {
        // Collect interior gradient values (skip 1-pixel border)
        int interiorCount = (width - 2) * (height - 2);
        float[] interior = new float[interiorCount];
        int k = 0;
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                interior[k++] = gradients[y * width + x];
            }
        }
        java.util.Arrays.sort(interior);

        int cutoffIdx = Math.min((int) (gradientThreshold * interiorCount), interiorCount - 1);
        return interior[cutoffIdx];
    }

	private static double[] toArray(ArrayList<Double> list) {
		double[] arr = new double[list.size()];

		for (int i = 0; i < list.size(); i++) {
			arr[i] = list.get(i);
		}

		return arr;
	}
}
