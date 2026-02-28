package com.shtengel.fib_sem.core;

import java.util.ArrayList;

import com.shtengel.fib_sem.data.*;
import ij.process.*;

/**
 * Computes the contrast metric for bimodal FIB-SEM intensity distributions.
 *
 * @see ContrastData
 */
public class ContrastAnalyzer {

	private ContrastAnalyzer() {}

	/**
	 * Computes contrast between high- and low-intensity phases of a FIB-SEM image.
	 *
	 * <p>Applies gradient filtering to retain low-texture pixels, smooths,
	 * then identifies intensity peaks via CDF thresholds to derive the contrast metric.</p>
	 *
	 * @param ip                source image processor (ROI is respected if set)
	 * @param darkCount         dark count I0 (intensity at zero variance)
	 * @param gradientThreshold fraction of pixels to retain by gradient magnitude (0–1)
	 * @param thrMinContrast    lower CDF quantile for low-intensity peak
	 * @param thrMaxContrast    upper CDF quantile for high-intensity peak
	 * @param nbins             number of histogram bins for threshold computation
	 * @return contrast analysis results
	 * @see ContrastData
	 */
	public static ContrastData computeContrast(ImageProcessor ip,
		double darkCount,
		double gradientThreshold,
		double thrMinContrast,
		double thrMaxContrast,
		int nbins) {

		// Work on the ROI-cropped region, if present
		ImageProcessor cropped = (ip.getRoi() != null) ? ip.crop() : ip;
		FloatProcessor fp = cropped.convertToFloatProcessor();
		int width = cropped.getWidth();
		int height = cropped.getHeight();
		int totalPixels = width * height;

		// Compute gradient magnitudes
		float[][] components = GradientMapAnalyzer.computeGradientComponents(fp.duplicate().convertToFloatProcessor(), false);
		float[] gradMag = GradientMapAnalyzer.computeGradientMagnitude(components[0], components[1], fp, false);

		// Determine gradient cutoff (from CDF)
		float gradientCutoff = computeGradientCutoff(gradMag, width, height, gradientThreshold);

		// Smooth image using the standard kernel
		FloatProcessor smoothedFp = fp.duplicate().convertToFloatProcessor();
		GradientMapAnalyzer.applyDefaultSmoothing(smoothedFp);
		float[] smoothed = (float[]) smoothedFp.getPixels();

		// Extract subset of smoothed pixels where gradient < cutoff (excluding the 1-pixel border)
		ArrayList<Double> subsetList = new ArrayList<>();  
		for (int y = 1; y < height - 1; y++) {
			for (int x = 1; x < width - 1; x++) {
				int idx = y * width + x;
				if (gradMag[idx] < gradientCutoff) {
					subsetList.add((double) smoothed[idx]);
				}
			}
		}
		double[] subset = toArray(subsetList);

		// Compute CDF-based thresholds on the subset to find I_low and I_high
		ThresholdData thresholds = ThresholdAnalyzer.computeThresholds(
			subset, thrMinContrast, thrMaxContrast, nbins
		);
		double iLow = thresholds.getMinThreshold();
		double iHigh = thresholds.getMaxThreshold();

		return new ContrastData(iLow, 
			iHigh,
			darkCount,
			gradientThreshold,
			totalPixels,
			subsetList.size(),
			subset,
			thrMinContrast,
			thrMaxContrast,
			thresholds
		);
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
    static float computeGradientCutoff(float[] gradients, int width, int height, double gradientThreshold) {
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
