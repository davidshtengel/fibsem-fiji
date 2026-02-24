package com.shtengel.fib_sem.core;

import com.shtengel.fib_sem.data.GradientMapData;

import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Utility class for computing gradient magnitude maps from ImageJ
 * {@link ImageProcessor} instances.
 *
 * <p>
 * Gradients are computed using central finite differences on a float-converted
 * copy of the input image. If an ROI is set on the input processor, computation
 * is restricted to the cropped ROI region and the output dimensions match the ROI.
 * </p>
 */
public class GradientMapAnalyzer {
	/**
     * Computes a gradient magnitude map from the given image.
     *
     * <p>
     * Processing pipeline:
     * <ol>
     *   <li>Crop to ROI if present</li>
     *   <li>Convert to float precision</li>
     *   <li>Optional smoothing via a normalized 3x3 kernel</li>
     *   <li>Central-difference gradient estimation (x and y)</li>
     *   <li>Gradient magnitude computation</li>
     *   <li>Optional normalization by local intensity</li>
     * </ol>
     * </p>
     *
     * @param ip                source {@link ImageProcessor}, optionally with ROI
     * @param performSmoothing  apply pre-gradient smoothing
     * @param normalize         normalize magnitude by local pixel intensity
     * @return gradient magnitude image wrapped in {@link GradientMapData}
     */
    public static GradientMapData computeGradientMap(ImageProcessor ip, 
    		boolean performSmoothing, 
    		boolean normalize) {

        ImageProcessor ipToProcess = (ip.getRoi() != null) ? ip.crop() : ip;
        FloatProcessor fp = ipToProcess.convertToFloatProcessor();
        int width = fp.getWidth();
        int height = fp.getHeight();
    	
        float[][] components = computeGradientComponents(fp, performSmoothing);
        float[] gradX = components[0];
        float[] gradY = components[1];
		float[] absGrad = computeGradientMagnitude(gradX, gradY, fp, normalize);
                
        FloatProcessor gradProcessor = new FloatProcessor(width, height, absGrad);
        
        return new GradientMapData(gradProcessor);
    }

    /**
     * Applies default 3x3 smoothing kernel to the entire image.
     * This respects ROIs in that only ROI pixels will be smoothed.
     *
     * If an ROI is present on the processor, convolution respects ROI bounds.
     *
     * @param ip ImageProcessor to smooth (modified in place)
     */
    public static void applyDefaultSmoothing(ImageProcessor ip) {
        double st = 1.0 / Math.sqrt(2.0);
        
        float[] kernel = new float[]{
            (float) st, 1f, (float) st,
                    1f, 1f,         1f,
            (float) st, 1f, (float) st
        };
        
        // Normalize kernel weights
        float sum = 0;
        for (float v : kernel) {
            sum += v;
        }
        for (int i = 0; i < kernel.length; i++) {
            kernel[i] /= sum;
        }
        
        new Convolver().convolve(ip, kernel, 3, 3);
    }

    /**
     * Computes gradient components (X and Y) from an image processor.
     * Returns separate X and Y gradient arrays for directional analysis.
     *
     * <p>
     * This method uses the same gradient computation logic as {@link #computeGradientMap}
     * but returns the directional components instead of magnitude.
     * </p>
     *
     * @param fp source image processor (will be cropped to ROI if present)
     * @param performSmoothing whether to apply default smoothing before gradient computation
     * @return float[2][] where [0] is gradX and [1] is gradY (flattened arrays)
     */
    public static float[][] computeGradientComponents(FloatProcessor fp, boolean performSmoothing) {
        int width = fp.getWidth();
        int height = fp.getHeight();
        
        // Optional smoothing
        if (performSmoothing) {
            applyDefaultSmoothing(fp);
        }

        // Gradient component buffers
        float[] gradX = new float[width * height];
        float[] gradY = new float[width * height];
        
        // Central finite differences (skip 1-pixel border)
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                int idx = y * width + x;

                float dx = fp.getf(x + 1, y) - fp.getf(x - 1, y);
                float dy = fp.getf(x, y + 1) - fp.getf(x, y - 1);
                
                gradX[idx] = dx * 0.5f;
                gradY[idx] = dy * 0.5f;
            }
        }
        
        return new float[][] {gradX, gradY};
    }

    /**
     * Computes gradient magnitude from pre-computed X and Y components.
     *
     * @param gradX X gradient component (flattened array)
     * @param gradY Y gradient component (flattened array)
	 * @param fp        source float processor (used for normalization)
	 * @param normalize whether to divide magnitude by local intensity
     * @return gradient magnitude array (same length as input arrays)
     */
    public static float[] computeGradientMagnitude(float[] gradX, float[] gradY, FloatProcessor fp, boolean normalize) {
        if (gradX.length != gradY.length) {
            throw new IllegalArgumentException("Gradient component arrays must have the same length.");
        }

        float[] magnitudes = new float[gradX.length];
        // Gradient magnitude (optionally normalized by intensity)
        for (int y = 0; y < fp.getHeight(); y++) {
            for (int x = 0; x < fp.getWidth(); x++) {
                int idx = y * fp.getWidth() + x;
                
                float magnitude = (float) Math.sqrt(
                    gradX[idx] * gradX[idx] + gradY[idx] * gradY[idx]
                );
                
                if (normalize) {
                    float v = fp.getf(x, y);
                    magnitudes[idx] = (v != 0) ? magnitude / v : 0f;
                } else {
                    magnitudes[idx] = magnitude;
                }
            }
        }
        return magnitudes;
    }
}
