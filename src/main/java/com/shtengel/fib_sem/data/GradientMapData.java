package com.shtengel.fib_sem.data;

import ij.process.ImageProcessor;

/**
 * Container for a gradient magnitude map computed from an image.
 *
 * <p>The gradient magnitude at each pixel represents the local rate of
 * intensity change, computed via central finite differences. Optionally
 * smoothed and/or normalized by local intensity.</p>
 *
 * @see com.shtengel.fib_sem.core.GradientMapAnalyzer
 */
public class GradientMapData {
	/** Gradient magnitude image (float precision). */
    private final ImageProcessor gradient;

    /**
     * @param gradient gradient magnitude image
     */
    public GradientMapData(ImageProcessor gradient) {
        this.gradient = gradient;
    }

    /**
     * @return gradient magnitude {@link ImageProcessor}
     */
    public ImageProcessor getGradient() {
        return gradient;
    }
}
