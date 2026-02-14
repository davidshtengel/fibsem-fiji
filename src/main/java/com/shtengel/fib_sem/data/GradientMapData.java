package com.shtengel.fib_sem.data;

import ij.process.ImageProcessor;

/**
 * Simple data container for gradient map results.
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
