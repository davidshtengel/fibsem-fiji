package com.shtengel.fib_sem.util;

import java.awt.Color;

/**
 * Static utility to tint grayscale pixel values, used for exclusion-mask overlays.
 *
 * <p>The {@link #tint(int, Color)} method takes a 0–255 grayscale value and a tint color,
 * returning a packed RGB integer that is a 50/50 blend of the tint with the grayscale value.
 * This preserves image detail while indicating exclusion categories.
 *
 * <p>Tint colors for each category are defined in {@link Col} (e.g. {@code Col.THR_MIN},
 * {@code Col.THR_MAX}, {@code Col.GRADIENT}).
 */
public class OverlayTint {

	private OverlayTint() {}

	/** Precomputed grayscale LUT: {@code GRAY_LUT[g] = (g << 16) | (g << 8) | g}. */
	public static final int[] GRAY_LUT = new int[256];
	static {
		for (int g = 0; g < 256; g++) {
			GRAY_LUT[g] = (g << 16) | (g << 8) | g;
		}
	}

	/**
	 * Maps a floating-point pixel value to a 0–255 grayscale index.
	 *
	 * @param value pixel intensity
	 * @param min   minimum intensity in the image
	 * @param range intensity range (max - min); must be &gt; 0
	 * @return clamped grayscale value in [0, 255]
	 */
	public static int toGray(float value, float min, float range) {
		return (int) Math.max(0, Math.min(255, 255.0f * (value - min) / range));
	}

	/** Returns a pure grayscale packed RGB for retained/included pixels. */
	public static int grayscale(int gray) {
		return GRAY_LUT[gray];
	}

	/**
	 * Returns a tinted packed RGB: 50/50 blend of the given color with grayscale.
	 *
	 * @param gray  grayscale value in [0, 255]
	 * @param color tint color (typically a {@link Col} constant)
	 * @return packed RGB integer
	 */
	public static int tint(int gray, Color color) {
		int r = (gray + color.getRed()) / 2;
		int g = (gray + color.getGreen()) / 2;
		int b = (gray + color.getBlue()) / 2;
		return (r << 16) | (g << 8) | b;
	}
}
