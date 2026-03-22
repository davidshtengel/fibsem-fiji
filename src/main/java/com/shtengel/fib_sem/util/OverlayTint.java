package com.shtengel.fib_sem.util;

/**
 * Static utility to tint grayscale pixel values, currently used for exclusion-mask overlays.
 *
 * <p>Each tint method takes a 0–255 grayscale value and returns a packed RGB integer (tint/grayscale blend is 50/50), 
 * to preserve image detail and indicate particular exclusion categories.
 *
 * <p>Color convention (matches both Contrast preview and Noise Statistics mask):
 * <ul>
 *   <li><b>Red</b> — below minimum intensity threshold</li>
 *   <li><b>Cyan</b> — above maximum intensity threshold</li>
 *   <li><b>Blue</b> — excluded by gradient magnitude (determined by gradient threshold)</li>
 *   <li><b>Grayscale</b> — retained / included pixels</li>
 * </ul>
 * </p>
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

	/** Returns a red-tinted packed RGB (below minimum intensity threshold). */
	public static int redTint(int gray) {
		return ((gray + 255) / 2 << 16) | (gray / 2 << 8) | gray / 2;
	}

	/** Returns a cyan-tinted packed RGB (above maximum intensity threshold). */
	public static int cyanTint(int gray) {
		return (gray / 2 << 16) | ((gray + 255) / 2 << 8) | (gray + 255) / 2;
	}

	/** Returns a blue-tinted packed RGB (excluded by high gradient magnitude). */
	public static int blueTint(int gray) {
		return (gray / 2 << 16) | (gray / 2 << 8) | (gray + 255) / 2;
	}
}
