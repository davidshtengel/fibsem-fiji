package com.shtengel.fib_sem.util;

import ij.IJ;
import ij.ImagePlus;
import java.util.LinkedHashMap;

/**
 * Utility class for persisting plugin parameters as {@link ImagePlus} properties.
 *
 * <p>Parameters are stored with a plugin-specific prefix to avoid key collisions
 * (e.g., {@code "R_lowerBound"} for EdgeTransitions, {@code "N_"} for NoiseStatistics,
 * {@code "C_"} for Contrast). Keys and their prefixes must be consistent between
 * {@code get} and {@code set} calls.</p>
 */
public class ParamPersister {

	private ParamPersister() {}

	public static double get(ImagePlus imp, String key, double def) {
		Object val = imp.getProperty(key);
		return val instanceof Number ? ((Number) val).doubleValue() : def;
	}

	public static int get(ImagePlus imp, String key, int def) {
		Object val = imp.getProperty(key);
		return val instanceof Number ? ((Number) val).intValue() : def;
	}

	public static boolean get(ImagePlus imp, String key, boolean def) {
		Object val = imp.getProperty(key);
		return val instanceof Boolean ? (Boolean) val : def;
	}

	public static float get(ImagePlus imp, String key, float def) {
		Object val = imp.getProperty(key);
		return val instanceof Number ? ((Number) val).floatValue() : def;
	}

	public static void set(ImagePlus imp, String key, Object value) {
		imp.setProperty(key, value);
	}

	/**
	 * Logs a labeled block of parameters to the ImageJ log window.
	 *
	 * @param sectionTitle header label shown at the top and bottom of the block
	 * @param params       ordered map of parameter label to value (use {@link LinkedHashMap} to preserve order)
	 */
	public static void logParams(String sectionTitle, LinkedHashMap<String, Object> params) {
		int dividerLen = Math.max(sectionTitle.length() + 10, 38);
		StringBuilder sb = new StringBuilder(dividerLen);
		for (int i = 0; i < dividerLen; i++) sb.append('-');
		String divider = sb.toString();
		IJ.log("--- " + sectionTitle + " ---");
		for (java.util.Map.Entry<String, Object> entry : params.entrySet()) {
			IJ.log(entry.getKey() + ": " + entry.getValue());
		}
		IJ.log(divider);
	}
}
