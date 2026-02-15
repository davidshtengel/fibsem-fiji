package com.shtengel.fib_sem.util;

import ij.ImagePlus;

public class ParamPersister {

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
}
