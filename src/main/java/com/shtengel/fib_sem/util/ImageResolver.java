package com.shtengel.fib_sem.util;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.process.ImageProcessor;

public class ImageResolver {
	private static ImagePlus lastSourceImage = null;
	
	private ImageResolver() {}
	
	public static ImagePlus resolveSourceImage() {
		ImagePlus imp = WindowManager.getCurrentImage();
		
		if(isValidSourceImage(imp)) {
			lastSourceImage = imp;
			return imp;
		}
		
		if(isValidSourceImage(lastSourceImage) && lastSourceImage.isVisible()) {
			return lastSourceImage;
		}
		
		IJ.error("No valid source image available");
		return null;
	}
	
	public static ImagePlus resolveSourceImage(ImagePlus candidate) {
	    if (isValidSourceImage(candidate)) {
	        lastSourceImage = candidate;
	        return candidate;
	    }

	    ImagePlus current = WindowManager.getCurrentImage();
	    if (isValidSourceImage(current)) {
	        lastSourceImage = current;
	        return current;
	    }

	    if (isValidSourceImage(lastSourceImage) && lastSourceImage.isVisible()) {
	        return lastSourceImage;
	    }

	    return null;
	}
	
	private static boolean isValidSourceImage(ImagePlus imp) {
        if (imp == null) return false;
        if (imp.getType() == ImagePlus.COLOR_RGB) return false;
        if (imp.getTitle().contains(":")) {return false;}
        return !imp.getTitle().startsWith("Plot");
    }
	
	/**
	 * Helper to get base filename without extension
	 */
	public static String getBaseName(ImagePlus imp) {
	    String title = imp.getTitle();
	    int dot = title.lastIndexOf('.');
	    return (dot > 0) ? title.substring(0, dot) : title;
	}
	
	public static ImageProcessor cropToRoiIfPresent(ImageProcessor ip) {
		return (ip.getRoi() != null ? ip.crop() : ip);
	}
}