package com.shtengel.fib_sem.util;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;

public class FigBuilder {

	public enum Type {
		THRESHOLD,
		RESOLUTION
	}
	
	private final static int TITLE_HEIGHT = 50;
	private final static int TITLE_FONT_SIZE = 32;
	private final static Font TITLE_FONT = new Font("SansSerif", Font.PLAIN, TITLE_FONT_SIZE);
	
	public static boolean save(ImagePlus figure, String filepath) {
		if (figure == null || filepath == null || filepath.isEmpty()) {
			return false;
		}
		
		try {
			IJ.saveAs(figure, getFormat(filepath), filepath);
			IJ.log("Figure saved: " + filepath);
			return true;
		} catch (Exception e) {
			return false;
		}
	}
	
	private static String getFormat(String filepath) {
		String lower = filepath.toLowerCase();
		
		if (lower.endsWith(".png")) {
			return "PNG";
		} else if (lower.endsWith(".tif") || lower.endsWith(".tiff")) {
			return "TIFF";
		} else if (lower.endsWith(".jpg") || lower.endsWith(".jpeg")) {
			return "JPEG";
		}
		
		// Default to PNG
		return "PNG";
	}
	
	public static void createAndSave(ImagePlus imp, String title, String filepath) {
		if (imp == null) {
			return;
		}
		
		// Get processor and create composite figure (with title)
		ImageProcessor ip = imp.getProcessor();
		ImageProcessor figureIp = createTitledProcessor(ip, title, TITLE_HEIGHT);
		
		// Create new ImagePlus with the figure
		ImagePlus figure = new ImagePlus(title, figureIp);
		
		save(figure, filepath);
	}

	public static void createAndSave(Plot plot, String title, String filepath) {
		if (plot == null) {
			return;
		}
		
		// Get processor and create composite figure (with title)
		ImagePlus plotImp = plot.getImagePlus();
		ImageProcessor plotIp = plotImp.getProcessor();
		ImageProcessor figureIp = createTitledProcessor(plotIp, title, TITLE_HEIGHT);
		
		// Create new ImagePlus with the figure
		ImagePlus figure = new ImagePlus(title, figureIp);
		
		save(figure, filepath);
	}
	
	private static ImageProcessor createTitledProcessor(ImageProcessor contentIp,
														String title,
														int titleHeight) {
		int contentWidth = contentIp.getWidth();
		int contentHeight = contentIp.getHeight();
		int totalHeight = contentHeight + titleHeight;
		
		// Create a color processor for the full figure
		ColorProcessor figureIp = new ColorProcessor(contentWidth, totalHeight);
		
		// Convert content to BufferedImage for drawing (title addition)
		BufferedImage contentImg = contentIp.getBufferedImage();
		BufferedImage figureImg = figureIp.getBufferedImage();
		
		Graphics2D g2d = figureImg.createGraphics();
		
		// Enable antialiasing for text
		g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, 
		                     RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		g2d.setRenderingHint(RenderingHints.KEY_RENDERING, 
		                     RenderingHints.VALUE_RENDER_QUALITY);
		
		// Draw title background and text
		g2d.setColor(Color.WHITE);
		g2d.fillRect(0, 0, contentWidth, titleHeight);
		g2d.setColor(Color.BLACK);
		g2d.setFont(TITLE_FONT);
		
		// Determine text position (centered horizontally, vertically centered in title area)
		java.awt.FontMetrics fm = g2d.getFontMetrics();
		int textWidth = fm.stringWidth(title);
		int textHeight = fm.getAscent();
		int textX = (contentWidth - textWidth) / 2;
		int textY = (titleHeight + textHeight) / 2 - fm.getDescent();
		
		g2d.drawString(title, textX, textY);
		
		// Draw content image
		g2d.drawImage(contentImg, 0, titleHeight, null);
		
		g2d.dispose();
		
		// Convert back to ColorProcessor
		ColorProcessor result = new ColorProcessor(figureImg);
		
		return result;
	}
}