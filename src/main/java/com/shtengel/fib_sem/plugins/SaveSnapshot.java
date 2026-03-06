package com.shtengel.fib_sem.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import ij.process.ImageProcessor;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.DatFileProcessor;
import com.shtengel.fib_sem.core.ThresholdAnalyzer;
import com.shtengel.fib_sem.data.FileData;
import com.shtengel.fib_sem.data.HeaderData;
import com.shtengel.fib_sem.data.ThresholdData;
import com.shtengel.fib_sem.util.FigBuilder;
import com.shtengel.fib_sem.util.ImageResolver;

/**
 * Generates a PNG snapshot summarizing a FIB-SEM .dat file.
 *
 * <p>The snapshot contains a parameter table, detector metadata,
 * and scaled preview images for each active channel, laid out
 * in a single ~3300x2400 image. No dialog is shown; the snapshot
 * is saved to the same directory as the source .dat file and
 * opened in ImageJ.</p>
 */
@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Save Snapshot (.dat)")
public class SaveSnapshot implements Command {

	@Parameter
	private LogService logService;

	// Layout constants
	private static final int WIDTH = 3300;

	private static final int MARGIN = 60;
	private static final int TITLE_HEIGHT = 80;
	private static final int NOTES_HEIGHT = 70;

	// Table layout
	private static final int TABLE_TOP = TITLE_HEIGHT + NOTES_HEIGHT + 20;
	private static final int TABLE_COLS = 3;
	private static final int TABLE_ROWS = 4;
	private static final int TABLE_COL_LABEL_WIDTH = 280;
	private static final int TABLE_COL_VALUE_WIDTH = 420;
	private static final int TABLE_CELL_WIDTH = TABLE_COL_LABEL_WIDTH + TABLE_COL_VALUE_WIDTH;
	private static final int TABLE_ROW_HEIGHT = 90;
	private static final int TABLE_WIDTH = TABLE_CELL_WIDTH * TABLE_COLS;
	private static final int TABLE_LEFT = (WIDTH - TABLE_WIDTH) / 2;
	private static final int TABLE_HEIGHT = TABLE_ROW_HEIGHT * TABLE_ROWS;

	// Detector sections start after table
	private static final int DETECTOR_LABEL_HEIGHT = 65;

	// Threshold CDF fractions for display range
	private static final double THR_MIN = 1.0e-3;
	private static final double THR_MAX = 1.0e-3;
	private static final int THR_BINS = 256;

	// Fonts
	private static final Font TITLE_FONT = new Font("SansSerif", Font.PLAIN, 32);
	private static final Font NOTES_FONT = new Font("SansSerif", Font.PLAIN, 36);
	private static final Font TABLE_LABEL_FONT = new Font("SansSerif", Font.PLAIN, 32);
	private static final Font TABLE_VALUE_FONT = new Font("SansSerif", Font.BOLD, 32);
	private static final Font DETECTOR_FONT = new Font("SansSerif", Font.PLAIN, 32);

	@Override
	public void run() {
		// Resolve source image
		ImagePlus sourceImp = ImageResolver.resolveSourceImage();
		if (sourceImp == null) {
			IJ.error("Save Snapshot", "No valid source image is open.");
			return;
		}

		// Get source .dat file from FileInfo
		FileInfo fi = sourceImp.getOriginalFileInfo();
		if (fi == null || fi.directory == null) {
			IJ.error("Save Snapshot", "Cannot determine source .dat file location.\n"
					+ "Please open a .dat file first using 'Read .dat'.");
			return;
		}

		// Find the .dat file in the source directory
		File sourceDir = new File(fi.directory);
		File datFile = findDatFile(sourceDir, sourceImp.getTitle());
		if (datFile == null) {
			IJ.error("Save Snapshot", "Cannot find the source .dat file in:\n" + fi.directory);
			return;
		}

		try {
			DatFileProcessor processor = new DatFileProcessor(logService);
			FileData fileData = processor.readDatFile(datFile);
			HeaderData header = fileData.getHeaderData();
			ImagePlus imp = processor.processImageData(fileData);

			BufferedImage snapshot = renderSnapshot(datFile, header, imp);

			// Save
			String baseName = ReadDatFile.getBaseName(datFile);
			String outputPath = sourceDir.getPath() + File.separator + baseName + "_snapshot.png";
			ImagePlus snapshotImp = new ImagePlus(baseName + "_snapshot", snapshot);
			FigBuilder.save(snapshotImp, outputPath);

			// Open
			snapshotImp.show();

		} catch (IOException e) {
			logService.error("Failed to create snapshot", e);
			IJ.error("Save Snapshot", "Error: " + e.getMessage());
		}
	}

	/**
	 * Finds the source .dat file based on the image title and source directory.
	 */
	private File findDatFile(File sourceDir, String imageTitle) {
		// Try direct match: title might be "filename.dat"
		File direct = new File(sourceDir, imageTitle);
		if (direct.exists() && imageTitle.toLowerCase().endsWith(".dat")) {
			return direct;
		}

		// Title is typically "basename_DetectorName" — strip detector suffix and add .dat
		// Try progressively shorter prefixes
		String name = imageTitle;
		int lastUnderscore = name.lastIndexOf('_');
		while (lastUnderscore > 0) {
			String candidate = name.substring(0, lastUnderscore) + ".dat";
			File f = new File(sourceDir, candidate);
			if (f.exists()) {
				return f;
			}
			name = name.substring(0, lastUnderscore);
			lastUnderscore = name.lastIndexOf('_');
		}

		// Fallback: find any .dat file in the directory
		File[] datFiles = sourceDir.listFiles((dir, n) -> n.toLowerCase().endsWith(".dat"));
		if (datFiles != null && datFiles.length == 1) {
			return datFiles[0];
		}

		return null;
	}

	/**
	 * Renders the full snapshot image.
	 */
	private BufferedImage renderSnapshot(File datFile, HeaderData header, ImagePlus imp) {
		// Compute channel info first so we know how many channels
		int channelCount = imp.getStackSize();
		ChannelInfo[] channels = new ChannelInfo[channelCount];
		for (int i = 0; i < channelCount; i++) {
			channels[i] = buildChannelInfo(header, imp, i);
		}

		// Calculate total height dynamically
		int headerAreaHeight = TABLE_TOP + TABLE_HEIGHT + 30;
		int availableForImages = computeAvailableImageHeight(channelCount);
		int totalHeight = headerAreaHeight;
		for (int i = 0; i < channelCount; i++) {
			totalHeight += DETECTOR_LABEL_HEIGHT;
			int imgH = computeScaledImageHeight(channels[i].ip, availableForImages);
			totalHeight += imgH + 10;
		}
		totalHeight = Math.max(totalHeight, 2400);

		// Recalculate available image height based on actual total
		int nonImageHeight = headerAreaHeight + channelCount * (DETECTOR_LABEL_HEIGHT + 10);
		int maxImgHeight = (totalHeight - nonImageHeight) / Math.max(channelCount, 1);

		BufferedImage image = new BufferedImage(WIDTH, totalHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = image.createGraphics();
		g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

		// White background
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, WIDTH, totalHeight);

		// 1. File path title
		drawTitle(g, datFile.getAbsolutePath());

		// 2. Notes line
		drawNotes(g, header.getNotes());

		// 3. Parameter table
		drawTable(g, header);

		// 4. Detector sections
		int y = headerAreaHeight;
		for (int i = 0; i < channelCount; i++) {
			y = drawDetectorSection(g, channels[i], y, maxImgHeight);
		}

		g.dispose();
		return image;
	}

	private int computeAvailableImageHeight(int channelCount) {
		int nonImageHeight = TABLE_TOP + TABLE_HEIGHT + 30
				+ channelCount * (DETECTOR_LABEL_HEIGHT + 15);
		return (2400 - nonImageHeight) / Math.max(channelCount, 1);
	}

	private int computeScaledImageHeight(ImageProcessor ip, int maxHeight) {
		int w = ip.getWidth();
		int h = ip.getHeight();
		int maxWidth = WIDTH - 2 * MARGIN;
		double scaleW = (double) maxWidth / w;
		double scaleH = (double) maxHeight / h;
		double scale = Math.min(scaleW, scaleH);
		return (int) (h * scale);
	}

	// --- Drawing methods ---

	private void drawTitle(Graphics2D g, String filePath) {
		g.setColor(new Color(0xF0F0F0));
		g.fillRect(0, 0, WIDTH, TITLE_HEIGHT);
		g.setColor(Color.DARK_GRAY);
		g.setFont(TITLE_FONT);
		FontMetrics fm = g.getFontMetrics();
		int textY = (TITLE_HEIGHT + fm.getAscent() - fm.getDescent()) / 2;
		// Right-align or truncate from left if too long
		int textWidth = fm.stringWidth(filePath);
		int textX = Math.max(MARGIN, WIDTH - MARGIN - textWidth);
		g.drawString(filePath, textX, textY);
	}

	private void drawNotes(Graphics2D g, String notes) {
		if (notes == null || notes.trim().isEmpty()) return;
		int y = TITLE_HEIGHT;
		g.setColor(Color.WHITE);
		g.fillRect(0, y, WIDTH, NOTES_HEIGHT);
		g.setColor(Color.BLACK);
		g.setFont(NOTES_FONT);
		FontMetrics fm = g.getFontMetrics();
		int textX = (WIDTH - fm.stringWidth(notes)) / 2;
		int textY = y + (NOTES_HEIGHT + fm.getAscent() - fm.getDescent()) / 2;
		g.drawString(notes, textX, textY);
	}

	private void drawTable(Graphics2D g, HeaderData h) {
		// Table data: [row][col] -> {label, value}
		String[][][] cells = {
			{
				{"Sample ID", h.getSampleId() != null ? h.getSampleId() : "N/A"},
				{"Frame Size", h.getXRes() + " x " + h.getYRes()},
				{"Scan Rate", String.format("%.3f MHz", h.getScanRate() / 1e6)}
			},
			{
				{"Machine ID", h.getMachineID() != null ? h.getMachineID() : "N/A"},
				{"Pixel Size", String.format("%.1f nm", h.getPixelSize())},
				{"Oversampling", String.valueOf(h.getOversampling())}
			},
			{
				{"FileVersion", String.valueOf(h.getFileVersion())},
				{"Working Dist.", String.format("%.3f mm", h.getWorkingDistance())},
				{"FIB Focus", String.format("%.1f  V", h.getFIBFocus())}
			},
			{
				{"Bit Depth", String.valueOf(h.getBitDepth())},
				{"EHT Voltage", String.format("%.3f kV", h.getEHT())},
				{"FIB Probe", String.valueOf(h.getFIBProbe())}
			}
		};

		// Also add SEM Current in the last row's middle cell area
		// (In the screenshot it appears as a second line under EHT Voltage)

		g.setStroke(new BasicStroke(1.0f));

		for (int row = 0; row < TABLE_ROWS; row++) {
			for (int col = 0; col < TABLE_COLS; col++) {
				int x = TABLE_LEFT + col * TABLE_CELL_WIDTH;
				int y = TABLE_TOP + row * TABLE_ROW_HEIGHT;

				// Draw cell border
				g.setColor(Color.LIGHT_GRAY);
				g.drawRect(x, y, TABLE_CELL_WIDTH, TABLE_ROW_HEIGHT);

				// Label (left side of cell)
				g.setColor(Color.DARK_GRAY);
				g.setFont(TABLE_LABEL_FONT);
				FontMetrics fmLabel = g.getFontMetrics();
				int labelX = x + 15;
				int labelY = y + (TABLE_ROW_HEIGHT + fmLabel.getAscent() - fmLabel.getDescent()) / 2;
				g.drawString(cells[row][col][0], labelX, labelY);

				// Value (right side of cell)
				g.setColor(Color.BLACK);
				g.setFont(TABLE_VALUE_FONT);
				FontMetrics fmValue = g.getFontMetrics();
				int valueX = x + TABLE_COL_LABEL_WIDTH + 10;

				if (row == 3 && col == 1) {
					// Special case: EHT Voltage + SEM Current stacked
					int lineH = TABLE_ROW_HEIGHT / 2;
					int valY1 = y + (lineH + fmValue.getAscent() - fmValue.getDescent()) / 2;
					g.drawString(cells[row][col][1], valueX, valY1);
					String semCurStr = String.format("%.3f nA", h.getSEMCurrent() * 1e9);
					int valY2 = y + lineH + (lineH + fmValue.getAscent() - fmValue.getDescent()) / 2;
					g.drawString(semCurStr, valueX, valY2);

					// Also draw "SEM Current" label for the second line
					g.setColor(Color.DARK_GRAY);
					g.setFont(TABLE_LABEL_FONT);
					g.drawString("SEM Current", labelX, valY2);
				} else {
					int valueY = y + (TABLE_ROW_HEIGHT + fmValue.getAscent() - fmValue.getDescent()) / 2;
					g.drawString(cells[row][col][1], valueX, valueY);
				}
			}
		}
	}

	private int drawDetectorSection(Graphics2D g, ChannelInfo ch, int y, int maxImgHeight) {
		// Detector label line
		g.setColor(Color.BLACK);
		g.setFont(DETECTOR_FONT);
		FontMetrics fm = g.getFontMetrics();

		String label = ch.label;
		int labelWidth = fm.stringWidth(label);
		int labelX = (WIDTH - labelWidth) / 2;
		int labelY = y + (DETECTOR_LABEL_HEIGHT + fm.getAscent() - fm.getDescent()) / 2;
		g.drawString(label, labelX, labelY);
		y += DETECTOR_LABEL_HEIGHT;

		// Scale and draw the image
		ImageProcessor ip = ch.ip;
		int imgW = ip.getWidth();
		int imgH = ip.getHeight();
		int maxWidth = WIDTH - 2 * MARGIN;
		double scaleW = (double) maxWidth / imgW;
		double scaleH = (double) maxImgHeight / imgH;
		double scale = Math.min(scaleW, scaleH);
		int scaledW = (int) (imgW * scale);
		int scaledH = (int) (imgH * scale);

		BufferedImage channelImg = createDisplayImage(ip, ch.displayMin, ch.displayMax);
		int imgX = (WIDTH - scaledW) / 2;
		g.drawImage(channelImg, imgX, y, scaledW, scaledH, null);

		y += scaledH + 10;
		return y;
	}

	/**
	 * Creates a grayscale BufferedImage from an ImageProcessor with the given display range.
	 */
	private BufferedImage createDisplayImage(ImageProcessor ip, double min, double max) {
		int w = ip.getWidth();
		int h = ip.getHeight();
		float[] pixels = (float[]) ip.convertToFloatProcessor().getPixels();
		BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_BYTE_GRAY);

		double range = max - min;
		if (range <= 0) range = 1;

		byte[] data = new byte[w * h];
		for (int i = 0; i < pixels.length; i++) {
			double val = (pixels[i] - min) / range;
			val = Math.max(0, Math.min(1, val));
			data[i] = (byte) ((1.0 - val) * 255); // inverted LUT for visual ease
		}
		img.getRaster().setDataElements(0, 0, w, h, data);
		return img;
	}

	// --- Channel info helper ---

	private static class ChannelInfo {
		String label;
		ImageProcessor ip;
		double displayMin;
		double displayMax;
	}

	private ChannelInfo buildChannelInfo(HeaderData header, ImagePlus imp, int channelIndex) {
		ChannelInfo info = new ChannelInfo();

		ImageProcessor ip;
		String detName;
		float brightness, contrast;

		if (imp.getStackSize() > 1) {
			ip = imp.getStack().getProcessor(channelIndex + 1);
			String sliceLabel = imp.getStack().getSliceLabel(channelIndex + 1);
			detName = (sliceLabel != null) ? sliceLabel : ("Channel " + (char)('A' + channelIndex));
		} else {
			ip = imp.getProcessor();
			String detA = header.getDetA();
			detName = (detA != null && !detA.trim().isEmpty()) ? detA.trim() : "Channel A";
		}

		if (channelIndex == 0) {
			brightness = header.getBrightnessA();
			contrast = header.getContrastA();
		} else {
			brightness = header.getBrightnessB();
			contrast = header.getContrastB();
		}

		// Compute thresholds for display range
		float[] pixels = (float[]) ip.convertToFloatProcessor().getPixels();
		ThresholdData thr = ThresholdAnalyzer.computeThresholds(pixels, THR_MIN, THR_MAX, THR_BINS);

		info.ip = ip;
		info.displayMin = thr.getMinThreshold();
		info.displayMax = thr.getMaxThreshold();

		// Build label string matching screenshot format
		String detLabel = (channelIndex == 0) ? "Detector A" : "Detector B";
		info.label = String.format(
				"%s:  %s,  Data Range:  %.1f \u00F7 %.1f with thr_min=%.1e, thr_max=%.1e    (Brightness: %.1f, Contrast: %.1f)",
				detLabel, detName,
				thr.getMinIntensity(), thr.getMaxIntensity(),
				THR_MIN, THR_MAX,
				brightness, contrast
		);

		return info;
	}
}
