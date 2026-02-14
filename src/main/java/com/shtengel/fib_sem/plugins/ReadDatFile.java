package com.shtengel.fib_sem.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.process.ImageProcessor;
import java.io.File;
import java.io.IOException;
import java.text.Normalizer;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import com.shtengel.fib_sem.core.DatFileProcessor;
import com.shtengel.fib_sem.data.HeaderData;
import com.shtengel.fib_sem.data.FileData;

/**
 * ImageJ/Fiji front-end plugin for reading and visualizing custom FIB-SEM .dat files.
 *
 * <ul>
 *   <li>Parses the file header and image data</li>
 *   <li>Creates one ImageJ window per detector/channel</li>
 *   <li>Applies spatial and temporal calibration</li>
 *   <li>Attaches full metadata to Image &gt; Show Info</li>
 * </ul>
 *
 * <p>
 * This class is intentionally lightweight and UI-focused. All file parsing
 * and image construction logic is delegated to {@link DatFileProcessor}.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins > FIB-SEM > Read .dat")
public class ReadDatFile implements Command {

	/** First output image produced by the plugin. **/
	@Parameter(type = ItemIO.OUTPUT)
	private ImagePlus outputImage;
	
	@Parameter
	private LogService logService;
	
	/**
     * Entry point invoked when the plugin is selected from the menu.
     * Will immediately prompt the user to select a .dat file.
     */
	@Override
	public void run() {
		// Prompt user for file immediately
        String path = IJ.getFilePath("Select FIB-SEM .dat file");
        if (path == null) {
            // User cancelled dialog
            logService.info("File selection cancelled by user");
            return;
        }
        
        File inputFile = new File(path);

        // Validate file
        if (!inputFile.exists() || !inputFile.isFile()) {
            IJ.error("Invalid File", "Selected file does not exist.");
            return;
        }
        if (!path.toLowerCase().endsWith(".dat")) {
            IJ.error("Unsupported File Type", "Please select a .dat file.");
            return;
        }
		
        try {
            DatFileProcessor processor = new DatFileProcessor(logService);
            processDatFile(processor, inputFile);

        } catch (IOException e) {
            logService.error("Failed to process file: " + inputFile.getName(), e);
            IJ.error("I/O Error", e.getMessage());

        } catch (Exception e) {
            logService.error("Unexpected error processing file", e);
            IJ.error("Unexpected Error", e.getMessage());
        }
	}
	
	/**
     * Reads a FIB-SEM .dat file and generates one or more ImageJ windows depending on the number of detector channels.
     *
     * @param processor  the core DAT file processor
     * @param inputFile  the selected .dat file
     * @throws IOException if file parsing fails
     */
	private void processDatFile(DatFileProcessor processor, File inputFile) throws IOException {
		logService.info("Processing .dat file: " + inputFile.getName());
		
		// Read the .dat file
		FileData fileData = processor.readDatFile(inputFile);
		HeaderData header = fileData.getHeaderData();
		logService.info(header.toString());
		
		ImagePlus imp = processor.processImageData(fileData);
		
		// Multi-channel data: open one window per channel
		if (imp.getStackSize() > 1) {
		
			for (int i = 1; i <= imp.getStackSize(); i++) {
				ImageProcessor ip = imp.getStack().getProcessor(i);
				String label = imp.getStack().getSliceLabel(i);
				
				String windowTitle = getBaseName(inputFile) + "_" + getSafeFileString(label);
				ImagePlus channelImage = new ImagePlus(windowTitle, ip);
				
				setCalibration(channelImage, header);
				addHeaderInfo(channelImage, header, label);
				setFileInfo(channelImage, inputFile);
				
				channelImage.show();
				
				if (i == 1) {
					outputImage = channelImage;
				}
			}
		} else {
			// Single channel data
			setCalibration(imp, header);
			addHeaderInfo(imp, header, "Single Channel");
			setFileInfo(imp, inputFile);
			imp.show();
			outputImage = imp;
		}
		
		logService.info("Successfully processed .dat file");
	}
	
	/**
     * Applies spatial and temporal calibration from the .dat header.
     *
     * @param imp     the ImagePlus to calibrate
     * @param header 	parsed .dat header
     */
	private void setCalibration(ImagePlus imp, HeaderData header) {
		Calibration cal = imp.getCalibration();
		
		// Set pixel size in nm
		double pixelSize = header.getPixelSize();
		if (pixelSize > 0) {
			cal.pixelWidth = pixelSize;
			cal.pixelHeight = pixelSize;
			cal.setUnit("nm");
		}
		
		// Set time calibration
		cal.setTimeUnit("ms");
		if (header.getTimeStep() > 0) {
			cal.frameInterval = header.getTimeStep() * 1000; // Convert to ms
		}
	}
	
	/**
     * Populates Image > Show Info with comprehensive metadata extracted from the .dat file header.
     *
     * @param imp          target image
     * @param header       .dat header data
     * @param channelLabel channel identifier
     */
	private void addHeaderInfo(ImagePlus imp, HeaderData header, String channelLabel) {
		StringBuilder info = new StringBuilder();
		
		// Channel and file information
		info.append("=== File ===\n");
		info.append("File: ").append(imp.getTitle()).append("\n");
		info.append("File Version: ").append(header.getFileVersion()).append("\n");
		info.append("Magic Number: ").append(header.getFileMagicNum()).append("\n");
		info.append("Sample: ").append(header.getSampleId()).append("\n");
		info.append("Notes: ").append(header.getNotes()).append("\n\n");
		
		info.append(header.toString());
		
		// Set the info string to be displayed in Image > Show Info
		imp.setProperty("Info", info.toString());
	}
	
	/**
	 * Get file base name, assuming file name contains a single "."
	 */
	public static String getBaseName(File file) {
		String filename = file.getName();
		int index = filename.lastIndexOf('.');
		if (index != -1) {
	        filename = filename.substring(0, index);
	    }
	    return filename;
	}
	
	/**
	 * Convert sample ID into a cross-platform compatible "sanitized" string
	 */
	public static String getSafeFileString(String id) {
		if (id== null) {
			return "null";
		}
		
		String normalized = Normalizer.normalize(id, Normalizer.Form.NFKD);
		String ascii = normalized.replaceAll("[^\\p{ASCII}]", "");
		
		// Replace invalid filename characters, then collapse whitespace and trim
		String sanitized = ascii
				.replaceAll("[\\\\/:*?\"<>|]", "_")
		        .replaceAll("\\s+", "")
		        .replaceAll("^\\.+", "")
		        .replaceAll("[. ]+$", "")
		        .trim();
		
		return (sanitized.isEmpty()) ? "file" : sanitized;
	}
	
	private void setFileInfo(ImagePlus imp, File sourceFile) {
		FileInfo fi = new FileInfo();
		fi.directory = sourceFile.getParent() + File.separator;
		fi.fileName = imp.getTitle();
		imp.setFileInfo(fi);
	}
}