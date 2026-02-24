package com.shtengel.fib_sem.core;

import com.shtengel.fib_sem.data.HeaderData;
import com.shtengel.fib_sem.data.FileData;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.scijava.log.LogService;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

/**
 * Processor for reading FIB-SEM .dat files.
 * 
 * This class handles the complete workflow of reading custom .dat file format used in FIB-SEM imaging,
 * including parsing binary headers, extracting metadata, and converting raw pixel data into ImageJ-compatible
 * image stacks
 * 
 * <h2>File Format Structure</h2>
 * <ul>
 *   <li><b>Header (1024 bytes)</b>: Contains metadata including resolution, bit depth, channel configuration,
 *   SEM/FIB parameters, stage positions, detector settings, and scaling factors</li>
 * 	 <li><b>Image Data (variable):</b> Raw pixel data stored as big-endian interleaved channels</li>
 * </ul>
 * 
 * <h2>Supported Data Types</h2>
 * <ul>
 *   <li>8-bit unsigned (0-255)</li>
 *   <li>16-bit signed (-32768 to 32767)</li>
 *   <li>32-bit floating point</li>
 * </ul>
 * 
 * <h2>Multi-Channel Support</h2>
 * <p>The processor handles up to 2 active channels (AI1/AI2) corresponding to different detectors
 * (e.g., In-lens SE, ESB). Channel data is stored interleaved in the file and de-interleaved
 * during processing.</p>
 */
public class DatFileProcessor {
	
	/** Magic number identifying valid .dat files (0xD3D3D3D2 in hex) */
    private static final long EXPECTED_MAGIC_NUMBER = 3555587570L;
    
    /** Size of the header section in bytes */
    private static final int HEADER_SIZE = 1024;
    
    private final LogService logService;
    
    /**
     * Constructs a new DatFileProcessor with the specified logging service.
     * 
     * @param logService the logging service for recording processing information and errors
     */
    public DatFileProcessor(LogService logService) {
    	this.logService = logService;
    }
    
    /**
     * Reads a complete .dat file including header metadata and image data.
     * 
     * <ol>
     * 	<li>Opens the file using RandomAccessFile FileChannel for efficient reading</li>
     * 	<li>Reads and parses the 1024-byte header</li>
     * 	<li>Validates the file format using the magic number</li>
     * 	<li>Reads the raw image data starting at byte 1024</li>
     * 	<li>Logs key file properties for debugging</li>
     * </ol>
     * 
     * @param file	the .dat file to read
     * @return a {@link FileData} object containing the parsed header and raw image data
     * @throws IOException if the file cannot be read, is not a valid .dat file (wrong magic number), or if the header/data sections are incomplete.
     * @throws IllegalArgumentException if the file is null
     */
    public FileData readDatFile(File file) throws IOException {
    	try(RandomAccessFile raf = new RandomAccessFile(file, "r"); 
    		FileChannel channel = raf.getChannel()) {
    		
    		// Read header (first 1024 bytes)
    		HeaderData header = readDatHeader(channel);
    		
    		// Verify magic number
    		if (header.getFileMagicNum() != EXPECTED_MAGIC_NUMBER) {
    			throw new IOException("Invalid .dat file: magic number mismatch. Expected 3555587570, got " + header.getFileMagicNum());
    		}
    		
    		// Read image data starting at byte 1024
    		channel.position(1024);
    		byte[] imageData = readImageData(channel, header.getDataSize());

    		// Log file info
    		logService.info("Read .dat file: " + file.getName());
    		logService.info("  Version: " + header.getFileVersion());
    		logService.info("  Dimensions: " + header.getXRes() + "x" + header.getYRes());
    		logService.info("  Bit depth: " + header.getBitDepth());
    		logService.info("  Channels: " + header.getChannelCount());
    		logService.info("  Active channels - AI1:" + header.getAi1() + " AI2:" + header.getAi2());
    			
    		return new FileData(header, imageData, file.getName());
    	}
    }
	
    /**
     * Processes raw image data into an ImageJ ImagePlus object with proper channel separation
     * 
     * <p>This method handles both single-channel and multi-channel data. For multi-channel data,
     * it de-interleaves the channels and creates separate slices in the ImageStack. The appropriate
     * ImageProcessor type (ByteProcessor, ShortProcessor, or FloatProcessor) is selected based on
     * the bit depth specified in the header.</p>
     * 
     * <p><b>NOTE:</b> 16-bit data is treated as signed shorts (-32768 to 32767) The display range is explicitly set to ensure
     * correct visualization of negative values.</p>
     * 
     * @param fileData 		the file data containing header information and raw pixel bytes
     * @return an ImagePlus object ready for display or further processing in ImageJ
     * @throws IllegalArgumentException if the header type is not DatHeaderData or if an unsupported bit depth is encountered (must be 8, 16, or 32)
     */
    public ImagePlus processImageData(FileData fileData) {
    	HeaderData header = fileData.getHeaderData();
    	int width, height; 
    	int depth = 1;  // Default depth for single slice
    	int bitDepth;
    	int channelCount = 1;
    	boolean ai1 = false, ai2 = false;
    	String detA = null, detB = null;
    		
		width = header.getXRes();
		height = header.getYRes();
		bitDepth = (header.getEightBit() == 1) ? 8 : 16;
		channelCount = header.getChannelCount();
		ai1 = header.getAi1();
		ai2 = header.getAi2();
		detA = header.getDetA();
		detB = header.getDetB();
    	
		ImageStack stack;
    	if (channelCount > 1) {
    		stack = createMultiChannelStack(fileData.getImageData(), width, height, bitDepth, ai1, ai2, detA, detB);
    	} else { 
    		stack = createImageStack(fileData.getImageData(), width, height, depth, bitDepth); 
    	}
    		    	
    	ImagePlus imp = new ImagePlus(fileData.getFileName(), stack);
    	return imp;
    }
    
    /**
     * Creates an ImageStack from multi-channel .dat data with proper channel de-interleaving.
     * 
     * <ol>
     * 	<li>Determines active channel names from detector configuration</li>
     * 	<li>Allocates separate pixel arrays for each channel</li>
     * 	<li>De-interleaves the raw data into channel-specific arrays</li>
     * 	<li>Creates appropriately configured ImageProcessors for each channel</li>
     * 	<li>Assembles the processors into a multi-slice ImageStack</li>
     * </ol>
     * 
     * @param data 		the raw interleaved pixel data bytes
     * @param width 	image width in pixels
     * @param height 	image height in pixels
     * @param bitDepth 	bit depth per pixel (8, 16, or 32)
     * @param ai1 		whether channel 1 (detector A) is active
     * @param ai2 		whether channel 2 (detector B) is active
     * @param detA 		name/identifier of detector A (e.g., "InLens")
     * @param detB 		name/identifier of detector B (e.g., "ESB")
     * @return an ImageStack containing one slice per active channel with appropriate labels
     * @throws IllegalArgumentException if bitDepth is not 8, 16, or 32
     */
    private ImageStack createMultiChannelStack(byte[] data, int width, int height, int bitDepth, boolean ai1, boolean ai2, String detA, String detB) {
	    ImageStack stack = new ImageStack(width, height);
	    	
	    String[] channelNames = getActiveChannelNames(ai1, ai2, detA, detB);
	    int channels = channelNames.length;
	    int pixelCount = width * height;
	    	
	    Object[] channelPixels = new Object[channels];
	    
        for (int c = 0; c < channels; c++) {
        	if (bitDepth == 8) {
        		channelPixels[c] = new byte[pixelCount];
        	} else if (bitDepth == 16) {
        		channelPixels[c] = new short[pixelCount];
        	} else if (bitDepth == 32) {
        		channelPixels[c] = new float[pixelCount];
        	} else {
        		throw new IllegalArgumentException("Unsupported bit depth: "+bitDepth);
        	}
        }
        ByteBuffer bb = ByteBuffer.wrap(data).order(ByteOrder.BIG_ENDIAN);
        	
        // De-interleave
        for (int i = 0; i < pixelCount; i++) {
        	for (int c = 0; c < channels; c++) {
        		if (bitDepth == 8) {
        			((byte[]) channelPixels[c])[i] = bb.get();
        		} else if (bitDepth == 16) {
        			((short[]) channelPixels[c])[i] = bb.getShort();
        		} else if (bitDepth == 32) {
        			((float[]) channelPixels[c])[i] = bb.getFloat();
        		}	
        	}
        }
        
        // Build ImageStack
        for (int c = 0; c < channels; c++) {
            ImageProcessor ip;
            if (bitDepth == 8) {
                ip = new ByteProcessor(width, height, (byte[]) channelPixels[c], null);
            } else if (bitDepth == 16) {
            	// Convert signed shorts to floats for proper display
                short[] signedPixels = (short[]) channelPixels[c];
                float[] floatPixels = new float[signedPixels.length];
                for (int i = 0; i < signedPixels.length; i++) {
                    floatPixels[i] = signedPixels[i];
                }
                ip = new FloatProcessor(width, height, floatPixels);
            } else {
                ip = new FloatProcessor(width, height, (float[]) channelPixels[c]);
            }
            stack.addSlice(channelNames[c], ip);
        }
                
        return stack;
    }
    
    /**
     * Reads and parses the 1024-byte binary header from a .dat file.
     * 
     * <p>The header contains comprehensive metadata organized into several sections:</p>
     * <ul>
     *   <li><b>File metadata:</b> magic number, version, file type, software date</li>
     *   <li><b>Image properties:</b> resolution, bit depth, channels, pixel size, scan rate</li>
     *   <li><b>SEM parameters:</b> magnification, working distance, EHT, current, aperture</li>
     *   <li><b>FIB parameters:</b> focus, probe, current, rotation, alignment, stigmation</li>
     *   <li><b>Stage position:</b> X, Y, Z, tilt, rotation, milling position</li>
     *   <li><b>Detector configuration:</b> active channels, detector names, brightness, contrast</li>
     *   <li><b>Scaling factors:</b> offset, gain, and higher-order calibration terms per channel</li>
     *   <li><b>Milling parameters:</b> resolution, size, line time, PID control (version 5+)</li>
     * </ul>
     * 
     * <p>The method handles multiple file format versions (1-9+) with version-specific field layouts.
     * All multi-byte values are read in big-endian byte order.</p>
     * 
     * @param channel the FileChannel positioned at the start of the file
     * @return a DatHeaderData object containing all parsed metadata
     * @throws IOException if the header cannot be fully read (file too short) or if reading fails
     */
    private HeaderData readDatHeader(FileChannel channel) throws IOException {
    	ByteBuffer buffer = ByteBuffer.allocate(1024);
        buffer.order(ByteOrder.BIG_ENDIAN);

        // Ensure complete header read
        int bytesRead = 0;
        while (bytesRead < HEADER_SIZE) {
            int read = channel.read(buffer);
            if (read == -1) {
                throw new IOException(String.format(
                    "Header incomplete: expected %d bytes, got %d (file size=%d)",
                    HEADER_SIZE, bytesRead, channel.size()
                ));
            }
            bytesRead += read;
        }
        buffer.flip();

        HeaderData header = new HeaderData();

        // File metadata
        header.setFileMagicNum(buffer.getInt(0) & 0xFFFFFFFFL);
        short fileVersion = buffer.getShort(4);
        header.setFileVersion(fileVersion);
        header.setFileType(buffer.getShort(6));

        byte[] swDateBytes = new byte[10];
        buffer.get(swDateBytes);
        header.setSWDate(readNullTerminatedString(swDateBytes));
        
        header.setTimeStep(buffer.getDouble(24));
        int chanNum = buffer.get(32) & 0xFF; 
        header.setChannelCount(chanNum);
        header.setEightBit(buffer.get(33) & 0xFF);

        // Scaling factors for each channel
        buffer.position(36);
        float[][] scaling;
        if (fileVersion == 1) {
        	scaling = new float[chanNum][4];
        	for (int c = 0; c < chanNum; c++) {
        		for (int r = 0; r < 4; r++) {
                    scaling[c][r] = (float) buffer.getDouble();
                }
            }
        } else if (fileVersion >=2 && fileVersion <= 6) {
        	scaling = new float[chanNum][4];
        	for (int c = 0; c < chanNum; c++) {
        		for (int r = 0; r < 4; r++) {
                    scaling[c][r] = buffer.getFloat();
                }
            }
        } else {
        	// Version 7+ stores scaling flat as [coeff0_ch0, coeff1_ch0, ..., coeff3_ch1];
        	// read the 8 floats then transpose into [channel][coefficient] layout
        	float[] rawScaling = new float[8];
        	for (int i = 0; i < 8; i++) {
        		rawScaling[i] = buffer.getFloat();
        	}
        	scaling = new float[chanNum][4];
        	for (int c = 0; c < chanNum; c++) {
        		for (int r = 0; r < 4; r++) {
        			scaling[c][r] = rawScaling[r + c * 4];
        		}
        	}
        }
        if (chanNum == 2 && (fileVersion == 8 || fileVersion == 9)) {
        	scaling[0][3] = (float) 0.1;
        	scaling[1][3] = (float) 0.1;
        }
        header.setScaling(scaling);

        // Flags and first pixel coordinates (version 9+)
        if (fileVersion > 8) {
            header.setRestartFlag(buffer.get(68) & 0xFF);
            header.setStageMove(buffer.get(69) & 0xFF);
            header.setFirstPixelX(buffer.getInt(70));
            header.setFirstPixelY(buffer.getInt(74));
        }

        // Resolution acquisition parameters
        header.setXRes(buffer.getInt(100));
        header.setYRes(buffer.getInt(104));

        if (fileVersion >= 1 && fileVersion <= 3) {
        	// Early version format
            header.setOversampling(buffer.get(108) & 0xFF);
            header.setAiDelay(buffer.getShort(109));
            header.setZeissScanSpeed(buffer.get(111) & 0xFF);
            header.setScanRate(buffer.getFloat(112));
            header.setFramelineRampdown(buffer.getDouble());
            header.setXmin(buffer.getDouble(128));
            header.setXmax(buffer.getDouble(136));
            header.setDetMin(-10);
            header.setDetMax(10);
        } else {
        	// Version 4+ format
    		header.setOversampling(buffer.getShort(108) & 0xFFFF);
    		header.setScanRate(buffer.getFloat(112));
    		header.setFramelineRampdown(buffer.getFloat(116));
    		header.setXmin(buffer.getFloat(120));
    		header.setXmax(buffer.getFloat(124));
    		header.setDetMin(buffer.getFloat(128));
    		header.setDetMax(buffer.getFloat(132));
    		header.setDecimatingFactor(buffer.getShort(136));
        }
        
        header.setSaveOversamples(buffer.get(138) != 0);
        header.setAi1(buffer.get(151) != 0);
        header.setAi2(buffer.get(152) != 0);
        
        // Notes field
        buffer.position(180);
        byte[] notes = new byte[200];
        buffer.get(notes);
        header.setNotes(readNullTerminatedString(notes));
        
        // Sample ID
        if (fileVersion > 8) {
            buffer.position(155);
            byte[] sampleIdBytes = new byte[25];
            buffer.get(sampleIdBytes);
            header.setSampleId(readNullTerminatedString(sampleIdBytes).trim());
        } else {
        	// Extract from notes field for older versions
    		String idFromNotes = readNullTerminatedString(notes).split(",")[0];
    		header.setSampleId(idFromNotes);
        }
        
        // Detector names
        buffer.position(380);
        byte[] rawDetA = new byte[10];
        byte[] rawDetB = new byte[18];
        buffer.get(rawDetA);
        buffer.get(rawDetB);
        header.setDetA(readNullTerminatedString(rawDetA));
        header.setDetB(readNullTerminatedString(rawDetB));
        
        // SEM/FIB parameters (version-dependent field sizes)
        if (fileVersion == 1 || fileVersion == 2) {
    		header.setMagnification(buffer.getDouble(408));
    		header.setPixelSize(buffer.getDouble(416));
    		header.setWorkingDistance(buffer.getDouble(424));
    		header.setEHT(buffer.getDouble(432));
    		
    		header.setSEMApr(buffer.get(440));
    		header.setHighCurrent(buffer.get(441) != 0);
    		header.setSEMCurrent(buffer.getDouble(448));
    		header.setSEMRotation(buffer.getDouble(456));
    		header.setChamberVacuum(buffer.getDouble(464));
    	    header.setGunVacuum(buffer.getDouble(472));
    	    
    	    header.setSEMStigX((float) buffer.getDouble(480));
    	    header.setSEMStigY((float) buffer.getDouble(488));
    	    header.setSEMAlignX((float) buffer.getDouble(496));
    	    header.setSEMAlignY((float) buffer.getDouble(504));
    	    
    	    header.setStageX((float) buffer.getDouble(512));
    	    header.setStageY((float) buffer.getDouble(520));
    	    header.setStageZ((float) buffer.getDouble(528));
    	    header.setStageT((float) buffer.getDouble(536));
    	    header.setStageR((float) buffer.getDouble(544));
    	    header.setStageM((float) buffer.getDouble(552));
    	    
    	    header.setBrightnessA((float) buffer.getDouble(560));
    	    header.setContrastA((float) buffer.getDouble(568));
    	    header.setBrightnessB((float) buffer.getDouble(576));
    	    header.setContrastB((float) buffer.getDouble(584));
    	    header.setMode(buffer.get(600));
    	    
    	    header.setFIBFocus((float) buffer.getDouble(608));
    	    header.setFIBProbe(buffer.get(616));
    	    header.setFIBCurrent((float) buffer.getDouble(624));
    	    header.setFIBRotation((float) buffer.getDouble(632));
    	    header.setFIBAlignX((float) buffer.getDouble(640));
    	    header.setFIBAlignY((float) buffer.getDouble(648));
    	    header.setFIBStigX((float) buffer.getDouble(656));
    	    header.setFIBStigY((float) buffer.getDouble(664));
    	    header.setFIBShiftX((float) buffer.getDouble(672));
    	    header.setFIBShiftY((float) buffer.getDouble(680));
        } else {
    		header.setMagnification(buffer.getFloat(460));
            header.setPixelSize(buffer.getFloat(464));
            header.setWorkingDistance(buffer.getFloat(468));
            header.setEHT(buffer.getFloat(472));

            header.setSEMApr(buffer.get(480));
            header.setHighCurrent(buffer.get(481) != 0);
            header.setSEMCurrent(buffer.getFloat(490));
            header.setSEMRotation(buffer.getFloat(494));
            header.setChamberVacuum(buffer.getFloat(498));
            header.setGunVacuum(buffer.getFloat(502));

            header.setSEMShiftX(buffer.getFloat(510));
            header.setSEMShiftY(buffer.getFloat(514));
            header.setSEMStigX(buffer.getFloat(518));
            header.setSEMStigY(buffer.getFloat(522));
            header.setSEMAlignX(buffer.getFloat(526));
            header.setSEMAlignY(buffer.getFloat(530));

            header.setStageX(buffer.getFloat(534));
            header.setStageY(buffer.getFloat(538));
            header.setStageZ(buffer.getFloat(542));
            header.setStageT(buffer.getFloat(546));
            header.setStageR(buffer.getFloat(550));
            header.setStageM(buffer.getFloat(554));

            header.setBrightnessA(buffer.getFloat(560));
            header.setContrastA(buffer.getFloat(564));
            header.setBrightnessB(buffer.getFloat(568));
            header.setContrastB(buffer.getFloat(572));
            header.setMode(buffer.get(600));

            header.setFIBFocus(buffer.getFloat(604));
            header.setFIBProbe(buffer.get(608));
            header.setFIBCurrent(buffer.getFloat(620));
            header.setFIBRotation(buffer.getFloat(624));
            header.setFIBAlignX(buffer.getFloat(628));
            header.setFIBAlignY(buffer.getFloat(632));
            header.setFIBStigX(buffer.getFloat(636));
            header.setFIBStigY(buffer.getFloat(640));
            header.setFIBShiftX(buffer.getFloat(644));
            header.setFIBShiftY(buffer.getFloat(648));
        }
        
        // Milling parameters (version 5+)
        if (fileVersion > 4) {
            header.setMillingXResolution(Integer.toUnsignedLong(buffer.getInt(652)));
            header.setMillingYResolution(Integer.toUnsignedLong(buffer.getInt(656)));
            header.setMillingXSize(buffer.getFloat(660));
            header.setMillingYSize(buffer.getFloat(664));
            header.setMillingULAngle(buffer.getFloat(668));
            header.setMillingURAngle(buffer.getFloat(672));
            header.setMillingLineTime(buffer.getFloat(676));
            header.setFIBFOV(buffer.getFloat(680));
            header.setMillingLinesPerImage(Short.toUnsignedInt(buffer.getShort(684)));
            header.setMillingPID(buffer.get(686) != 0);
            header.setMillingPIDMeasured(buffer.get(689));
            header.setMillingPIDTarget(buffer.getFloat(690));
            header.setMillingPIDTargetSlope(buffer.getFloat(694));
            header.setMillingPIDP(buffer.getFloat(698));
            header.setMillingPIDI(buffer.getFloat(702));
            header.setMillingPIDD(buffer.getFloat(706));

            byte[] machineId = new byte[30];
            buffer.position(800);
            buffer.get(machineId);
            header.setMachineID(readNullTerminatedString(machineId));
        }

        // Additional parameters (version 6+)
        if (fileVersion > 5) {
            header.setTemperature(buffer.getFloat(850));
            header.setFaradayCupCurrent(buffer.getFloat(854));
            header.setFIBSpecimenCurrent(buffer.getFloat(858));
            header.setBeamDump1Current(buffer.getFloat(862));
            header.setSEMSpecimenCurrent(buffer.getFloat(866));
            header.setMillingYVoltage(buffer.getFloat(870));
            header.setFocusIndex(buffer.getFloat(874));
            header.setFIBSliceNum(Integer.toUnsignedLong(buffer.getInt(878)));
        }
        
        // Extended parameters (version 8+)
        if (fileVersion > 7) {
            header.setBeamDump2Current(buffer.getFloat(882));
            header.setMillingCurrent(buffer.getFloat(886));
        }

        return header;
    }
    
    /**
     * Determines active channel names from detector configuration.
     * 
     * <p>Generates channel names based on which detectors are active and their
     * associated detector names from the header. If detector names are not specified or are empty,
     * default names "Channel A" and "Channel B" are used.</p>
     * 
     * <p><b>Note:</b> Channels C and D (DetC/DetD) are deprecated in current file versions and
     * are not supported by this method.</p>
     * 
     * @param ai1 	whether detector A is active
     * @param ai2 	whether detector B is active
     * @param detA 	detector A name/identifier (e.g., "InLens", etc.)
     * @param detB 	detector B name/identifier (e.g., "ESB", etc.)
     * @return array of channel names for active channels only
     */
    private String[] getActiveChannelNames(boolean ai1, boolean ai2, String detA, String detB) {
        java.util.List<String> names = new java.util.ArrayList<>();
        if (ai1) { names.add((detA != null && !detA.trim().isEmpty()) ? detA.trim() : "Channel A"); }
        if (ai2) { names.add((detB != null && !detB.trim().isEmpty()) ? detB.trim() : "Channel B"); }
        
        return names.toArray(new String[0]);
    }
    
    /**
     * Reads the raw image data section from the file channel.
     * 
     * <p>Reads the complete image data payload into memory. The data size is determined from
     * the header metadata. This method ensures all bytes are read even if the channel returns
     * data in multiple chunks.</p>
     * 
     * @param channel 		the FileChannel positioned at the start of image data (byte 1024)
     * @param dataSize 		the number of bytes to read as specified in the header
     * @return byte array containing the raw image data
     * @throws IOException if the file size is larger than Integer.MAX_VALUE (>2GB),
     *                     if the file ends before all expected data is read,
     *                     or if a read error occurs
     */
    private byte[] readImageData(FileChannel channel, long dataSize) throws IOException {
    	if (dataSize > Integer.MAX_VALUE) {
            throw new IOException("File too large to load into memory: " + dataSize);
        }
    	byte[] data = new byte[(int) dataSize];
        ByteBuffer buffer = ByteBuffer.wrap(data);

        // Read until buffer is full
        while (buffer.hasRemaining()) {
            int n = channel.read(buffer);
            if (n == -1) {
                throw new EOFException("Unexpected end of file while reading image data");
            }
        }
        
        return data;
    }
    
    /**
	 * Creates an ImageStack from single-channel or non-interleaved multi-slice data.
	 * 
	 * <p>Processes image data that is stored sequentially (non-interleaved), creating one slice
	 * per z-position. Each slice is processed independently using the appropriate bit depth.</p>
	 * 
	 * @param data 		the raw sequential pixel data
	 * @param width 	image width in pixels
	 * @param height 	image height in pixels
	 * @param depth 	number of z-slices (typically 1 for 2D FIB-SEM data)
	 * @param bitDepth 	bit depth per pixel (8, 16, or 32)
	 * @return an ImageStack containing the processed slices
	 */
    private ImageStack createImageStack(byte[] data, int width, int height, int depth, int bitDepth) {
    	ImageStack stack = new ImageStack(width, height);
        
    	int bytesPerPixel = bitDepth / 8;
    	int sliceSize = width * height * bytesPerPixel;
    		
    	for (int z = 0; z < depth; z++) {
   			ImageProcessor ip = createProcessor(data, z * sliceSize, width, height, bitDepth);
   			stack.addSlice("Slice " + (z + 1), ip);
   		}       
    	return stack;
   	}
    
    /**
	 * Creates an ImageProcessor of the appropriate type based on bit depth.
	 * 
	 * @param data 		the raw pixel data bytes
	 * @param offset 	byte offset to start reading from
	 * @param width 	image width in pixels
	 * @param height 	image height in pixels
	 * @param bitDepth 	bit depth (8, 16, or 32)
	 * @return the appropriate ImageProcessor subclass for the bit depth
	 * @throws IllegalArgumentException if bitDepth is not 8, 16, or 32
	 */    	
   	private ImageProcessor createProcessor(byte[] data, int offset, int width, int height, int bitDepth) {
   		if (bitDepth == 8) {
   			return create8BitProcessor(data, offset, width, height);
   		} else if (bitDepth == 16) {
       		return create16BitProcessor(data, offset, width, height);
        } else if (bitDepth == 32) {
       		return create32BitProcessor(data, offset, width, height);
        } else {
       		throw new IllegalArgumentException("Unsupported bit depth: " + bitDepth);
        }
   	}
    
   	/**
	 * Creates an 8-bit ByteProcessor from raw data.
	 * 
	 * @param data 		the raw pixel data
	 * @param offset 	byte offset to start reading from
	 * @param width 	image width in pixels
	 * @param height 	image height in pixels
	 * @return a ByteProcessor containing the 8-bit unsigned pixel data
	 */
    private ImageProcessor create8BitProcessor(byte[] data, int offset, int width, int height) {
    	byte[] pixels = new byte[width * height];
        System.arraycopy(data, offset, pixels, 0, pixels.length);
        
        return new ij.process.ByteProcessor(width, height, pixels, null);
    }
    
   	/**
	 * Creates an 16-bit ShortProcessor from raw data with proper signed value handling.
	 * 
	 * <p>Reads 16-bit values in big-endian byte order and treats them as signed integers.</p>
	 * 
	 * @param data 		the raw pixel data
	 * @param offset 	byte offset to start reading from
	 * @param width 	image width in pixels
	 * @param height 	image height in pixels
	 * @return a ShortProcessor configured for signed 16-bit data display
	 */ 
    private ImageProcessor create16BitProcessor(byte[] data, int offset, int width, int height) {
    	short[] signedPixels = new short[width * height];
    	ByteBuffer.wrap(data, offset, signedPixels.length * 2)
   				.order(ByteOrder.BIG_ENDIAN)  // .dat files use big-endian
   				.asShortBuffer()
   				.get(signedPixels);
    	
    	// Convert to float to preserve signed values
        float[] floatPixels = new float[signedPixels.length];
        for (int i = 0; i < signedPixels.length; i++) {
            floatPixels[i] = signedPixels[i];
        }
        
        return new FloatProcessor(width, height, floatPixels);
    }
    
    /**
	 * Creates an 32-bit FloatProcessor from raw data.
	 * 
	 * <p>Reads 32-bit floating-point values in big-endian byte order. This format may be
	 * used for pre-calibrated or processed data.</p>
	 * 
	 * @param data 		the raw pixel data
	 * @param offset 	byte offset to start reading from
	 * @param width 	image width in pixels
	 * @param height 	image height in pixels
	 * @return a FloatProcessor containing the 32-bit floating-point pixel data
	 */
    private ImageProcessor create32BitProcessor(byte[] data, int offset, int width, int height) {
    	float[] pixels = new float[width * height];
    	ByteBuffer.wrap(data, offset, pixels.length * 4)
    			.order(ByteOrder.BIG_ENDIAN)  // .dat files use big-endian
   				.asFloatBuffer()
  				.get(pixels);
    	
   		return new FloatProcessor(width, height, pixels);
    }
    	
    /**
     * Extracts a null-terminated string from a byte array.
     * 
     * <p>Reads bytes until a null terminator (0x00) is encountered or the end of the array
     * is reached. This is used for reading C-style strings from the binary header.</p>
     * 
     * @param bytes the byte array containing a null-terminated string
     * @return the extracted string (excluding the null terminator)
     */
   	private String readNullTerminatedString(byte[] bytes) {
   		int length = 0;
   		while (length < bytes.length && bytes[length] != 0) {
    		length++;
   		}
   		
   		return new String(bytes, 0, length);
    }
}