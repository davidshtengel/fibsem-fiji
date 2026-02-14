package com.shtengel.fib_sem.data;

import java.util.Objects;

/**
 * Container for file data, including
 * - Header information
 * - Raw image data
 */
public class FileData {
	
	private final HeaderData datHeader;
	private final byte[] imageData;
	private final String fileName;
	
	/**
	 * Instantiate a FileData instance.
	 * 
	 * @param headerData		the parsed header information
	 * @param imageData		the raw image data bytes
	 * @param fileName		the source file name
	 */
	public FileData(HeaderData datHeader, byte[] imageData, String fileName) {
		this.datHeader = Objects.requireNonNull(datHeader, "Header data cannot be null");
        this.imageData = Objects.requireNonNull(imageData, "Image data cannot be null");
        this.fileName = Objects.requireNonNull(fileName, "File name cannot be null");
	}
	
	public HeaderData getHeaderData() { return datHeader; }
	public byte[] getImageData() { return imageData; }
	public String getFileName() { return fileName; }
	
	/**
	 * Gets the size of the image data in bytes.
     * 
     * @return the data size
     */
    public long getDataSize() {
        return imageData.length;
    }
    
    @Override
    public String toString() {
        return String.format("FileData[fileName=%s, dataSize=%d bytes, header=%s]",
                fileName, imageData.length, datHeader);
    }

}