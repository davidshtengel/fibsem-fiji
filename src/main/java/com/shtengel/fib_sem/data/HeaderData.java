package com.shtengel.fib_sem.data;

/**
 * Header data for custom .dat file format used in FIB-SEM imaging.
 * Contains metadata about image dimensions, scaling, and acquisition parameters.
 * Based on Shan Xu's binary format with 1024-byte header.
 */
public class HeaderData {
    
	// File/format metadata
	private long fileMagicNum;
    private short fileVersion;
    private short fileType;
    private String machineID;
    private long fileLength;

    public long getFileMagicNum() { return fileMagicNum; }
    public void setFileMagicNum(long fileMagicNum) { this.fileMagicNum = fileMagicNum; }

    public short getFileVersion() { return fileVersion; }
    public void setFileVersion(short fileVersion) { this.fileVersion = fileVersion; }

    public short getFileType() { return fileType; }
    public void setFileType(short fileType) { this.fileType = fileType; }

    public String getMachineID() { return machineID; }
    public void setMachineID(String machineID) { this.machineID = machineID; }

    public long getFileLength() { return fileLength; }
    public void setFileLength(long fileLength) { this.fileLength = fileLength; }
    
    // Acquisition timing/control
    private String swDate;
    private double timeStep;
    private int restartFlag;
    private int stageMove;
    private byte mode;	// 0=SEM, 1=FIB, 2=Milling, etc.

    public String getSWDate() { return swDate; }
    public void setSWDate(String swDate) { this.swDate = swDate; }

    public double getTimeStep() { return timeStep; }
    public void setTimeStep(double timeStep) { this.timeStep = timeStep; }

    public int getRestartFlag() { return restartFlag; }
    public void setRestartFlag(int restartFlag) { this.restartFlag = restartFlag; }

    public int getStageMove() { return stageMove; }
    public void setStageMove(int stageMove) { this.stageMove = stageMove; }

    public byte getMode() { return mode; }
    public void setMode(byte mode) { this.mode = mode; }
	
    // Channel/signal configuration
    private int channelCount;
    private int eightBit;
    private boolean saveOversamples;
    private int decimatingFactor;
    private float[][] scaling;

    public int getChannelCount() { return channelCount; }
    public void setChannelCount(int channelCount) { this.channelCount = channelCount; }

    public int getEightBit() { return eightBit; }
    public void setEightBit(int eightBit) { this.eightBit = eightBit; }

    public boolean isSaveOversamples() { return saveOversamples; }
    public void setSaveOversamples(boolean saveOversamples) { this.saveOversamples = saveOversamples; }

    public int getDecimatingFactor() { return decimatingFactor; }
    public void setDecimatingFactor(int decimatingFactor) { this.decimatingFactor = decimatingFactor; }

    public float[][] getScaling() {
        if (scaling == null) return null;
        float[][] copy = new float[scaling.length][];
        for (int i = 0; i < scaling.length; i++) {
            copy[i] = scaling[i].clone();
        }
        return copy;
    }

    public void setScaling(float[][] scaling) {
        this.scaling = scaling;
    }
    
    // Scan geometry/resolution
    private int firstPixelX;
    private int firstPixelY;
    private int xResolution;
    private int yResolution;
    private int oversampling;
    private double xmin, xmax, detMin, detMax;
    public int getXRes() { return xResolution; }
    public int getYRes() { return yResolution; }
    public int getOversampling() { return oversampling; }
    public int getFirstPixelX() { return firstPixelX; }
    public int getFirstPixelY() { return firstPixelY; }
    public double getXmin() { return xmin; }
    public double getXmax() { return xmax; }
    public double getDetMin() { return detMin; }
    public double getDetMax() { return detMax; }
    public void setXRes(int xResolution) { this.xResolution = xResolution; }
    public void setYRes(int yResolution) { this.yResolution = yResolution; }
    public void setOversampling(int oversampling) { this.oversampling = oversampling; }
    public void setFirstPixelX(int firstPixelX) { this.firstPixelX = firstPixelX; }
    public void setFirstPixelY(int firstPixelY) { this.firstPixelY = firstPixelY; }
    public void setXmin(double xmin) { this.xmin = xmin; }
    public void setXmax(double xmax) { this.xmax = xmax; }
    public void setDetMin(double detMin) { this.detMin = detMin; }
    public void setDetMax(double detMax) { this.detMax = detMax; }
    
    // Scan timing
    private int zeissScanSpeed;
    private double scanRate;
    private double framelineRampdown;
    public int getZeissScanSpeed() { return zeissScanSpeed; }
    public double getScanRate() { return scanRate; }
    public double getFramelineRampdown() { return framelineRampdown; }
    public void setZeissScanSpeed(int zeissScanSpeed) { this.zeissScanSpeed = zeissScanSpeed; }
    public void setScanRate(double scanRate) { this.scanRate = scanRate; }
    public void setFramelineRampdown(double framelineRampdown) { this.framelineRampdown = framelineRampdown; }

    // Inputs
    private boolean ai1;
    private boolean ai2;
    private short aiDelay;
    public boolean getAi1() { return ai1; }
    public boolean getAi2() { return ai2; }
    public short getAiDelay() { return aiDelay; }
    public void setAi1(boolean ai1) { this.ai1 = ai1; }
    public void setAi2(boolean ai2) { this.ai2 = ai2; }
    public void setAiDelay(short aiDelay) { this.aiDelay = aiDelay; }
    
    // Sample metadata
    private String sampleId;
    private String notes;
    public String getSampleId() { return sampleId; }
    public String getNotes() { return notes; }    
    public void setSampleId(String sampleId) { this.sampleId = sampleId; }
    public void setNotes(String notes) { this.notes = notes; }
    
    // SEM parameters
    private double magnification;
    private double pixelSize;
    private double workingDistance;
    private double eht;
    private int semAperture;
    private boolean highCurrent;
    private double semCurrent;
    private double semRotation;
    public double getMagnification() { return magnification; }
    public double getPixelSize() { return pixelSize; }
    public double getWorkingDistance() { return workingDistance; }
    public double getEHT() { return eht; }
    public int getSEMApr() { return semAperture; }
    public boolean isHighCurrent() { return highCurrent; }
    public double getSEMCurrent() { return semCurrent; }
    public double getSEMRotation() { return semRotation; }
    public void setMagnification(double magnification) { this.magnification = magnification; }
    public void setPixelSize(double pixelSize) { this.pixelSize = pixelSize; }
    public void setWorkingDistance(double workingDistance) { this.workingDistance = workingDistance; }
    public void setEHT(double eht) { this.eht = eht; }
    public void setSEMApr(int semAperture) { this.semAperture = semAperture; }
    public void setHighCurrent(boolean highCurrent) { this.highCurrent = highCurrent; }
    public void setSEMCurrent(double semCurrent) { this.semCurrent = semCurrent; }
    public void setSEMRotation(double semRotation) { this.semRotation = semRotation; }
    
    // SEM alignment, shift, and stigmation
    private float semAlignX;
    private float semAlignY;
    private float semShiftX;
    private float semShiftY;
    private float semStigX;
    private float semStigY;
    public void setSEMAlignX(float semAlignX) { this.semAlignX = semAlignX; }
	public void setSEMAlignY(float semAlignY) { this.semAlignY = semAlignY; }
	public void setSEMShiftX(float semShiftX) { this.semShiftX = semShiftX; }
	public void setSEMShiftY(float semShiftY) { this.semShiftY = semShiftY; }
	public void setSEMStigX(float semStigX) { this.semStigX = semStigX; }
	public void setSEMStigY(float semStigY) { this.semStigY = semStigY; }
    
    // Vacuum and environment
    private double chamberVacuum;
    private double gunVacuum;
    private double temperature;	// (F)
    public void setChamberVacuum(double chamberVacuum) { this.chamberVacuum = chamberVacuum; }
	public void setGunVacuum(double gunVacuum) { this.gunVacuum = gunVacuum; }
	public void setTemperature(double temperature) { this.temperature = temperature; }
    
    // Stage position/orientation
    private float stageX; 	// (mm)
    private float stageY;
    private float stageZ;
    private float stageT;	// (degrees)
    private float stageR;
    private float stageM;
	public void setStageX(float stageX) { this.stageX = stageX; }
	public void setStageY(float stageY) { this.stageY = stageY; }
	public void setStageZ(float stageZ) { this.stageZ = stageZ; }
	public void setStageT(float stageT) { this.stageT = stageT; }
	public void setStageR(float stageR) { this.stageR = stageR; }
	public void setStageM(float stageM) { this.stageM = stageM; }
    
    // Detectors
    private String detA;
    private String detB;
    public String getDetA() { return detA; }
    public String getDetB() { return detB; }
    public void setDetA(String detA) { this.detA = detA; }
	public void setDetB(String detB) { this.detB = detB; }
    
    // Detector settings
    private float brightnessA;
    private float contrastA;
    private float brightnessB;
    private float contrastB;
    public void setBrightnessA(float brightnessA) { this.brightnessA = brightnessA; }
	public void setContrastA(float contrastA) { this.contrastA = contrastA; }
	public void setBrightnessB(float brightnessB) { this.brightnessB = brightnessB; }
	public void setContrastB(float contrastB) { this.contrastB = contrastB; }

    // FIB beam parameters
    private float fibFocus;
    private byte fibProbe;
    private float fibCurrent;
    private float fibRotation;
    public void setFIBFocus(float fibFocus) { this.fibFocus = fibFocus; }
	public void setFIBProbe(byte fibProbe) { this.fibProbe = fibProbe; }
	public void setFIBCurrent(float fibCurrent) { this.fibCurrent = fibCurrent; }
	public void setFIBRotation(float fibRotation) { this.fibRotation = fibRotation; }

    private float fibAlignX;
    private float fibAlignY;
    private float fibStigX;
    private float fibStigY;
    private float fibShiftX;
    private float fibShiftY;
	public void setFIBAlignX(float fibAlignX) { this.fibAlignX = fibAlignX; }
	public void setFIBAlignY(float fibAlignY) { this.fibAlignY = fibAlignY; }
	public void setFIBShiftX(float fibShiftX) { this.fibShiftX = fibShiftX; }
	public void setFIBShiftY(float fibShiftY) { this.fibShiftY = fibShiftY; }
	public void setFIBStigX(float fibStigX) { this.fibStigX = fibStigX; }
	public void setFIBStigY(float fibStigY) { this.fibStigY = fibStigY; }

    private long fibSliceNum;
    private float fibFOV;           // (um)

    public void setFIBSliceNum(long fibSliceNum) { this.fibSliceNum = fibSliceNum; }
	public void setFIBFOV(float fibFOV) { this.fibFOV = fibFOV; }
    
    // Milling geometry and timing
    private long millingXResolution;
    private long millingYResolution;
	public void setMillingXResolution(long millingXResolution) { this.millingXResolution = millingXResolution; }
	public void setMillingYResolution(long millingYResolution) { this.millingYResolution = millingYResolution; }

    private float millingXSize;     // in um
    private float millingYSize;    
    public void setMillingXSize(float millingXSize) { this.millingXSize = millingXSize; }
    public void setMillingYSize(float millingYSize) { this.millingYSize = millingYSize; }

    private float millingULAngle;
    private float millingURAngle;
    public void setMillingULAngle(float millingULAngle) { this.millingULAngle = millingULAngle; }
    public void setMillingURAngle(float millingURAngle) { this.millingURAngle = millingURAngle; }
    
    private float millingLineTime;
    private float millingLinesPerImage;
    private float millingYVoltage;
    public void setMillingLineTime(float millingLineTime) { this.millingLineTime = millingLineTime; }
    public void setMillingLinesPerImage(float millingLinesPerImage) { this.millingLinesPerImage = millingLinesPerImage; }	
    public void setMillingYVoltage(float millingYVoltage) { this.millingYVoltage = millingYVoltage; }
    
    // Milling PID control
    private boolean millingPID;
    private byte millingPIDMeasured;
    private float millingPIDTarget;
    private float millingPIDTargetSlope;
    
    private float millingPID_P;
    private float millingPID_I;
    private float millingPID_D;
    
    // Electrical & beam current monitoring
    private float semSpecimenCurrent;
    private float fibSpecimenCurrent;
    private float faradayCupCurrent;
    
    private float beamDump1Current;
    private float beamDump2Current;
    private float millingCurrent;
    
    // Focus/indexing
    private float focusIndex;
        
    public HeaderData() {
        // Default constructor
    }
    
    
    public int getBitDepth() { return (eightBit == 1) ? 8 : 16; }

    public long getDataSize() {
        long elementsPerChannel = (long) xResolution * yResolution;
        long totalElements = elementsPerChannel * channelCount;
        if (saveOversamples) {
            totalElements *= oversampling;
        }
        return totalElements * ((eightBit == 1 )? 1 : 2);
    }

	public void setFIBSpecimenCurrent(float fibSpecimenCurrent) { this.fibSpecimenCurrent = fibSpecimenCurrent; }


	// Milling PID
	public void setMillingPID(boolean millingPID) { this.millingPID = millingPID; }
	public void setMillingPIDMeasured(byte millingPIDMeasured) { this.millingPIDMeasured = millingPIDMeasured; }
	public void setMillingPIDTarget(float millingPIDTarget) { this.millingPIDTarget = millingPIDTarget; }
	public void setMillingPIDTargetSlope(float millingPIDTargetSlope) { this.millingPIDTargetSlope = millingPIDTargetSlope; }
	public void setMillingPIDP(float millingPID_P) { this.millingPID_P = millingPID_P; }
	public void setMillingPIDI(float millingPID_I) { this.millingPID_I = millingPID_I; }
	public void setMillingPIDD(float millingPID_D) { this.millingPID_D = millingPID_D; }

	public void setFaradayCupCurrent(float faradayCupCurrent) { this.faradayCupCurrent = faradayCupCurrent; }
	public void setBeamDump1Current(float beamDump1Current) { this.beamDump1Current = beamDump1Current; }
	public void setBeamDump2Current(float beamDump2Current) { this.beamDump2Current = beamDump2Current; }
	public void setMillingCurrent(float millingCurrent) { this.millingCurrent = millingCurrent; }
	
	public void setFocusIndex(float focusIndex) { this.focusIndex = focusIndex; }
		
	public void setSEMSpecimenCurrent(float semSpecimenCurrent) { this.semSpecimenCurrent = semSpecimenCurrent; }
    
	
	@Override
	public String toString() {
	    StringBuilder info = new StringBuilder();

	    // Image properties
	    info.append("=== Image Properties ===\n");
	    info.append("Dimensions: ").append(xResolution).append(" x ").append(yResolution).append(" pixels\n");
	    info.append("Bit Depth: ").append((eightBit == 1) ? 8 : 16).append(" bit\n");
	    info.append("Number of Channels: ").append(channelCount).append("\n");
	    info.append("Pixel Size: ").append(String.format("%.2f", pixelSize)).append(" nm\n");
	    info.append("Oversampling: ").append(oversampling).append("\n");
	    info.append("Scan Rate: ").append(String.format("%.3f", scanRate / 1e6)).append(" MHz\n");
	    info.append("Time Step: ").append(String.format("%.6f", timeStep)).append(" s\n\n");

	    // Scaling factors
	    if (scaling != null && scaling.length > 0) {
	        info.append("=== Scaling Factors ===\n");
	        for (int c = 0; c < scaling.length; c++) {
	        	info.append("Channel ").append(c + 1).append(":\n");
            	info.append("    ").append(String.format("%.6f", scaling[c][0])).append("\n");
                info.append("    ").append(String.format("%.6f", scaling[c][1])).append("\n");
                info.append("    ").append(String.format("%.6f", scaling[c][2])).append("\n");
                info.append("    ").append(String.format("%.6f", scaling[c][3])).append("\n");	                
	        }
	        info.append("\n");
	    }

	    // SEM parameters
	    info.append("=== SEM Parameters ===\n");
	    info.append("SEMApr= ").append(semAperture).append("\n");
	    info.append("HighCurrent= ").append(highCurrent).append("\n");
	    info.append("SEMCurr= ").append(semCurrent).append("\n");
	    info.append("SEMSpecimenI= ").append(semSpecimenCurrent).append("\n");
	    info.append("SEMRot= ").append(semRotation).append("\n");
	    info.append("ChamVac= ").append(chamberVacuum).append("\n");
	    info.append("GunVac= ").append(gunVacuum).append("\n");
	    info.append("SEMShiftX= ").append(semShiftX).append("\n");
	    info.append("SEMShiftY= ").append(semShiftY).append("\n");
	    info.append("SEMStiX= ").append(semStigX).append("\n");
	    info.append("SEMStiY= ").append(semStigY).append("\n");
	    info.append("SEMAlnX= ").append(semAlignX).append("\n");
	    info.append("SEMAlnY= ").append(semAlignY).append("\n");
	    info.append("BrightnessA= ").append(brightnessA).append("\n");
	    info.append("ContrastA= ").append(contrastA).append("\n");
	    info.append("BrightnessB= ").append(brightnessB).append("\n");
	    info.append("ContrastB= ").append(contrastB).append("\n");
	    info.append("FocusIndex= ").append(focusIndex).append("\n\n");

	    // Stage position
	    info.append("=== Stage Position ===\n");
	    info.append("StageX= ").append(stageX).append("\n");
	    info.append("StageY= ").append(stageY).append("\n");
	    info.append("StageZ= ").append(stageZ).append("\n");
	    info.append("StageT= ").append(stageT).append("\n");
	    info.append("StageR= ").append(stageR).append("\n");
	    info.append("StageM= ").append(stageM).append("\n\n");

	    // FIB parameters
	    info.append("=== FIB Parameters ===\n");
	    info.append("FIBFocus= ").append(fibFocus).append("\n");
	    info.append("FIBProb= ").append(fibProbe).append("\n");
	    info.append("FIBCurr= ").append(fibCurrent).append("\n");
	    info.append("FIBRot= ").append(fibRotation).append("\n");
	    info.append("FIBAlnX= ").append(fibAlignX).append("\n");
	    info.append("FIBAlnY= ").append(fibAlignY).append("\n");
	    info.append("FIBStiX= ").append(fibStigX).append("\n");
	    info.append("FIBStiY= ").append(fibStigY).append("\n");
	    info.append("FIBShiftX= ").append(fibShiftX).append("\n");
	    info.append("FIBShiftY= ").append(fibShiftY).append("\n");
	    info.append("MillingXResolution= ").append(millingXResolution).append("\n");
	    info.append("MillingYResolution= ").append(millingYResolution).append("\n");
	    info.append("MillingXSize= ").append(millingXSize).append("\n");
	    info.append("MillingYSize= ").append(millingYSize).append("\n");
	    info.append("MillingULAng= ").append(millingULAngle).append("\n");
	    info.append("MillingURAng= ").append(millingURAngle).append("\n");
	    info.append("MillingLineTime= ").append(millingLineTime).append("\n");
	    info.append("MillingLinesPerImage= ").append(millingLinesPerImage).append("\n");
	    info.append("FIBFOV (um)= ").append(fibFOV).append("\n");
	    info.append("MillingPIDOn= ").append(millingPID).append("\n");
	    info.append("MillingPIDMeasured= ").append(millingPIDMeasured).append("\n");
	    info.append("MillingPIDTarget= ").append(millingPIDTarget).append("\n");
	    info.append("MillingPIDTargetSlope= ").append(millingPIDTargetSlope).append("\n");
	    info.append("MillingPIDP= ").append(millingPID_P).append("\n");
	    info.append("MillingPIDI= ").append(millingPID_I).append("\n");
	    info.append("MillingPIDD= ").append(millingPID_D).append("\n");
	    info.append("SEMSpecimenI= ").append(semSpecimenCurrent).append("\n");
	    info.append("Temperature= ").append(temperature).append("\n");
	    info.append("FaradayCupI= ").append(faradayCupCurrent).append("\n");
	    info.append("FIBSpecimenI= ").append(fibSpecimenCurrent).append("\n");
	    info.append("BeamDump1I= ").append(beamDump1Current).append("\n");
	    info.append("BeamDump2I= ").append(beamDump2Current).append("\n");
	    info.append("MillingYVoltage= ").append(millingYVoltage).append("\n");
	    info.append("FIBSliceNum= ").append(fibSliceNum).append("\n");
	    info.append("MillingI= ").append(millingCurrent).append("\n\n");

	    // Detectors
	    info.append("=== Detectors ===\n");
	    info.append("Detector A: ").append(detA != null ? detA.trim() : "N/A").append("\n");
	    info.append("Detector B: ").append(detB != null ? detB.trim() : "N/A").append("\n");
	    info.append("AI1 (Channel A): ").append(ai1 ? "Active" : "Inactive").append("\n");
	    info.append("AI2 (Channel B): ").append(ai2 ? "Active" : "Inactive").append("\n");

	    return info.toString();
	}
}