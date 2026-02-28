package com.shtengel.fib_sem.core;

import ij.IJ;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.Arrays;

import com.shtengel.fib_sem.data.EdgeTransitionData;
import com.shtengel.fib_sem.data.ThresholdData;
import com.shtengel.fib_sem.util.ImageResolver;

/**
 * Detects edges in FIB-SEM images and measures transition distances
 * along gradient directions.
 *
 * <p>Edge points are identified by gradient magnitude thresholding,
 * filtered by neighbor exclusion, and each transition is measured
 * by finding where the intensity profile crosses configurable
 * lower/upper bounds between local min and max.</p>
 */
public class EdgeTransitionAnalyzer {
	/**
	 * Measures the edge transition distance at a single edge point.
	 *
	 * <p>Extracts an intensity profile along the gradient direction using bilinear
	 * interpolation, locates local min/max within a configurable aperture, validates
	 * against threshold criteria, then interpolates the exact positions where the
	 * profile crosses the lower and upper transition bounds.</p>
	 *
	 * @return a {@link TransitionResult} with errorFlag == 0 on success, or a
	 *         non-zero flag indicating the failure mode:
	 *         <ul>
	 *           <li>1 — local minimum exceeds the minimum intensity criterion</li>
	 *           <li>2 — local maximum falls below the maximum intensity criterion</li>
	 *           <li>4 — crossing point not found or transition out of range</li>
	 *           <li>8 — transition distance below low limit</li>
	 *           <li>16 — transition distance is NaN</li>
	 *         </ul>
   	 */
    public static EdgeTransitionData analyzeEdgeTransitions(
            ImageProcessor ip,
            double lowerBound,
            double upperBound,
            double pixelSize,
            int subsetSize,
            double sectionLength,
            int minMaxAperture,
            double transitionLowLimit,
            double transitionHighLimit,
            double neighbourExclusionRadius,
            double thrMinCriterion,
            double thrMaxCriterion,
            boolean excludeCenter,
            double centerExclusionRadius,
            double gradThreshold) {
    	
        // Crop to ROI, if present
        ImageProcessor ipToProcess = ImageResolver.cropToRoiIfPresent(ip);
        FloatProcessor fp = ipToProcess.convertToFloatProcessor();
        
        int width = fp.getWidth();
        int height = fp.getHeight();

        // Compute gradients
        float[][] gradComponents = GradientMapAnalyzer.computeGradientComponents(fp, false);
        float[] gradX = gradComponents[0];
        float[] gradY = gradComponents[1];
        float[] absGrad = GradientMapAnalyzer.computeGradientMagnitude(gradX, gradY, fp, false);
        float gradMax = computeGradientThreshold(absGrad, width, height, gradThreshold);
        
        ArrayList<Integer> xPtsList = new ArrayList<>();
        ArrayList<Integer> yPtsList = new ArrayList<>();
        ArrayList<Float> gradPtsList = new ArrayList<>();

        int halfSubset = subsetSize / 2;
        for (int y = halfSubset; y < height - halfSubset; y++) {
            for (int x = halfSubset; x < width - halfSubset; x++) {
                int idx = y * width + x;
                if(absGrad[idx] > gradMax) {
                    xPtsList.add(x);
                    yPtsList.add(y);
                    gradPtsList.add(absGrad[idx]);
                }
            }
        }
        
        int initialEdgeCount = xPtsList.size();
        IJ.log("Initial edge count is "+initialEdgeCount);
        if (initialEdgeCount == 0) { return createEmptyResult(0, 0, 0, 0); }

        // Sort by gradient magnitude (descending)
        Integer[] indices = new Integer[initialEdgeCount];
        for (int i = 0; i < initialEdgeCount; i++) {
        	indices[i] = i;
        }
        Arrays.sort(indices, (a, b) -> Float.compare(gradPtsList.get(b), gradPtsList.get(a)));
        
        int[] xPtsSorted = new int[initialEdgeCount];
        int[] yPtsSorted = new int[initialEdgeCount];
        for (int i = 0; i < initialEdgeCount; i++) {
            xPtsSorted[i] = xPtsList.get(indices[i]);
            yPtsSorted[i] = yPtsList.get(indices[i]);
        }

        // Apply center exclusion if requested
        if (excludeCenter) {
        	int xCenter = width / 2;
        	int yCenter = height / 2;
            ArrayList<Integer> xFiltered = new ArrayList<>();
            ArrayList<Integer> yFiltered = new ArrayList<>();
            
            for (int i = 0; i < xPtsSorted.length; i++) {
            	double dist = Math.sqrt(
                    Math.pow(xPtsSorted[i] - xCenter, 2) + 
                    Math.pow(yPtsSorted[i] - yCenter, 2)
                );
            	
            	if (dist > centerExclusionRadius) {
            		xFiltered.add(xPtsSorted[i]);
            		yFiltered.add(yPtsSorted[i]);
            	}
            }
        	
            xPtsSorted = toIntArray(xFiltered);
            yPtsSorted = toIntArray(yFiltered);
        }
        
        // Apply neighbor exclusion to find edge centers
        ArrayList<Integer> xCenters = new ArrayList<>();
        ArrayList<Integer> yCenters = new ArrayList<>();
        ArrayList<float[]> centerGrads = new ArrayList<>();
        boolean[] used = new boolean[xPtsSorted.length];

        for (int i = 0; i < xPtsSorted.length; i++) {
            if (used[i]) continue;

            int xPeak = xPtsSorted[i];
            int yPeak = yPtsSorted[i];

            xCenters.add(xPeak);
            yCenters.add(yPeak);

            float sumGradX = 0;
            float sumGradY = 0;
            int count = 0;

            // Compute average gradient in 7x7 window around peak
            for (int dy = -3; dy <= 3; dy++) {
                for (int dx = -3; dx <= 3; dx++) {
                    int nx = xPeak + dx;
                    int ny = yPeak + dy;
                    
                    if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                        int nidx = ny * width + nx;
                        sumGradX += gradX[nidx];
                        sumGradY += gradY[nidx];
                        count++;
                    }
                }
            }

            float avgGradX = sumGradX / count;
            float avgGradY = sumGradY / count;
            centerGrads.add(new float[] {avgGradX, avgGradY});

            // Mark neighbors as used
            for (int j = i + 1; j < xPtsSorted.length; j++) {
                double dist = Math.sqrt(
                    Math.pow(xPtsSorted[j] - xPeak, 2) + 
                    Math.pow(yPtsSorted[j] - yPeak, 2)
                );
                
                if (dist <= neighbourExclusionRadius) { used[j] = true; }
            }
        }

        int filteredEdgeCount = xCenters.size();
        IJ.log("Filter edge count is "+filteredEdgeCount);
        
        if (filteredEdgeCount == 0) {
            return createEmptyResult(initialEdgeCount, 0, 0, 0);
        }

        // Analyze transitions at each edge point
        ArrayList<Integer> xValid = new ArrayList<>();
        ArrayList<Integer> yValid = new ArrayList<>();
        ArrayList<double[]> imgValsAll = new ArrayList<>();
        ArrayList<Double> transitionsValid = new ArrayList<>();
        ArrayList<Double> cosXValid = new ArrayList<>();
        ArrayList<Double> cosYValid = new ArrayList<>();
        int[] xCentersArr = toIntArray(xCenters);
        int[] yCentersArr = toIntArray(yCenters);

        for (int i = 0; i < filteredEdgeCount; i++) {
            int cx = xCentersArr[i];
            int cy = yCentersArr[i];
            float[] gradVec = centerGrads.get(i);
            
            TransitionResult result = estimateEdgeTransitions(
                fp,
                cx,
                cy,
                gradVec,
                lowerBound,
                upperBound,
                subsetSize,
                sectionLength,
                thrMinCriterion,
                thrMaxCriterion,
                minMaxAperture,
                transitionLowLimit,
                transitionHighLimit
            );

            if (result.errorFlag == 0) {
                xValid.add(cx);
                yValid.add(cy);
                imgValsAll.add(result.intensityProfile);
                transitionsValid.add(result.transitionDistance);
                
                double gradMag = Math.sqrt(gradVec[0] * gradVec[0] + gradVec[1] * gradVec[1]);
                cosXValid.add(gradVec[0] / gradMag);  // X gradient
                cosYValid.add(gradVec[1] / gradMag);  // Y gradient
            }
        }
        
        int validCount = xValid.size();
        IJ.log("'Valid' transition count is "+ validCount);
        
        if (validCount == 0) {
            return createEmptyResult(initialEdgeCount, filteredEdgeCount, filteredEdgeCount, 0);
        }

        // Compute statistics
        double[] transitionsArr = toDoubleArray(transitionsValid);
        double meanTrans = Arrays.stream(transitionsArr).average().getAsDouble();
        double stdTrans = std(transitionsArr, meanTrans);
        double minTrans = Arrays.stream(transitionsArr).min().getAsDouble();
        double maxTrans = Arrays.stream(transitionsArr).max().getAsDouble();

        double[][] imgValsArray = imgValsAll.toArray(new double[0][]);
        
        return new EdgeTransitionData(
            toIntArray(xValid),
            toIntArray(yValid),
            toDoubleArray(cosXValid),
            toDoubleArray(cosYValid),
            transitionsArr,
            meanTrans,
            stdTrans,
            minTrans,
            maxTrans,
            initialEdgeCount,
            filteredEdgeCount,
            filteredEdgeCount,
            validCount,
            imgValsArray
        );
    }

    public static float computeGradientThreshold(float[] absGrad, int width, int height, double gradThreshold) {
    	// Create a FloatProcessor from the gradient array to use ThresholdAnalyzer
    	FloatProcessor gradFp = new FloatProcessor(width, height, absGrad);
    	
    	// Use ThresholdAnalyzer with thr_min=0.005, thr_max=gradThreshold
        ThresholdData thresholds = ThresholdAnalyzer.computeThresholds(gradFp, 0.005, gradThreshold, 256);
        
        return (float) thresholds.getMaxThreshold();
    }

    private static TransitionResult estimateEdgeTransitions(
            FloatProcessor fp,
            int centerX,
            int centerY,
            float[] gradVec,
            double lowerBound,
            double upperBound,
            int subsetSize,
            double sectionLength,
            double thrMinCriterion,
            double thrMaxCriterion,
            int minMaxAperture,
            double transitionLowLimit,
            double transitionHighLimit) {
    	
        int width = fp.getWidth();
        int height = fp.getHeight();
        
        // Normalize gradient vector
        double gradMag = Math.sqrt(gradVec[0] * gradVec[0] + gradVec[1] * gradVec[1]);
        if (gradMag < 1e-6) {
            return new TransitionResult(1, 0.0, null);  // Error: no gradient
        }
        
        double cosX = gradVec[0] / gradMag;
        double cosY = gradVec[1] / gradMag;
        
        // Extract subset for min/max criteria
        int halfSubset = subsetSize / 2;
        
        int subsetXi = Math.max(0, centerX - halfSubset);
        int subsetXa = Math.min(width, subsetXi + subsetSize);
        subsetXi = Math.max(0, subsetXa - subsetSize);
        
        int subsetYi = Math.max(0, centerY - halfSubset);
        int subsetYa = Math.min(height, subsetYi + subsetSize);
        subsetYi = Math.max(0, subsetYa - subsetSize);
        
        // Compute subset min/max (via Threshold Analyzer)
        FloatProcessor subsetFp = new FloatProcessor(subsetSize, subsetSize);
        for (int y = 0; y < subsetSize; y++) {
        	for (int x = 0; x < subsetSize; x++) {
        		int srcX = subsetXi + x;
        		int srcY = subsetYi + y;
        		
        		if (srcX >= 0 && srcX < width && srcY >= 0 && srcY < height) {
        			subsetFp.setf(x, y, fp.getf(srcX, srcY));
        		}
        	}
        }
        
        ThresholdData thresholds = ThresholdAnalyzer.computeThresholds(subsetFp, thrMinCriterion, thrMaxCriterion, 256);
        double minCriterion = thresholds.getMinThreshold();
        double maxCriterion = thresholds.getMaxThreshold();
        
        // Create distance array: -section_length // 2 + 1 to section_length // 2 + 1
        int secLen = (int) sectionLength;
        int halfLen = secLen / 2;
        int[] distPix = new int[secLen];
        for (int i = 0; i < secLen; i++) {
        	distPix[i] = -halfLen + 1 + i;
        }
        
        // Extract profile along gradient direction with bilinear interpolation
        double[] imgVal = new double[secLen];
        for (int i = 0; i < secLen; i++) {
        	double xPos = centerX + cosX * distPix[i];
        	double yPos = centerY + cosY * distPix[i];
        	imgVal[i] = bilinearInterpolate(fp, xPos, yPos);
        }
        
        int jc = (secLen - 1) / 2; // Center point index
        
        // Find local min and max positions
        int locMin = 0;
        double minVal = Double.MAX_VALUE;
        for (int i = 0; i <= jc; i++) {
        	if (imgVal[i] < minVal) {
        		minVal = imgVal[i];
        		locMin = i;
        	}
        }
        
        int locMax = jc;
        double maxVal = -Double.MAX_VALUE;
        for (int i = jc; i < secLen; i++) {
        	if (imgVal[i] > maxVal) {
        		maxVal = imgVal[i];
        		locMax = i;
        	}
        }
        
        // Compute mean values with aperture around min and max
        int halfAperture = minMaxAperture / 2;
        // Min aperture indices
        int minStart = locMin - halfAperture + 1;
        int minEnd = locMin + halfAperture + 1;
        if (minStart < 0) {
            int shift = -minStart;
            minStart += shift;
            minEnd += shift;
        }
        
        double minValMean = 0;
        int minCount = 0;
        for (int i = minStart; i < minEnd && i < secLen; i++) {
            minValMean += imgVal[i];
            minCount++;
        }
        minValMean /= minCount;
        
        // Max aperture indices
        int maxStart = locMax - halfAperture + 1;
        int maxEnd = locMax + halfAperture + 1;
        if (maxEnd > secLen) {
            int shift = maxEnd - secLen;
            maxStart -= shift;
            maxEnd -= shift;
        }
        
        double maxValMean = 0;
        int maxCount = 0;
        for (int i = maxStart; i < maxEnd && i < secLen; i++) {
            maxValMean += imgVal[i];
            maxCount++;
        }
        maxValMean /= maxCount;
        
        // Check criteria
        int errorFlag = 0;
        if (minValMean > minCriterion) { errorFlag += 1; }
        if (maxValMean < maxCriterion) { errorFlag += 2; }
        if (errorFlag != 0) { return new TransitionResult(errorFlag, 0.0, null); }
        
        // Calculate transition thresholds
        double transitionMin = minValMean + lowerBound * (maxValMean - minValMean);
        double transitionMax = minValMean + upperBound * (maxValMean - minValMean);
        
        // Scan forward from center to find where profile crosses transitionMax
        int ja = jc;
        while (ja < (secLen - 1) && imgVal[ja] < transitionMax) { ja++; }
        
        if (ja >= secLen - 1) {
        	return new TransitionResult(4, 0.0, null);
        }
        
        // Interpolate exact crossing position
        double xa;
        if (Math.abs(imgVal[ja] - imgVal[ja - 1]) < 1e-9) {
            xa = distPix[ja];
        } else {
            xa = distPix[ja - 1] + (transitionMax - imgVal[ja - 1]) * 
                 (distPix[ja] - distPix[ja - 1]) / (imgVal[ja] - imgVal[ja - 1]);
        }
        
        // Scan backward from center to find where profile crosses transitionMin
        int ji = jc;
        while (ji > 0 && imgVal[ji] > transitionMin) { ji--; }
        
        if (ji <= 0) {
            return new TransitionResult(4, 0.0, null);
        }
        
        // Interpolate exact crossing position
        double xi;
        if (Math.abs(imgVal[ji] - imgVal[ji + 1]) < 1e-9) {
            xi = distPix[ji];
        } else {
            xi = distPix[ji + 1] + (transitionMin - imgVal[ji + 1]) * 
                 (distPix[ji] - distPix[ji + 1]) / (imgVal[ji] - imgVal[ji + 1]);
        }
        
        double transitionDistance = Math.abs(xa - xi);
        
        // Check limits
        if (transitionDistance < transitionLowLimit) { errorFlag += 4; }
        if (transitionDistance > transitionHighLimit) { errorFlag += 8; }
        if (Double.isNaN(transitionDistance)) { errorFlag += 16; }
        
        return new TransitionResult(errorFlag, transitionDistance, imgVal);
    }
    
    /**
     * Bilinear interpolation for subpixel-accurate intensity values
     */
    private static double bilinearInterpolate(FloatProcessor fp, double x, double y) {
    	int width = fp.getWidth();
    	int height = fp.getHeight();
    	
    	// Clamp to image bounds
    	if (x < 0) { x = 0; }
    	if (x >= width - 1) { x = width - 1 - 1e-6; }
    	if (y < 0) { y = 0; }
    	if (y >= height - 1) { y = height - 1 - 1e-6; }
	   
    	int x0 = (int) x;
    	int y0 = (int) y;
    	int x1 = x0 + 1;
    	int y1 = y0 + 1;
   
    	double dx = x - x0;
    	double dy = y - y0;
   
    	double v00 = fp.getf(x0, y0);
    	double v10 = fp.getf(x1, y0);
    	double v01 = fp.getf(x0, y1);
    	double v11 = fp.getf(x1, y1);
   
    	double v0 = v00 * (1 - dx) + v10 * dx;
    	double v1 = v01 * (1 - dx) + v11 * dx;
   
    	return v0 * (1 - dy) + v1 * dy;
    }

    /**
     * Creates empty result for cases with no valid data.
     */
    private static EdgeTransitionData createEmptyResult(int initialCount, int filteredCount, int analyzedCount, int validCount) {
        return new EdgeTransitionData(
            new int[0], new int[0],
            new double[0], new double[0],
            new double[0],
            0.0, 0.0, 0.0, 0.0,
            initialCount, filteredCount, analyzedCount, validCount,
            null
        );
    }

    private static double std(double[] arr, double mean) {
        if (arr.length == 0) return 0.0;
        double sumSq = 0;
        for (double v : arr) {
            double diff = v - mean;
            sumSq += diff * diff;
        }
        return Math.sqrt(sumSq / (arr.length - 1));
    }

    private static int[] toIntArray(ArrayList<Integer> list) {
        int[] arr = new int[list.size()];
        for (int i = 0; i < list.size(); i++) {
            arr[i] = list.get(i);
        }
        return arr;
    }
    
    private static double[] toDoubleArray(ArrayList<Double> list) {
        double[] arr = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            arr[i] = list.get(i);
        }
        return arr;
    }

    /**
     * Simple container for transition estimation result.
     */
    private static class TransitionResult {
        final int errorFlag;
        final double transitionDistance;
        final double[] intensityProfile;
        
        TransitionResult(int errorFlag, double transitionDistance, double[] intensityProfile) {
            this.errorFlag = errorFlag;
            this.transitionDistance = transitionDistance;
            this.intensityProfile = intensityProfile;
        }
    }
}
