package com.shtengel.fib_sem.util;

import java.util.Arrays;
import org.apache.commons.math3.fitting.leastsquares.GaussNewtonOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ij.IJ;

/**
 * Fits a centered ellipse to angular transition distance data using
 * nonlinear least-squares optimization.
 *
 * <p>Model: r(θ) = 1 / √((cos(θ-φ)/a)² + (sin(θ-φ)/b)²),
 * where a, b are semi-axes and φ is the rotation angle.</p>
 */
public class EllipseFitter {
	
	private static final int N_BINS = 36;
	
	/**
	 * Fits a centered ellipse to transition distance data parameterized by direction.
	 *
	 * @param cosX   x-component of gradient direction cosines
	 * @param cosY   y-component of gradient direction cosines
	 * @param transitions  measured transition distances (must be positive)
	 * @return fitted parameters [a, b, phi] with a >= b, or null if fitting fails
	 */
	public static double[] fitCenteredEllipse(double[] cosX, double[] cosY, double[] transitions) {
		if (!validateInput(cosX, cosY, transitions)) {
			return null;
		}
		
		final double[] theta = toPolarCoordinates(cosX, cosY);
		
		// Initial parameter estimates
	    double meanTrans = Arrays.stream(transitions).average().getAsDouble();
	    double variance = 0.0;
	    for (double t : transitions) {
	        variance += (t - meanTrans) * (t - meanTrans);
	    }
	    double stdTrans = Math.sqrt(variance / transitions.length);
	    
	    // Validate initial parameters
	    if (meanTrans <= 0 || Double.isNaN(meanTrans) || Double.isNaN(stdTrans)) {
	        IJ.log(String.format("fitCenteredEllipse ERROR: invalid initial parameters - mean:%f, std:%f", 
	            meanTrans, stdTrans));
	        return null;
	    }
	    
	    IJ.log(String.format("Initial estimates: mean=%.3f, std=%.3f", meanTrans, stdTrans));
	    
	    final double[] observed = transitions.clone();
	    
	    // Model: r(theta) = 1 / sqrt((cos(theta-phi)/a)^2 + (sin(theta-phi)/b)^2)
	    // 	a, b = semi-major, semi-minor axes
	    //	phi = rotation angle 
	    MultivariateJacobianFunction model = createEllipseModel(theta);
	    
	    // Estimate initial guess from angular binning
	    double[] initialGuess = estimateInitialGuess(theta, transitions, meanTrans);
	    IJ.log(String.format("Initial guess: a=%.3f, b=%.3f, phi=%.3f rad (%.1f\u00B0)", 
	        initialGuess[0], initialGuess[1], initialGuess[2], Math.toDegrees(initialGuess[2])));
	    
	 // Attempt fit with Levenberg-Marquardt, fall back to Gauss-Newton
	    double[] fittedParams = fitWithLM(model, observed, initialGuess);
	    
	    if (fittedParams == null) {
	        fittedParams = fitWithGN(model, observed, initialGuess);
	    }
	    
	    if (fittedParams == null) {
	        return null;
	    }
	    
	    // Normalize: ensure a >= b, take absolute values, and wrap phi
	    fittedParams = normalizeParameters(fittedParams);
	    
	    double ellipticity = computeEllipticity(fittedParams[0], fittedParams[1]);
	    
	    IJ.log(String.format("Ellipse fit parameters: a=%.3f, b=%.3f, phi=%.3f rad (%.1f\u00B0), ellipticity=%.3f", 
	        fittedParams[0], fittedParams[1], fittedParams[2], 
	        Math.toDegrees(fittedParams[2]), ellipticity));
	    
	    return fittedParams;
	}
	
	/**
	 * Computes ellipticity from semi-major and semi-minor axes.
	 *
	 * @param a semi-major axis (must be >= b)
	 * @param b semi-minor axis
	 * @return ellipticity in [0, 1), where 0 = circle
	 */
	public static double computeEllipticity(double a, double b) {
	    if (a <= 0) return 0.0;
	    return 1.0 - (Math.min(a, b) / Math.max(a, b));
	}
	
	/**
	 * Creates the ellipse model function with analytical Jacobian.
	 */
	private static MultivariateJacobianFunction createEllipseModel(final double[] theta) {
	    return new MultivariateJacobianFunction() {
	        @Override
	        public org.apache.commons.math3.util.Pair<RealVector, RealMatrix> value(RealVector point) {
	            double a = point.getEntry(0);
	            double b = point.getEntry(1);
	            double phi = point.getEntry(2);
	            
	            double[] values = new double[theta.length];
	            double[][] jacobian = new double[theta.length][3];
	            
	            for (int i = 0; i < theta.length; i++) {
	                double thetaRot = theta[i] - phi;
	                double cosT = Math.cos(thetaRot);
	                double sinT = Math.sin(thetaRot);
	                
	                double term1 = cosT / a;
	                double term2 = sinT / b;
	                double denom = Math.sqrt(term1 * term1 + term2 * term2);
	                
	                values[i] = 1.0 / denom;
	                
	                // Partial derivatives: dr/da, dr/db, dr/dphi
	                double denom3 = denom * denom * denom;
	                jacobian[i][0] = (cosT * cosT) / (a * a * a * denom3);
	                jacobian[i][1] = (sinT * sinT) / (b * b * b * denom3);
	                jacobian[i][2] = (sinT * cosT * (1.0 / (b * b) - 1.0 / (a * a))) / denom3;
	            }
	            
	            return new org.apache.commons.math3.util.Pair<>(
	                new ArrayRealVector(values),
	                new Array2DRowRealMatrix(jacobian)
	            );
	        }
	    };
	}
	
	private static double[] estimateInitialGuess(double[] theta, double[] transitions, double meanTrans) {
		double[] binSums = new double[N_BINS];
	    int[] binCounts = new int[N_BINS];

	    for (int i = 0; i < transitions.length; i++) {
	        // Map theta from [-pi, pi] to [0, 2*pi]
	        double angle = theta[i];
	        if (angle < 0) angle += 2 * Math.PI;
	        
	        int bin = (int) ((angle / (2 * Math.PI)) * N_BINS);
	        if (bin >= N_BINS) bin = N_BINS - 1;
	        
	        binSums[bin] += transitions[i];
	        binCounts[bin]++;
	    }

	    // Calculate mean for each bin (use overall mean for empty bins)
	    double[] binMeans = new double[N_BINS];
	    for (int i = 0; i < N_BINS; i++) {
	        binMeans[i] = (binCounts[i] > 0) ? binSums[i] / binCounts[i] : meanTrans;
	    }

	    // Find the bin with maximum mean (approximate major axis direction)
	    double maxBinMean = Double.NEGATIVE_INFINITY;
	    int maxBin = 0;

	    for (int i = 0; i < N_BINS; i++) {
	        if (binCounts[i] > 0 && binMeans[i] > maxBinMean) {
	            maxBinMean = binMeans[i];
	            maxBin = i;
	        }
	    }

	    // Initial rotation angle from the max bin
	    double phiInitial = (maxBin * 2 * Math.PI / N_BINS) - Math.PI;

	    // Use 80th and 20th percentiles of bin means for robust axis estimates
	    double[] sortedBinMeans = binMeans.clone();
	    Arrays.sort(sortedBinMeans);
	    double aInitial = sortedBinMeans[(int) (N_BINS * 0.8)];
	    double bInitial = sortedBinMeans[(int) (N_BINS * 0.2)];

	    // Ensure a >= b (semi-major >= semi-minor)
	    if (aInitial < bInitial) {
	        double temp = aInitial;
	        aInitial = bInitial;
	        bInitial = temp;
	        phiInitial += Math.PI / 2;
	    }

	    // Normalize phi to [-pi, pi]
	    phiInitial = normalizePhi(phiInitial);
		
	    return new double[] {aInitial, bInitial, phiInitial};
	}
	
	/**
	 * Attempts fit using Levenberg-Marquardt optimizer
	 * 
	 * @returns fitted [a, b, phi] (null on failure)
	 */
	private static double[] fitWithLM(MultivariateJacobianFunction model, double[] observed, double[] initialGuess) {
		try {
	        LeastSquaresProblem problem = new LeastSquaresBuilder()
	            .model(model)
	            .target(observed)
	            .start(initialGuess)
	            .maxEvaluations(2000)
	            .maxIterations(2000)
	            .build();
	        
	        LevenbergMarquardtOptimizer lmOptimizer = new LevenbergMarquardtOptimizer();
	        LeastSquaresOptimizer.Optimum optimum = lmOptimizer.optimize(problem);
	        
	        return optimum.getPoint().toArray();
	        
	    } catch (Exception e) {
	        IJ.log("Levenberg-Marquardt fitting failed: " + e.getMessage());
	        return null;
	    }
	}
	
	/**
	 * Attempts fit using Gauss-Newton optimizer (fallback).
	 *
	 * @return fitted [a, b, phi] or null on failure
	 */
	private static double[] fitWithGN(MultivariateJacobianFunction model, double[] observed, double[] initialGuess) {
	    try {
	        IJ.log("Attempting fallback with Gauss-Newton optimizer...");
	        
	        LeastSquaresProblem problem = new LeastSquaresBuilder()
	            .model(model)
	            .target(observed)
	            .start(initialGuess)
	            .maxEvaluations(1000)
	            .maxIterations(1000)
	            .build();
	        
	        GaussNewtonOptimizer gnOptimizer = new GaussNewtonOptimizer();
	        LeastSquaresOptimizer.Optimum optimum = gnOptimizer.optimize(problem);
	        
	        double[] fittedParams = optimum.getPoint().toArray();
	        IJ.log(String.format("Fallback fit successful: a=%.4f, b=%.4f, phi=%.4f", 
	            fittedParams[0], fittedParams[1], fittedParams[2]));
	        
	        return fittedParams;
	        
	    } catch (Exception e) {
	        IJ.log("Fallback fitting also failed: " + e.getMessage());
	        return null;
	    }
	}
	
	/**
	 * Normalizes fitted parameters: takes absolute values of axes,
	 * ensures a >= b, and wraps phi to [-pi, pi].
	 */
	private static double[] normalizeParameters(double[] params) {
	    double a = Math.abs(params[0]);
	    double b = Math.abs(params[1]);
	    double phi = params[2];
	    
	    // Ensure a >= b; if swapped, rotate by 90 degrees
	    if (a < b) {
	        double temp = a;
	        a = b;
	        b = temp;
	        phi += Math.PI / 2;
	    }
	    
	    phi = normalizePhi(phi);
	    
	    return new double[] { a, b, phi };
	}
	
	/**
	 * Wraps an angle to the range [-pi, pi].
	 */
	private static double normalizePhi(double phi) {
		while (phi > Math.PI) phi -= 2 * Math.PI;
		while (phi < -Math.PI) phi += 2 * Math.PI;
		return phi;
	}
	
	private static boolean validateInput(double[] cosX, double[] cosY, double[] transitions) {
		// Comprehensive input validation
	    if (cosX == null || cosY == null || transitions == null) {
	    	IJ.log("fitCenteredEllipse ERROR: null input array detected");
	    	return false;
	    }
	    if (cosX.length == 0 || cosY.length == 0 || transitions.length == 0) {
	    	IJ.log("fitCenteredEllipse ERROR: empty input arrays");
	    	return false;
	    }
	    if (cosX.length != cosY.length || cosX.length != transitions.length) {
	    	IJ.log(String.format("fitCenteredEllipse ERROR: array length mismatch - cosX:%d, cosY:%d, transitions:%d", cosX.length, cosY.length, transitions.length));
	    	return false;
	    }
	    
	    // Check for invalid values
	    for (int i = 0; i < transitions.length; i++) {
	        if (Double.isNaN(transitions[i]) || Double.isInfinite(transitions[i]) || transitions[i] <= 0) {
	            IJ.log(String.format("fitCenteredEllipse ERROR: invalid transition value at index %d: %f", i, transitions[i]));
	            return false;
	        }
	        if (Double.isNaN(cosX[i]) || Double.isInfinite(cosX[i]) || Double.isNaN(cosY[i]) || Double.isInfinite(cosY[i])) {
	            IJ.log(String.format("fitCenteredEllipse ERROR: invalid cosX/cosY value at index %d", i));
	            return false;
	        }
	    }
	    
	    return true;
	}
	
	private static double[] toPolarCoordinates(double[] cosX, double[] cosY) {
		final double[] theta = new double[cosX.length];
	    for (int i = 0; i < cosX.length; i++) {
	        theta[i] = Math.atan2(cosY[i], cosX[i]);
	    }
	    
	    return theta;
	}
	
}