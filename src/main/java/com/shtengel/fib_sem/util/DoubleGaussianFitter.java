package com.shtengel.fib_sem.util;

import org.apache.commons.math3.fitting.leastsquares.*;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.Pair;

import ij.IJ;

/**
 * Fits a sum of two Gaussians (no offset) to a 1-D histogram / PDF.
 *
 * <h3>Fitting strategy</h3>
 * <ol>
 *   <li>Find the dominant peak of the PDF and fit a single Gaussian.</li>
 *   <li>Subtract that Gaussian to obtain a residual; fit a second Gaussian to it.</li>
 *   <li>Use both single-Gaussian results as the initial guess for a full
 *       6-parameter Levenberg–Marquardt refinement (via Commons Math).</li>
 * </ol>
 */
public class DoubleGaussianFitter {

    private static final int MAX_EVALUATIONS = 1000;
    private static final int MAX_ITERATIONS  = 500;

    private DoubleGaussianFitter() {}

    /**
     * Holds the fitted parameters {@code [A1, μ1, σ1, A2, μ2, σ2]}
     * with μ1 &lt; μ2 guaranteed, plus the coefficient of determination (R²).
     */
    public static class FitResult {

        private final double[] params;
        private final double rSquared;

        public FitResult(double[] params, double rSquared) {
            this.params = params;
            this.rSquared = rSquared;
        }

        public double[] getParams() { return params; }
        public double rSquared()    { return rSquared; }

        /** Evaluate the fitted double-Gaussian at {@code x}. */
        public double evaluate(double x) {
            return doubleGaussian(x, params);
        }

        public double getA1()     { return params[0]; }
        public double getMu1()    { return params[1]; }
        public double getSigma1() { return params[2]; }
        public double getA2()     { return params[3]; }
        public double getMu2()    { return params[4]; }
        public double getSigma2() { return params[5]; }
    }

    /**
     * Fits a double-Gaussian model to the supplied PDF histogram.
     *
     * @param binCenters  x-coordinates (intensity bin centres)
     * @param pdf         y-coordinates (probability density values)
     * @return the fit result, with μ1 &lt; μ2
     */
    public static FitResult fit(double[] binCenters, double[] pdf) {
        int n = binCenters.length;

        // Fit dominant peak
        int peakIdx = argmax(pdf, 0, n);
        double[] g1 = fitSingleGaussian(binCenters, pdf,
            pdf[peakIdx], binCenters[peakIdx], estimateSigma(binCenters, pdf, peakIdx));

        // Subtract to get fit residual
        double[] residual = new double[n];
        for (int i = 0; i < n; i++) {
            residual[i] = Math.max(0, pdf[i] - singleGaussian(binCenters[i], g1));
        }

        // Search for residual peak OUTSIDE the dominant peak's core (±1σ).
        // Without this, the residual peak lands in the overlap region
        // between two peaks rather than at the true secondary peak centre.
        double dominantMu    = g1[1];
        double dominantSigma = Math.abs(g1[2]);
        int resPeakIdx = -1;
        double bestResidual  = -1;
        for (int i = 0; i < n; i++) {
            if (Math.abs(binCenters[i] - dominantMu) > dominantSigma
                    && residual[i] > bestResidual) {
                bestResidual = residual[i];
                resPeakIdx = i;
            }
        }
        // Fallback: if all bins are within 1σ (very narrow distribution), use global max
        if (resPeakIdx < 0) {
            resPeakIdx = argmax(residual, 0, n);
        }
        double[] g2 = fitSingleGaussian(binCenters, residual,
            residual[resPeakIdx], binCenters[resPeakIdx], estimateSigma(binCenters, residual, resPeakIdx));

        // Joint 6-param refinement
        double[] init = {
            g1[0], g1[1], Math.abs(g1[2]),
            g2[0], g2[1], Math.abs(g2[2])
        };

        double[] fitted = fitDoubleGaussianLM(binCenters, pdf, init);

        // Ensure σ > 0
        fitted[2] = Math.abs(fitted[2]);
        fitted[5] = Math.abs(fitted[5]);

        // Ensure μ1 < μ2
        if (fitted[1] > fitted[4]) {
            swap(fitted, 0, 3);
            swap(fitted, 1, 4);
            swap(fitted, 2, 5);
        }

        double r2 = computeRSquared(binCenters, pdf, fitted);

        IJ.log(String.format(
            "Double-Gaussian fit:  μ1=%.2f  σ1=%.2f  A1=%.4f  |  μ2=%.2f  σ2=%.2f  A2=%.4f  |  R²=%.5f",
            fitted[1], fitted[2], fitted[0], fitted[4], fitted[5], fitted[3], r2));

        return new FitResult(fitted, r2);
    }

    /**
     * Fits a double-Gaussian model while keeping both means fixed. The supplied
     * sigmas are initial guesses only; both sigmas and amplitudes are optimized.
     * The returned parameter order matches {@link FitResult}.
     */
    public static FitResult fitWithFixedMeansAndSigmas(
            double[] x, double[] y,
            double mu1, double sigma1, double mu2, double sigma2) {
        if (x == null || y == null || x.length != y.length || x.length < 3) {
            throw new IllegalArgumentException("x and y must have the same length (at least 3)");
        }
        if (!(sigma1 > 0.0) || !(sigma2 > 0.0)
                || !Double.isFinite(sigma1) || !Double.isFinite(sigma2)) {
            throw new IllegalArgumentException("Initial Gaussian sigmas must be finite and positive");
        }
        if (mu1 > mu2) {
            double t = mu1; mu1 = mu2; mu2 = t;
            t = sigma1; sigma1 = sigma2; sigma2 = t;
        }

        // Estimate initial amplitudes by linear least squares at given means and initial sigmas.
        double c11 = 0.0, c22 = 0.0, c12 = 0.0;
        double b1 = 0.0, b2 = 0.0;
        double[] g1 = new double[x.length];
        double[] g2 = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double d1 = (x[i] - mu1) / sigma1;
            double d2 = (x[i] - mu2) / sigma2;
            g1[i] = Math.exp(-0.5 * d1 * d1);
            g2[i] = Math.exp(-0.5 * d2 * d2);
            c11 += g1[i] * g1[i];
            c22 += g2[i] * g2[i];
            c12 += g1[i] * g2[i];
            b1 += g1[i] * y[i];
            b2 += g2[i] * y[i];
        }
        double det = c11 * c22 - c12 * c12;
        double a1 = 0.0, a2 = 0.0;
        if (det > 1e-30) {
            a1 = Math.max(0.0, (b1 * c22 - b2 * c12) / det);
            a2 = Math.max(0.0, (b2 * c11 - b1 * c12) / det);
        } else {
            a1 = c11 > 0.0 ? Math.max(0.0, b1 / c11) : 0.0;
            a2 = c22 > 0.0 ? Math.max(0.0, b2 / c22) : 0.0;
        }
        // Keep positive start so sigma derivatives are available when
		// initial linear estimate assigns a component zero amplitude.
        double maxY = 0.0;
        for (double v : y) maxY = Math.max(maxY, Math.abs(v));
        double amplitudeFloor = Math.max(1e-12, maxY * 1e-9);

        // Optimize [log(A1), log(sigma1), log(A2), log(sigma2)]. 
		// The log parameterization keeps both amplitudes and sigmas positive.
        final double fixedMu1 = mu1;
        final double fixedMu2 = mu2;
        MultivariateJacobianFunction model = point -> {
            double[] p = point.toArray();
            double amp1 = Math.exp(p[0]);
            double sig1 = Math.exp(p[1]);
            double amp2 = Math.exp(p[2]);
            double sig2 = Math.exp(p[3]);
            double[] values = new double[x.length];
            double[][] jacobian = new double[x.length][4];
            for (int i = 0; i < x.length; i++) {
                double d1 = (x[i] - fixedMu1) / sig1;
                double d2 = (x[i] - fixedMu2) / sig2;
                double e1 = Math.exp(-0.5 * d1 * d1);
                double e2 = Math.exp(-0.5 * d2 * d2);
                values[i] = amp1 * e1 + amp2 * e2;
                jacobian[i][0] = amp1 * e1;
                jacobian[i][1] = amp1 * e1 * d1 * d1;
                jacobian[i][2] = amp2 * e2;
                jacobian[i][3] = amp2 * e2 * d2 * d2;
            }
            return new Pair<>(new ArrayRealVector(values), new Array2DRowRealMatrix(jacobian));
        };

        LeastSquaresProblem problem = new LeastSquaresBuilder()
            .start(new double[] {
                Math.log(Math.max(amplitudeFloor, a1)), Math.log(sigma1),
                Math.log(Math.max(amplitudeFloor, a2)), Math.log(sigma2)
            })
            .model(model)
            .target(y)
            .maxEvaluations(MAX_EVALUATIONS)
            .maxIterations(MAX_ITERATIONS)
            .lazyEvaluation(false)
            .build();

        LeastSquaresOptimizer.Optimum optimum =
            new LevenbergMarquardtOptimizer().optimize(problem);
        double[] fitted = optimum.getPoint().toArray();
        double[] params = {
            Math.exp(fitted[0]), fixedMu1, Math.exp(fitted[1]),
            Math.exp(fitted[2]), fixedMu2, Math.exp(fitted[3])
        };
        double r2 = computeRSquared(x, y, params);
        IJ.log(String.format(
            "Fixed-mean double-Gaussian fit: μ1=%.2f σ1=%.2f A1=%.4f | μ2=%.2f σ2=%.2f A2=%.4f | R²=%.5f",
            params[1], params[2], params[0], params[4], params[5], params[3], r2));
        return new FitResult(params, r2);
    }

    /**
     * 6-parameter double-Gaussian fit using Commons Math LevenbergMarquardtOptimizer.
     */
    private static double[] fitDoubleGaussianLM(double[] x, double[] y, double[] start) {
        MultivariateJacobianFunction model = point -> {
            double[] p = point.toArray();
            int n = x.length;
            double[] values = new double[n];
            double[][] jacobian = new double[n][6];

            for (int i = 0; i < n; i++) {
                values[i] = doubleGaussian(x[i], p);
                fillJacobianRow(x[i], p, jacobian[i]);
            }

            return new Pair<>(
                new ArrayRealVector(values),
                new Array2DRowRealMatrix(jacobian)
            );
        };

        LeastSquaresProblem problem = new LeastSquaresBuilder()
            .start(start)
            .model(model)
            .target(y)
            .maxEvaluations(MAX_EVALUATIONS)
            .maxIterations(MAX_ITERATIONS)
            .lazyEvaluation(false)
            .build();

        LeastSquaresOptimizer.Optimum optimum =
            new LevenbergMarquardtOptimizer().optimize(problem);

        return optimum.getPoint().toArray();
    }

    /**
     * 3-parameter single-Gaussian fit (used for initialisation).
     */
    private static double[] fitSingleGaussian(double[] x, double[] y,
                                               double a0, double mu0, double sigma0) {
        double[] start = { a0, mu0, Math.abs(sigma0) + 1e-6 };

        MultivariateJacobianFunction model = point -> {
            double[] p = point.toArray();
            int n = x.length;
            double[] values = new double[n];
            double[][] jacobian = new double[n][3];

            for (int i = 0; i < n; i++) {
                double a = p[0], mu = p[1], sigma = p[2];
                double dx = x[i] - mu;
                double s2 = sigma * sigma;
                double expTerm = Math.exp(-0.5 * dx * dx / s2);

                values[i] = a * expTerm;
                jacobian[i][0] = expTerm;                                  // ∂f/∂A
                jacobian[i][1] = a * (dx / s2) * expTerm;                 // ∂f/∂μ
                jacobian[i][2] = a * (dx * dx / (s2 * sigma)) * expTerm;  // ∂f/∂σ
            }

            return new Pair<>(
                new ArrayRealVector(values),
                new Array2DRowRealMatrix(jacobian)
            );
        };

        LeastSquaresProblem problem = new LeastSquaresBuilder()
            .start(start)
            .model(model)
            .target(y)
            .maxEvaluations(MAX_EVALUATIONS)
            .maxIterations(MAX_ITERATIONS)
            .lazyEvaluation(false)
            .build();

        LeastSquaresOptimizer.Optimum optimum =
            new LevenbergMarquardtOptimizer().optimize(problem);

        return optimum.getPoint().toArray();
    }

    static double doubleGaussian(double x, double[] p) {
        return singleGaussian(x, p[0], p[1], p[2])
             + singleGaussian(x, p[3], p[4], p[5]);
    }

    private static double singleGaussian(double x, double a, double mu, double sigma) {
        double dx = x - mu;
        return a * Math.exp(-0.5 * dx * dx / (sigma * sigma));
    }

    private static double singleGaussian(double x, double[] p) {
        return singleGaussian(x, p[0], p[1], p[2]);
    }

    private static void fillJacobianRow(double x, double[] p, double[] out) {
        for (int g = 0; g < 2; g++) {
            int off = g * 3;
            double a     = p[off];
            double mu    = p[off + 1];
            double sigma = p[off + 2];
            double dx    = x - mu;
            double s2    = sigma * sigma;
            double expTerm = Math.exp(-0.5 * dx * dx / s2);

            out[off]     = expTerm;                                  // ∂f/∂A
            out[off + 1] = a * (dx / s2) * expTerm;                 // ∂f/∂μ
            out[off + 2] = a * (dx * dx / (s2 * sigma)) * expTerm;  // ∂f/∂σ
        }
    }

    private static double computeRSquared(double[] x, double[] y, double[] p) {
        double mean = 0;
        for (double v : y) mean += v;
        mean /= y.length;

        double ssTot = 0, ssRes = 0;
        for (int i = 0; i < x.length; i++) {
            ssTot += (y[i] - mean) * (y[i] - mean);
            double r = y[i] - doubleGaussian(x[i], p);
            ssRes += r * r;
        }
        return 1.0 - ssRes / (ssTot + 1e-30);
    }

    /** FWHM-based σ estimate around a peak. */
    private static double estimateSigma(double[] x, double[] y, int peakIdx) {
        double halfMax = y[peakIdx] / 2.0;

        int left = peakIdx;
        while (left > 0 && y[left] > halfMax) left--;
        int right = peakIdx;
        while (right < x.length - 1 && y[right] > halfMax) right++;

        double fwhm = x[right] - x[left];
        return Math.max(fwhm / 2.355, (x[1] - x[0]) * 2);
    }

    private static int argmax(double[] data, int from, int to) {
        int best = from;
        for (int i = from + 1; i < to; i++) {
            if (data[i] > data[best]) best = i;
        }
        return best;
    }

    private static void swap(double[] arr, int i, int j) {
        double tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}