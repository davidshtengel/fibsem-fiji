package com.shtengel.fib_sem.data;

/**
 * Container for edge transition analysis results.
 * 
 * <p> This class stores edge point locations, gradient information, 
 * transition distances, and statistical summaries from edge transition analysis. </p>
 * 
 * <p><b>NOTE:</b> Arrays with "Selected" suffix contain only valid transitions (error flag = 0)</p>
 */
public class EdgeTransitionData {
    private final int[] xSelected;
    private final int[] ySelected;
    private final double[] cosXSelected;
    private final double[] cosYSelected;

    private final double[] transitionDistances;
    private final double[][] imageValuesAll;
    
    private final double meanTransition;
    private final double stdTransition;
    private final double minTransition;
    private final double maxTransition;

    private final int initialEdgeCount;
    private final int filteredEdgeCount;
    private final int analyzedCount;
    private final int validCount;

    public EdgeTransitionData(
            int[] xSelected,
            int[] ySelected,
            double[] cosXSelected,
            double[] cosYSelected,
            double[] transitionDistances,
            double meanTransition,
            double stdTransition,
            double minTransition,
            double maxTransition,
            int initialEdgeCount,
            int filteredEdgeCount,
            int analyzedCount,
            int validCount,
            double[][] imageValuesAll) {
        this.xSelected = xSelected;
        this.ySelected = ySelected;
        this.cosXSelected = cosXSelected;
        this.cosYSelected = cosYSelected;
        this.transitionDistances = transitionDistances;
        this.meanTransition = meanTransition;
		this.stdTransition = stdTransition;
        this.minTransition = minTransition;
        this.maxTransition = maxTransition;
        this.initialEdgeCount = initialEdgeCount;
        this.filteredEdgeCount = filteredEdgeCount;
        this.analyzedCount = analyzedCount;
        this.validCount = validCount;
        this.imageValuesAll = imageValuesAll;
    }

    // Getters
    public int[] getXSelected() { return xSelected; }
    public int[] getYSelected() { return ySelected; }
    public double[] getCosXSelected() { return cosXSelected; }
    public double[] getCosYSelected() { return cosYSelected; }
    public double[] getTransitionDistances() { return transitionDistances; }
    public double getMeanTransition() { return meanTransition; }
    public double getStdTransition() { return stdTransition; }
    public double getMinTransition() { return minTransition; }
    public double getMaxTransition() { return maxTransition; }
    public int getInitialEdgeCount() { return initialEdgeCount; }
    public int getFilteredEdgeCount() { return filteredEdgeCount; }
    public int getAnalyzedCount() { return analyzedCount; }
    public int getValidCount() { return validCount; }
    public double[][] getImageValuesAll() { return imageValuesAll; }

    @Override
    public String toString() {
        return String.format(
            "EdgeTransitionData[valid=%d/%d, mean=%.3f\u00B1%.3f px, range=[%.3f, %.3f]]",
            validCount, analyzedCount, meanTransition, stdTransition, minTransition, maxTransition
        );
    }
}
