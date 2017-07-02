/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialRobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 13, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

/**
 * Listener to be notified of events such as when estimation starts, ends or
 * when progress changes.
 */
public interface PolynomialRobustEstimatorListener {
    
    /**
     * Called when an estimator starts the polynomial estimation process.
     * @param estimator reference to a polynomial robust estimator.
     */
    void onEstimateStart(PolynomialRobustEstimator estimator);
    
    /**
     * Called when an estimator ends the polynomial estimation process.
     * @param estimator reference to a polynomial robust estimator.
     */
    void onEstimateEnd(PolynomialRobustEstimator estimator);
    
    /**
     * Called when estimator iterates to refine a possible solution.
     * @param estimator reference to a polynomial robust estimator.
     * @param iteration current iteration.
     */
    void onEstimateNextIteration(PolynomialRobustEstimator estimator, 
            int iteration);
    
    /**
     * Called when estimation progress changes significantly.
     * @param estimator reference to a polynomial robust estimator.
     * @param progress progress of estimation expressed as a value between 0.0
     * and 1.0.
     */
    void onEstimateProgressChange(PolynomialRobustEstimator estimator,
            float progress);
}
