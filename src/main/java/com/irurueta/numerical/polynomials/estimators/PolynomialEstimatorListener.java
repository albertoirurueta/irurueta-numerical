/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

/**
 * Listener to be notified when estimation starts, finishes or any progress
 * changes.
 */
public interface PolynomialEstimatorListener {
    
    /**
     * Called when an estimator starts the polynomial estimation process.
     * @param estimator reference to a polynomial estimator.
     */
    void onEstimateStart(PolynomialEstimator estimator);
    
    /**
     * Called when an estimator ends the polynomial estimation process.
     * @param estimator reference to a polynomial estimator.
     */
    void onEstimateEnd(PolynomialEstimator estimator);
}
