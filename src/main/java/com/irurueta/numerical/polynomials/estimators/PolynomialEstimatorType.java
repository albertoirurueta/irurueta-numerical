/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialEstimatorType
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

/**
 * Polynomial estimator types.
 */
public enum PolynomialEstimatorType {
    /**
     * Polynomial estimator using LMSE (Least Mean Square Error) solutions.
     */
    LMSE_POLYNOMIAL_ESTIMATOR,
    
    /**
     * Polynomial estimator using weighted evaluations.
     */
    WEIGHTED_POLyNOMIAL_ESTIMATOR
}
