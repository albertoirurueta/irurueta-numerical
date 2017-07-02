/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialEvaluationType
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

/**
 * Determines different types of polynomial evaluations that can be used
 * to estimate a polynomial.
 */
public enum PolynomialEvaluationType {
    /**
     * A direct evaluation of a polynomial.
     */
    DIRECT_EVALUATION,
    
    /**
     * Evaluation of the nth-derivative of a polynomial.
     */
    DERIVATIVE_EVALUATION,
    
    /**
     * Evaluation of the nth-integral of a polynomial (assuming a given 
     * constant value)
     */
    INTEGRAL_EVALUATION,
    
    /**
     * Interval nth-integration of a polynomial.
     */
    INTEGRAL_INTERVAL
}
