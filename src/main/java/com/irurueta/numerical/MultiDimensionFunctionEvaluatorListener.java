/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Interface to define how multi dimension functions can be evaluated.
 * This interface is used in several algorithms to provide methods to evaluate
 * functions and retrieve their minima/maxima, etc.
 */
public interface MultiDimensionFunctionEvaluatorListener {
    
    /**
     * Evaluates a multi dimension function such as f([x1, x2, ..., xn]) at 
     * provided multidimensional point and returns the result as a scalar value.
     * @param point Multidimensional point where function will be evaluated. 
     * This must be an array of length equal to the dimensionality of the 
     * function
     * @return Value returned by the function.
     * @throws Throwable Raised if something failed during the evaluation.
     */    
    public double evaluate(double[] point) throws Throwable;    
}
