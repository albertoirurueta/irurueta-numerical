/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Interface to define how single dimension functions can be evaluated.
 * This interface is used in several algorithms to provide methods to evaluate
 * functions and retrieve their minima/maxima, etc.
 */
public interface SingleDimensionFunctionEvaluatorListener {
    
    /**
     * Evaluates a single dimension function such as f(x) at provided point and
     * returns the result.
     * @param point Point where function will be evaluated.
     * @return Value returned by the function.
     * @throws Throwable Raised if something failed during the evaluation.
     */
    public double evaluate(double point) throws Throwable;
}
