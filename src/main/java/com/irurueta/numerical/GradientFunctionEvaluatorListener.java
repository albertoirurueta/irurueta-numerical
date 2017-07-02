/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.GradientFunctionEvaluatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Listener to evaluate/retrieve a multidimensional function's gradient
 */
public interface GradientFunctionEvaluatorListener {
    
    /**
     * Computes/retrieves a multidimensional function's gradient.
     * @param params Array of input parameters of a multidimensional function
     * @param result Array containing estimated gradient. This parameter must 
     * contain an already instantiated array having the same length as params.
     * The values of this gradient will be rewritten after executing this method
     * and this array will contain the estimated or evaluated gradient
     * @throws Throwable Raised if something fails.
     */
    public void evaluateGradient(double[] params, double[] result) 
            throws Throwable;
}
