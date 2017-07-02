/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.MultiVariateFunctionEvaluatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 31, 2015
 */
package com.irurueta.numerical;

/**
 * Interface to define how multivariate functions can be evaluated.
 */
public interface MultiVariateFunctionEvaluatorListener {
    /**
     * Evaluates a multi variate function such as f1(x1, x2, ...), 
     * f2(x1, x2, ...) at providd multidimensional point and returns the result
     * as a vectorial value
     * @param point multidimensional point where function will be evaluated
     * @param result vector where function evaluation will be stored
     * @throws Throwable if something failed during the evaluation
     */
    public void evaluate(double[] point, double[] result) throws Throwable;
    
    /**
     * Number of variables of function f. This is equal to the length of the
     * array obtained as function evaluations. Hence, a function f can
     * be expressed as f = [f1, f2, ... fN], and the number of variables would
     * be N
     * @return number of variables of function f
     */
    public int getNumberOfVariables();
}
