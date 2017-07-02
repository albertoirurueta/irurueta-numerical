/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fittin.LevenbergMarquardtSingleDimensionFunctionEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 27, 2015
 */
package com.irurueta.numerical.fitting;

/**
 * Interface to evaluate non-linear single dimensional functions.
 * Evaluation of functions requires both function value at provided point x and
 * function gradient respect to its parameters (i.e. derivatives respect to its
 * parameters)
 */
public interface LevenbergMarquardtSingleDimensionFunctionEvaluator {
    
    /**
     * Creates array where estimated parameters will be stored.
     * This array MUST contain the initial guessed solution for the Levenberg-
     * Marquardt algorithm
     * @return array where estimated parameters will be stored
     */
    public double[] createInitialParametersArray();
    
    /**
     * Evaluates a non-linear single dimension function at provided point using 
     * provided parameters and returns its evaluation and derivatives of the 
     * function respect the function parameters
     * @param i number of sample being evaluated
     * @param point point where function is evaluated
     * @param params initial parameters estimation to be tried. These will 
     * change as the Levenberg-Marquard algorithm iterates to the best solution.
     * These are used as input parameters along with point to evaluate function
     * @param derivatives partial derivatives of the function respect to each
     * provided parameter
     * @return function evaluation at provided point and using provided 
     * parameters
     * @throws Throwable raised if something failed during the evaluation
     */
    public double evaluate(int i, double point, double[] params, 
            double[] derivatives) throws Throwable;
}
