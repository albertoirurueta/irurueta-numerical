/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fitting.LevenbergMarquardtMultiVariateFunctionEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date June 1, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.Matrix;

/**
 * Interface to evaluate non-linear multi variate and multi dimensional 
 * functions.
 * Evaluation of functions requires both function value at provided point x and
 * function jacopiab respect to its parameters (i.e. derivatives respect to its
 * parameters for each function output or variable)
 */
public interface LevenbergMarquardtMultiVariateFunctionEvaluator {
    
    /**
     * Number of dimensions of points (i.e. length of input points arrays) 
     * evaluated by this function evaluator
     * @return number of dimensions of points
     */
    public int getNumberOfDimensions();
    
    /**
     * Number of variables of function f. This is equal to the length of the
     * array obtained as function evaluations. Hence, a function f can
     * be expressed as f = [f1, f2, ... fN], and the number of variables would
     * be N
     * @return number of variables of function f
     */
    public int getNumberOfVariables();
    
    /**
     * Creates array where estimated parameters will be stored.
     * This array MUST contain the initial guessed solution for the LEvenberg-
     * Marquardt algorithm
     * @return array where estimated parameters will be stored
     */
    public double[] createInitialParametersArray();
    
    /**
     * Evaluates a non-linear multi variate function at provided point using 
     * provided parameters and returns its evaluation and jacobian of the 
     * function respect the function parameters
     * @param i number of sample being evaluated
     * @param point point where function will be evaluated
     * @param result result of function evaluation. Its length is equal to the
     * number of function variables
     * @param params initial parameters estimation to be tried. These will 
     * change as the Levenberg-Marquard algorithm iterates to the best solution.
     * These are used as input parameters along with point to evaluate function
     * @param jacobian jacobian containing partial derivatives of the function 
     * respect to each provided parameter for each function output or variable
     * @throws Throwable raised if something failed during the evaluation
     */
    public void evaluate(int i, double[] point, double[] result, 
            double[] params, Matrix jacobian) throws Throwable;    
    
}
