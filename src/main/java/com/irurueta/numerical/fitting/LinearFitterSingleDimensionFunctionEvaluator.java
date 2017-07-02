/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.LinearFitterSingleDimensionFunctionEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 24, 2015
 */
package com.irurueta.numerical.fitting;

/**
 * Interface to evaluate linear single dimensional functions 
 * f(x) = a * f0(x) + b * f1(x) + ...
 * Where the linear function is composed of a linear combination of a basis of
 * functions f0, f1, ... fM
 * For each evaluation at a given point x, this interface will return an array
 * containing the evaluations of the basis functions at such point f0(x), f1(x),
 * ..., fM(x)
 */
public interface LinearFitterSingleDimensionFunctionEvaluator {
    
    /**
     * Creates array where basis function results will be stored
     * @return array where basis function results will be stored
     */
    public double[] createResultArray();
    
    /**
     * Evaluates a linear single dimension function at provided point and
     * returns the evaluations of the basis functions at such point
     * @param point point where function will be evaluated
     * @param result array where result of evaluation of basis functions is 
     * stored
     * @throws Throwable raised if something failed during the evaluation
     */
    public void evaluate(double point, double[] result) throws Throwable;
    
}
