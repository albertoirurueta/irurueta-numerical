/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fitting.LinearFitterMultiDimensionFunctionEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 26, 2015
 */
package com.irurueta.numerical.fitting;

/**
 * Interface to evaluate linear multi dimensional functions
 * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
 * Where the linear function is composed of a linear combination of a basis of
 * functions f0, f1, ... fM
 * For each evaluation at a given point (x1, x2, ...), this interface will 
 * return an array containing the evaluations of the basis functions at such 
 * point f0(x1, x2, ...), f1(x1, x2, ...), ..., fM(x1, x2, ...)
 */
public interface LinearFitterMultiDimensionFunctionEvaluator {
    
    /**
     * Number of dimensions of points (i.e. length of arrays) evaluated by
     * this function evaluator
     * @return number of dimensions of points
     */
    public int getNumberOfDimensions();
    
    /**
     * Creates array where basis function results will be stored
     * @return array where basis function results will be stored
     */
    public double[] createResultArray();
    
    /**
     * Evaluates a linear multi dimension function at provided point and
     * returns the evaluations of the basis functions at such point
     * @param point point where function will be evaluated
     * @param result array where result of evaluation of basis functions is 
     * stored
     * @throws Throwable raised if something failed during the evaluation
     */
    public void evaluate(double[] point, double[] result) throws Throwable;
    
}
