/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.RealPolynomialEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 8, 2016
 */
package com.irurueta.numerical;

/**
 * Utility class to evaluate real polynomials.
 * This class is useful when the same real polynomial needs to be evaluated 
 * multiple times.
 */
public class RealPolynomialEvaluator extends PolynomialEvaluator{
    
    /**
     * Polynomial coefficients.
     * A polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * Hence, the array of polynomial coefficients is [a0, a1, ... a(n-1), an]
     */
    private double[] mPolyParams;
    
    /**
     * Constructor.
     * @param polyParams polynomial coefficients
     * @throws IllegalArgumentException if provided array is null or has length 
     * 0.
     */
    public RealPolynomialEvaluator(double[] polyParams) 
            throws IllegalArgumentException{
        if(polyParams == null || polyParams.length == 0){
            throw new IllegalArgumentException();
        }
        
        mPolyParams = polyParams;
    }
    
    /**
     * Gets polynomial parameters.
     * @return polynomial parameters.
     */
    public double[] getPolyParams(){
        return mPolyParams;
    }
    
    /**
     * Evaluates polynomial at provided point x.
     * @param x point where polynomial is evaluated.
     * @return result of evaluation.
     */
    public double evaluate(double x){
        return evaluate(mPolyParams, x);
    }
}
