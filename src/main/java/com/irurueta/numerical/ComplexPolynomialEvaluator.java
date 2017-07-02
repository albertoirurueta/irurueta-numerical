/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.ComplexPolynomialEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 8, 2016
 */
package com.irurueta.numerical;

import com.irurueta.algebra.Complex;

/**
 * Utility class to evaluate complex polynomials.
 * This class is useful when the same complex polynomial needs to be evaluated
 * multiple times.
 */
public class ComplexPolynomialEvaluator extends PolynomialEvaluator{
    
    /**
     * Polynomial coefficients.
     * A polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * Hence, the array of polynomial coefficients is [a0, a1, ... a(n-1), an]
     */
    private Complex[] mPolyParams;
    
    /**
     * Constructor.
     * @param polyParams polynomial coefficients
     * @throws IllegalArgumentException if provided array is null or has length 
     * 0.
     */
    public ComplexPolynomialEvaluator(Complex[] polyParams) 
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
    public Complex[] getPolyParams(){
        return mPolyParams;
    }
    
    /**
     * Evaluates polynomial at provided point x.
     * @param x point where polynomial is evaluated.
     * @return result of evaluation.
     */
    public Complex evaluate(Complex x){
        return evaluate(mPolyParams, x);
    }    
}
