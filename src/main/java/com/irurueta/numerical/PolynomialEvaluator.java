/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.PolynomialEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 8, 2016
 */
package com.irurueta.numerical;

import com.irurueta.algebra.Complex;

/**
 * Utility class to evaluate polynomials having either real or complex 
 * coefficients.
 */
public class PolynomialEvaluator {
    
    /**
     * Empty constructor.
     */
    protected PolynomialEvaluator(){}
    
    /**
     * Evaluates polynomial formed by provided polynomial parameters at provided
     * point x. The polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * hence, the array of parameters is [a0, a1, ... a(n-1), an].
     * @param polyParams array of polynomial parameters.
     * @param x point where polynomial is evaluated.
     * @return result of evaluation.
     * @throws IllegalArgumentException if provided array is null or has length 
     * 0.
     */
    public static double evaluate(double[] polyParams, double x) 
            throws IllegalArgumentException{
        if(polyParams == null || polyParams.length == 0){
            throw new IllegalArgumentException();
        }
        
        int length = polyParams.length;
        
        double result = 0.0;
        double powX = 1.0;
        for(int i = length - 1; i >= 0; i--){
            result += polyParams[i]*powX;            
            powX *= x;
        }
        
        return result;
    }

    /**
     * Evaluates polynomial formed by provided polynomial parameters at provided
     * point x. The polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * hence, the array of parameters is [a0, a1, ... a(n-1), an].
     * @param polyParams array of polynomial parameters.
     * @param x point where polynomial is evaluated.
     * @return result of evaluation.
     * @throws IllegalArgumentException if provided array is null or has length 
     * 0.
     */
    public static Complex evaluate(Complex[] polyParams, Complex x) 
            throws IllegalArgumentException{
        if(polyParams == null || polyParams.length == 0){
            throw new IllegalArgumentException();
        }
        
        int length = polyParams.length;
        
        Complex result = new Complex();
        Complex powX = new Complex(1.0, 0.0);
        Complex tmp = new Complex();
        for(int i = length - 1; i >= 0; i--){
            
            polyParams[i].multiply(powX, tmp);
            result.add(tmp);
            
            powX.multiply(x);
        }
                
        return result;
    }
}
