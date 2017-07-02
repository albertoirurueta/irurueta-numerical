/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.SavitzkyGolayDerivativeEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 8, 2012
 * Copyright 2012. Visual Engineering. All rights reserved.
 */
package com.irurueta.numerical;

import com.irurueta.algebra.*;

/**
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is more robust against small discontinuities than 
 * SymmetricDerivativeEstimator, but it is also slower to compute.
 * This method interpolates the sampled function values into a polynomial of
 * 2nd degree (parabolic), whose derivative is known.
 * Because a linear system of equations has to be solved to determine such
 * polynomial, this method might be less accurate when large values are involved
 * due to limited machine precision.
 */
public class SavitzkyGolayDerivativeEstimator extends DerivativeEstimator{

    /**
     * Number of required point to evaluate to compute derivative
     */
    public static final int N_POINTS = 3;
    
    /**
     * Constructor
     * @param listener listener to evaluate a single dimension function
     */    
    public SavitzkyGolayDerivativeEstimator(
            SingleDimensionFunctionEvaluatorListener listener){
        super(listener);
    }
    
    /**
     * Computes the function derivative at provided point x.
     * @param x Point where derivative is estimated
     * @return Derivative of function at provided point
     * @throws EvaluationException Raised if function cannot be properly 
     * evaluated
     */    
    @Override
    public double derivative(double x) throws EvaluationException{
        //fit a polynomial of degree 2 by evaluating function at x-h, x and x+h
        double h = EPS * Math.abs(x);
        if(h == 0.0) h = EPS; //Trick to reduce finite-precision error
        
        double xh1 = x + h;
        double xh2 = x - h;
        
        double f, fh1, fh2;
        try{
            f = listener.evaluate(x);
            fh1 = listener.evaluate(xh1);
            fh2 = listener.evaluate(xh2);
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
        
        //express the problem as:
        //a * x^2 + b * x + c = f(x)
        //b * xh1^2 + b * xh1 + c = f(xh1)
        //c * xh2^2 + b * xh2 + c = f(xh2)
        
        Matrix a = null;
        try{
            a = new Matrix(N_POINTS, N_POINTS);
        }catch(WrongSizeException ignore){}
        
        double[] b = new double[N_POINTS];
        
        a.setElementAt(0, 0, x * x);
        a.setElementAt(1, 0, xh1 * xh1);
        a.setElementAt(2, 0, xh2 * xh2);
        
        a.setElementAt(0, 1, x);
        a.setElementAt(1, 1, xh1);
        a.setElementAt(2, 1, xh2);
        
        a.setElementAt(0, 2, 1.0);
        a.setElementAt(1, 2, 1.0);
        a.setElementAt(2, 2, 1.0);
        
        //normalize to increase accuracy
        double normA = Utils.normF(a);
        a.multiplyByScalar(1.0 / normA);
        
        b[0] = f;
        b[1] = fh1;
        b[2] = fh2;
        
        //normalize to increase accuracy
        ArrayUtils.multiplyByScalar(b, 1.0 / normA, b);
        
        SingularValueDecomposer decomposer = new SingularValueDecomposer(a);
        
        double aParam, bParam;
        try{
            decomposer.decompose();
            
            //now solve the system of equations in Least Mean Squared Error 
            //because SVD allows the system of equations to be solved using the 
            //pseudo-inverse
            double[] params = decomposer.solve(b);
            aParam = params[0];
            bParam = params[1];
            
        }catch(Throwable t){
            return Double.NaN;
        }
        
        //and c = params[2], but we don't need it
        
        //because we have fitted the function into a polynomial that has
        //expression: a * x^2 + b * x + c, then its derivative is:
        //2.0 * a * x + b, therefore:
        return 2.0 * aParam * x + bParam;
    }
}
