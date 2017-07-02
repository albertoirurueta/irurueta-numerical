/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.SymmetricDerivativeEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 11, 2012
 */
package com.irurueta.numerical;

/**
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is more accurate although slower than 
 * DerivativeEstimator
 */
public class SymmetricDerivativeEstimator extends DerivativeEstimator{
    
    /**
     * Constructor
     * @param listener listener to evaluate a single dimension function
     */    
    public SymmetricDerivativeEstimator(
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
        try{
            double h = EPS * Math.abs(x);
            if(h == 0.0) h = EPS; //Trick to reduce finite-precision error
            
            double xh1 = x + h;
            double xh2 = x - h;
            //because of machine precision h could be different in both cases
            double h1 = xh1 - x;
            double h2 = x - xh2;
            
            double hh = h1 + h2; //this is more or less equal to 2.0 * h
            
            double fh1 = listener.evaluate(xh1);
            double fh2 = listener.evaluate(xh2);
            
            return (fh1 - fh2) / hh;
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
    }
}
