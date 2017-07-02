/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.DerivativeEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is faster although less accurate than 
 * SymmetricDerivativeEstimator
 */
public class DerivativeEstimator {
    
    /**
     * Constant defining machine precision for this algorithm.
     */
    public static final double EPS = 1e-8;
    
    /**
     * Listener to evaluate a single dimension function.
     */
    protected SingleDimensionFunctionEvaluatorListener listener;
    
    /**
     * Constructor
     * @param listener listener to evaluate a single dimension function
     */
    public DerivativeEstimator(
            SingleDimensionFunctionEvaluatorListener listener){
        this.listener = listener;
    }
    
    /**
     * Computes the function derivative at provided point x.
     * @param x Point where derivative is estimated
     * @return Derivative of function at provided point
     * @throws EvaluationException Raised if function cannot be properly 
     * evaluated
     */
    public double derivative(double x) throws EvaluationException{
        try{
            double fold = listener.evaluate(x);
        
            double h = EPS * Math.abs(x);
            if(h == 0.0) h = EPS; //Trick to reduce finite-precision error
            double xh = x + h;
            h = xh - x;
            double fh = listener.evaluate(xh);
        
            return (fh - fold) / h;
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
    }
}
