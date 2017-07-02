/**
 * @file
 * This file contains implementation of 
 * com.irurueta.numerical.SymmetricGradientEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 11, 2012
 */
package com.irurueta.numerical;

/**
 * Class to estimate the gradient of a multidimensional function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the gradient at such point
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is more accurate although slower than 
 * GradientEstimator
 */
public class SymmetricGradientEstimator extends GradientEstimator{
    
    /**
     * Internal array containing one point to sample close to the original one
     */
    private double[] xh1;
    
    /**
     * Internal array containing one point to sample close to the original one
     */    
    private double[] xh2;
    
    /**
     * Constructor
     * @param listener Listener to evaluate a multidimensional function
     */   
    public SymmetricGradientEstimator(
            MultiDimensionFunctionEvaluatorListener listener){
        super(listener);
    }
    
    /**
     * Sets estimated gradient in provided result array of a multidimensional
     * function at provided point.
     * This method is preferred respect to gradient(double[]) because result
     * array can be reused and hence is more memory efficient
     * @param point Input point
     * @param result Output parameter containing estimated array. This parameter
     * must be an array of length equal to point
     * @throws EvaluationException Raised if function cannot be evaluated
     * @throws IllegalArgumentException Raised if length of result and point are
     * not equal.
     */    
    @Override
    public void gradient(double[] point, double[] result) 
            throws EvaluationException, IllegalArgumentException{
        int n = point.length;
        if(result.length != n) throw new IllegalArgumentException();
        
        if(xh1 == null || xh1.length != n){ 
            xh1 = new double[n];
            System.arraycopy(point, 0, xh1, 0, n);
        }
        if(xh2 == null || xh2.length != n){ 
            xh2 = new double[n];
            System.arraycopy(point, 0, xh2, 0, n);
        }
                
        try{
            double temp, h, h1, h2, hh, fh1, fh2;
            //double fold = listener.evaluate(point);
            for(int j = 0; j < n; j++){
                temp = point[j];
                h = EPS * Math.abs(temp);
                if(h == 0.0) h = EPS; //Trich to reduce finite-precision error
                xh1[j] = temp + h;
                xh2[j] = temp - h;
                //because of machine precision h could be different in both cases
                
                h1 = xh1[j] - temp;
                h2 = temp - xh2[j];
                
                hh = h1 + h2; //this is more or less equal to 2.0 * h
                
                fh1 = listener.evaluate(xh1);
                fh2 = listener.evaluate(xh2);
                
                xh1[j] = temp;
                xh2[j] = temp;
                
                result[j] = (fh1 - fh2) / hh;
            }
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
    }    
}
