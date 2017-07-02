/**
 * This file contains implementation of
 * com.irurueta.numerical.GradientEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 4, 2012
 */
package com.irurueta.numerical;

/**
 * Class to estimate the gradient of a multidimensional function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the gradient at such point
 */
public class GradientEstimator {
    
    /**
     * Constant considered as machine precision
     */
    public static final double EPS = 1e-8;

    /**
     * Listener to evaluate a multidimensional function.
     */
    public MultiDimensionFunctionEvaluatorListener listener;
    
    /**
     * Internal array to hold input parameter's values.
     */
    private double[] xh;
    
    /**
     * Constructor
     * @param listener Listener to evaluate a multidimensional function
     */
    public GradientEstimator(MultiDimensionFunctionEvaluatorListener listener){
        this.listener = listener;
    }
    
    /**
     * Returns the gradient of a multidimensional function at provided point
     * @param point Input point
     * @return Gradient
     * @throws EvaluationException Raised if function cannot be evaluated.
     */
    public double[] gradient(double[] point) throws EvaluationException{
        double[] result = new double[point.length];
        gradient(point, result);
        return result;
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
    public void gradient(double[] point, double[] result) 
            throws EvaluationException, IllegalArgumentException{
        int length = point.length;
        if(result.length != length) throw new IllegalArgumentException();
        
        if(xh == null || xh.length != length){
            xh = new double[length];
        }
        System.arraycopy(point, 0, xh, 0, length);
        
        try{
            double temp, h, fh;
            double fold = listener.evaluate(point);
            for(int j = 0; j < length; j++){
                temp = point[j];
                h = EPS * Math.abs(temp);
                if(h == 0.0) h = EPS;
                xh[j] = temp + h;
                h = xh[j] - temp;
                fh = listener.evaluate(xh);
                xh[j] = temp;
                result[j] = (fh - fold) / h;
            }
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
    }
}
