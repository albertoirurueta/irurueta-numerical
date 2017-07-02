/**
 * This file contains implementation of
 * com.irurueta.numerical.JacobianEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 31, 2015
 */
package com.irurueta.numerical;

import com.irurueta.algebra.Matrix;

/**
 * Class to estimate the Jacobian of a multi variate and multidimensional 
 * function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the Jacobian at such point
 */
public class JacobianEstimator {
    
    /**
     * Constant considered as machine precision
     */
    public static final double EPS = 1e-8;
    
    /**
     * Listener to evaluate a multivariate function.
     */
    public MultiVariateFunctionEvaluatorListener listener;
    
    /**
     * Internal array to hold input parameter's values
     */
    private double[] xh;
    
    /**
     * Constructor
     * @param listener listener to evaluate a multivariate function
     */
    public JacobianEstimator(MultiVariateFunctionEvaluatorListener listener){
        this.listener = listener;
    }
    
    /**
     * Returns the Jacobian of a multivariate function at provided point
     * @param point input point
     * @return jacobian
     * @throws EvaluationException raised if function cannot be evaluated
     */
    public Matrix jacobian(double[] point) throws EvaluationException{
        try{
            Matrix result = new Matrix(listener.getNumberOfVariables(), 
                    point.length);
            jacobian(point, result);
            return result;
        }catch(EvaluationException e){
            throw e;
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
    }
    
    /**
     * Sets estimated jacobian in provided result matrix of a multivariate 
     * function at provided point.
     * This method is preferred respect to jacobian(double[]) because result
     * matrix can be reused and hence is more memory efficient
     * @param point input point
     * @param result output matrix containing estimated jacobian. This parameter
     * must be a matrix having the same number of columns as the length of the
     * point and the same number of rows as the number of variables returned by
     * function evaluations
     * @throws EvaluationException raised if function cannot be evaluated
     * @throws IllegalArgumentException if size of result is not valid
     */
    public void jacobian(double[] point, Matrix result) 
            throws EvaluationException, IllegalArgumentException{
        int numdims = point.length;
        int numvars = listener.getNumberOfVariables();
        if(result.getColumns() != numdims) throw new IllegalArgumentException();
        if(result.getRows() != numvars) throw new IllegalArgumentException();
        
        double[] fold = new double[numvars];
        double[] fh = new double[numvars];
        if(xh == null || xh.length != numdims){
            xh = new double[numdims];            
        }
        System.arraycopy(point, 0, xh, 0, numdims);
        
        try{
            double temp, h;
            listener.evaluate(point, fold);
            for(int j = 0; j < numdims; j++){            
                temp = point[j];
                h = EPS * Math.abs(temp);
                if(h == 0.0) h = EPS;
                xh[j] = temp + h;
                h = xh[j] - temp;
                listener.evaluate(xh, fh);
                xh[j] = temp;
                for(int i = 0; i < numvars; i++){
                    result.setElementAt(i, j, (fh[i] - fold[i]) / h);
                }
            }
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
    }
}
