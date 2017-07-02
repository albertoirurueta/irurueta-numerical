/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.DirectionalDerivativeEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * This class evaluates a multidimensional function and obtains its gradient 
 * along a line; such line is defined by an input point and a given direction.
 * Using provided input point and direction, the multidimensional function's 
 * input parameters are determined so they all lay on a line.
 */
public class DirectionalDerivativeEvaluator extends DirectionalEvaluator {
    
    /**
     * Listener to evaluate a multidimensional function's gradient.
     */    
    private GradientFunctionEvaluatorListener gradientListener;
    
    /**
     * Array containing gradient at point p. This is used internally.
     */
    double[] dft;
    
    /**
     * Constructor
     * @param listener Listener to evaluate a multidimensional function
     * @param gradientListener Listener to evaluate a multidimensional 
     * function's gradient
     * @param point Point used as a reference to determine the function's input
     * parameters along a line.
     * @param direction Vector indicating the direction of the line where the 
     * function is evaluated.
     * @throws IllegalArgumentException Raised if point and direction don't have
     * the same length
     */
    public DirectionalDerivativeEvaluator(
            MultiDimensionFunctionEvaluatorListener listener,
            GradientFunctionEvaluatorListener gradientListener, double[] point, 
            double[] direction) throws IllegalArgumentException{
        super(listener, point, direction);
        
        this.gradientListener = gradientListener;        
        dft = new double[point.length];
    }
    
    /**
     * Returns gradient listener that evaluates a multidimensional function 
     * gradient.
     * If the gradient expression is not known (e.g. is not a closed 
     * expression), then a GradientEstimator can be used internally in the 
     * listener implementation
     * @return Gradient listener
     */
    public GradientFunctionEvaluatorListener getGradientListener(){
        return gradientListener;
    }
    
    /**
     * Sets gradient listener that evaluates a multidimensional function 
     * gradient.
     * If the gradient expression is not known (e.g. is not a closed 
     * expression), then a GradientEstimator can be used internally in the 
     * listener implementation being provided.
     * @param gradientListener Gradient listener
     */
    public void setGradientListener(
            GradientFunctionEvaluatorListener gradientListener){
        this.gradientListener = gradientListener;
    }
    
    /**
     * Computes derivative on current direction of a function at distance x from
     * current point and using current listener and gradient listener.
     * @param x Distance from current point using current direction where 
     * function is being evaluated
     * @return Result of evaluating function
     * @throws EvaluationException Thrown if function evaluation fails 
     */
    public double differentiateAt(double x) throws EvaluationException{
        int length = point.length;
        
        for(int i = 0; i < length; i++){
            p[i] = point[i] + x * direction[i];
        }
        
        //Compute gradient at such point
        try{
            gradientListener.evaluateGradient(p, dft);
        }catch(Throwable t){
            throw new EvaluationException(t);
        }
        
        //Obtain 1D derivative on corresponding direction
        double df1 = 0.0;
        for(int i = 0; i < length; i++){
            df1 += dft[i] * direction[i];
        }        
        return df1;
    }
}
