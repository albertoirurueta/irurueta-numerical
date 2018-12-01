/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.irurueta.numerical.optimization;

import com.irurueta.numerical.*;

/**
 * Class to find a minimum on a multidimensional function along a given line
 * of input values. The difference between this abstract class and 
 * LineMultiOptimizer is that subclasses of this class will use the function
 * gradient information. By using gradient information, typically convergence is
 * faster, although, if gradient does not have a closed expression, then a
 * GradientEstimator will be needed in the gradient listener provided to this 
 * class. Notice that a GradientEstimator might estimate gradients with a 
 * certain error, so depending on the function topology, LineMultiOptimizer
 * subclasses might obtain greater accuracy than subclasses of this class.
 */
@SuppressWarnings("WeakerAccess")
public abstract class DerivativeLineMultiOptimizer extends MultiOptimizer {
    /**
     * n-dimensional point containing a minimum in a given line.
     */
    protected double[] p;
    
    /**
     * Direction to make the search.
     */
    protected double[] xi;
    
    /**
     * Number of dimensions on function being evaluated.
     */    
    private int n;
    
    /**
     * Listener to evaluate the function's gradient. If the function's 
     * gradient is not know (e.g. does not have a closed expression), then
     * a GradientEstimator might be used inside the listener implementation.
     */    
    protected GradientFunctionEvaluatorListener gradientListener;
    
    /**
     * Class in charge of evaluating a function through a given line.
     */    
    private DirectionalDerivativeEvaluator evaluator;
    
    /**
     * Internal optimizer to find a minimum of a function along aline of
     * input values. Hence, input is converted to a single dimension using a
     * DirectionalEvaluator.
     */
    private DerivativeBrentSingleOptimizer dbrent;
    
    /**
     * Empty constructor.
     */
    public DerivativeLineMultiOptimizer() {
        super();
        p = xi = null;
        n = 0;
        gradientListener = null;
        evaluator = null;
        dbrent = null;
    }  
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     */
    public DerivativeLineMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener) {
        super(listener);
        p = xi = null;
        n = 0;
        gradientListener = null;
        evaluator = null;
        dbrent = null;        
    }
    
    /**
     * Constructor.
     * @param gradientListener Listener to evaluate the function's gradient.
     */
    public DerivativeLineMultiOptimizer(
            GradientFunctionEvaluatorListener gradientListener) {
        super();
        p = xi = null;
        n = 0;
        this.gradientListener = gradientListener;
        evaluator = null;
        dbrent = null;        
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     * @param gradientListener Listener to evaluate the function's gradient.
     */
    public DerivativeLineMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener,
            GradientFunctionEvaluatorListener gradientListener) {
        super(listener);
        p = xi = null;
        n = 0;
        this.gradientListener = gradientListener;
        evaluator = null;
        dbrent = null;        
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimension function.
     * @param gradientListener Listener to evaluate the function's gradient.
     * @param point Start point where algorithm will be started. Start point 
     * should be close to the local minimum to be found. Provided array must 
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param direction Direction to start looking for a minimum. Provided array
     * must have the same length as the number of dimensions of the function 
     * being evaluated. Provided direction is considered as a vector pointing
     * to the minimum to be found.
     * @throws IllegalArgumentException Raised if provided point and direction
     * don't have the same length.
     */    
    public DerivativeLineMultiOptimizer(
            MultiDimensionFunctionEvaluatorListener listener,
            GradientFunctionEvaluatorListener gradientListener,
            double[] point, double[] direction) 
            throws IllegalArgumentException {
        super(listener);
        internalSetStartPointAndDirection(point, direction);
        n = 0;
        this.gradientListener = gradientListener;
        evaluator = null;
        dbrent = null;        
    }
    
    /**
     * Returns boolean indicating whether start point has already been provided
     * and is ready for retrieval.
     * @return True if available, false otherwise.
     */

    public boolean isStartPointAvailable() {
        return p != null;
    }
    
    /**
     * Returns start point where algorithm will be started. Start point should 
     * be close to the local minimum to be found.
     * @return Start point where algorithm will be started.
     * @throws NotAvailableException Raised if start point has not yet been
     * provided and is not available.
     */    
    public double[] getStartPoint() throws NotAvailableException {
        if (!isStartPointAvailable()) {
            throw new NotAvailableException();
        }
        return p;
    }
    
    /**
     * Returns boolean indicating whether direction has already been provided 
     * and is ready for retrieval.
     * @return True if available, false otherwise.
     */    
    public boolean isDirectionAvailable() {
        return xi != null;
    }
    
    /**
     * Returns direction to start looking for a minimum. Provided array must 
     * have the same length as the number of dimensions of the function being 
     * evaluated. Provided direction is considered as a vector pointing to the 
     * minimum to be found.
     * @return Direction to start looking for a minimum.
     * @throws NotAvailableException Raised if direction has not yet been 
     * provided and is not available.
     */    
    public double[] getDirection() throws NotAvailableException {
        if (!isDirectionAvailable()) {
            throw new NotAvailableException();
        }
        return xi;
    }
    
    /**
     * Internal method to set start point and direction to start the search for
     * a local minimum.
     * This method does not check whether this instance is locked.
     * @param point Start point where algorithm will be started. Start point 
     * should be close to the local minimum to be found. Provided array must 
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param direction Direction to start looking for a minimum. Provided array
     * must have the same length as the number of dimensions of the function 
     * being evaluated. Provided direction is considered as a vector pointing
     * to the minimum to be found.
     * @throws IllegalArgumentException Raised if provided point and direction
     * don't have the same length.
     */    
    private void internalSetStartPointAndDirection(double[] point, 
            double[] direction) throws IllegalArgumentException {
        if (point.length != direction.length) {
            throw new IllegalArgumentException();
        }
        
        p = point;
        xi = direction;        
    }
    
    /**
     * Internal method to set start point and direction to start the search for
     * a local minimum.
     * @param point Start point where algorithm will be started. Start point 
     * should be close to the local minimum to be found. Provided array must 
     * have a length equal to the number of dimensions of the function being
     * evaluated, otherwise and exception will be raised when searching for the
     * minimum.
     * @param direction Direction to start looking for a minimum. Provided array
     * must have the same length as the number of dimensions of the function 
     * being evaluated. Provided direction is considered as a vector pointing
     * to the minimum to be found.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided point and direction
     * don't have the same length.
     */    
    public void setStartPointAndDirection(double[] point, double[] direction)
            throws LockedException, IllegalArgumentException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetStartPointAndDirection(point, direction);
    }
    
    /**
     * Returns boolean indicating whether this instance is considered to be
     * ready to start the estimation of a minimum.
     * This instance is considered to be ready once a listener, gradient 
     * listener, start point and direction are provided.
     * @return True if this instance is ready, false otherwise.
     */    
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isGradientListenerAvailable() &&
                isStartPointAvailable() && isDirectionAvailable();
    }
    
    /**
     * Returns gradient listener.
     * The gradient listener is used to evaluate the function's gradient. If the
     * function's gradient is not know (e.g. does not have a closed expression),
     * then a GradientEstimator might be used inside the listener 
     * implementation.
     * @return Gradient listener.
     * @throws NotAvailableException Raised if gradient listener is not 
     * available for retrieval.
     */
    public GradientFunctionEvaluatorListener getGradientListener()
            throws NotAvailableException {
        if (!isGradientListenerAvailable()) {
            throw new NotAvailableException();
        }
        return gradientListener;
    }
    
    /**
     * Sets gradient listener.
     * The gradient listener is used to evaluate the function's gradient. If the
     * function's gradient is not know (e.g. does not have a closed expression),
     * then a GradientEstimator might be used inside the listener 
     * implementation.
     * @param gradientListener Gradient listener.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setGradientListener(
            GradientFunctionEvaluatorListener gradientListener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.gradientListener = gradientListener;
    }
    
    /**
     * Returns boolean indicating whether the gradient listener has been 
     * provided and is available for retrieval.
     * @return True if available, false otherwise.
     */
    public boolean isGradientListenerAvailable() {
        return gradientListener != null;
    }       
    
    /**
     * Searches for a minimum along a given line of input values.
     * The line being searched is obtained by using a start point and direction.
     * @return Returns function evaluation at minimum that has been found.
     */
    @SuppressWarnings("Duplicates")
    protected double linmin() {
        double ax, xx, linxmin;
        n = p.length;
        
        if (evaluator == null) {
            //attempt to reuse evaluator
            evaluator = new DirectionalDerivativeEvaluator(listener, 
                    gradientListener, p, xi);
        }
        if (evaluator.getListener() != listener) {
            //update listener
            evaluator.setListener(listener);
        }
        if (evaluator.getGradientListener() != gradientListener) {
            //update gradient listener
            evaluator.setGradientListener(gradientListener);
        }
        if (evaluator.getPoint() != p || evaluator.getDirection() != xi) {
            evaluator.setPointAndDirection(p, xi);
        }
        
        ax = 0.0;
        xx = 1.0;
                
        try {
            if (dbrent == null) {
                //attempt to reuse brent single optimizer
                dbrent = new DerivativeBrentSingleOptimizer(
                        new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(double point) throws EvaluationException {
                        return evaluator.evaluateAt(point);
                    }
                }, new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(double point) throws EvaluationException {
                        return evaluator.differentiateAt(point);
                    }
                }, DerivativeBrentSingleOptimizer.DEFAULT_MIN_EVAL_POINT,
                    DerivativeBrentSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT,
                    DerivativeBrentSingleOptimizer.DEFAULT_MAX_EVAL_POINT,
                    DerivativeBrentSingleOptimizer.DEFAULT_TOLERANCE);
            }
            
            dbrent.computeBracket(ax, xx);
            dbrent.minimize();
            linxmin = dbrent.getResult();
            
            for (int j = 0; j < n; j++) {
                xi[j] *= linxmin;
                p[j] += xi[j];
            }
            
            return dbrent.getEvaluationAtResult();
        } catch (Throwable t) {
            //if minimization fails we assume that obtained result is the worst
            //possible one
            return Double.MAX_VALUE;
        }
    }
}
