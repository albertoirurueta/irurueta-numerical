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

import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * This class for a single dimensional function's local minimum.
 * This class is based in the Golden search algorithm found in
 * Numerical Recipes 3rd ed. Section 10.2 page 492.
 */
@SuppressWarnings({"WeakerAccess", "Duplicates"})
public class GoldenSingleOptimizer extends BracketedSingleOptimizer {
    
    /**
     * Golden ratio.
     */
    public static final double R = 0.61803399;
    
    /**
     * Golden ratio.
     */
    public static final double C = 1.0 - R;
    
    /**
     * Constant defining the default accuracy of the estimated minimum.
     */
    public static final double DEFAULT_TOLERANCE = 3e-8;
    
    /**
     * Minimum allowed tolerance.
     */    
    public static final double MIN_TOLERANCE = 0.0;
    
    /**
     * Tolerance value. The algorithm will iterate until the result converges
     * below this value of accuracy or until the maximum number of iterations is
     * achieved (and in such case, convergence will be assumed to have failed).
     */
    private double tolerance;
    
    /**
     * Empty constructor.
     */    
    public GoldenSingleOptimizer() {
        super();
        tolerance = DEFAULT_TOLERANCE;
    }
    
    /**
     * Constructor. Creates an instance with provided bracket of values and a
     * listener to get single dimension function evaluations.
     * @param listener Listener to evaluate a function.
     * @param minEvalPoint Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint Maximum bracket evaluation point.
     * @param tolerance Tolerance or accuracy to be obtained in estimated 
     * minimum.
     * @throws InvalidBracketRangeException Raised if the following condition is
     * not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public GoldenSingleOptimizer(
            SingleDimensionFunctionEvaluatorListener listener,
            double minEvalPoint, double middleEvalPoint, double maxEvalPoint,
            double tolerance) throws InvalidBracketRangeException, 
            IllegalArgumentException {
        super(listener, minEvalPoint, middleEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }
    
    /**
     * Returns tolerance value, which is the accuracy to be obtained when a
     * minimum is estimated.
     * The algorithm will iterate until the result converges below this value of
     * accuracy or until the maximum number of iterations is achieved (and in 
     * such case, convergence will be assumed to have failed).
     * @return Tolerance value.
     */
    public double getTolerance() {
        return tolerance;
    }
    
    /**
     * Sets algorithm's tolerance.
     * The algorithm will iterate until the result converges below this value of
     * accuracy or until the maximum number of iterations is achieved (an in 
     * such case, convergence will be assumed to have failed).
     * @param tolerance Tolerance or accuracy to be obtained in estimated 
     * minimum.
     * @throws LockedException Raised if this instance is locked. This instance
     * will be locked while doing some operations. Attempting to change any 
     * parameter while being locked will raise this exception.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */    
    public void setTolerance(double tolerance) throws LockedException,
            IllegalArgumentException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }
    
    /**
     * This function estimates a function minimum within provided or computed
     * bracket of values.
     * Given a function f, and given a bracketing triplet of abscissas ax, bx,
     * cx (such that bx is between ax and cx, and f(bx) is less than both f(ax)
     * and f(cx), this routine isolates the minimum to a fractional prevision of
     * about tolerance using Brent's method. The abscissa of the minimum is 
     * returned as xmin, and the function value of the minimum is returned as
     * fmin, the returned function value.
     * @throws LockedException Raised if this instance is locked, because
     * estimation is being computed.
     * @throws NotReadyException Raised if this instance is not ready because
     * either a listener or a bracket has not yet been provided or computed.
     * @throws OptimizationException Raised if the algorithm failed because of
     * lack of convergence or because function couldn't be evaluated.
     */    
    @Override
    public void minimize() throws LockedException, NotReadyException,
            OptimizationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        
        double[] v1 = new double[1];
        double[] v2 = new double[2];
        double[] v3 = new double[3];

        
        try {
            //At any given time we will keep track of four points x0, x1, x2, x3
            double x1, x2;
            double x0 = ax;
            double x3 = cx;
            
            //Make x0 to x1 the smaller segment, and fill in the new point to be
            //tried
            if (Math.abs(cx - bx) > Math.abs(bx - ax)) {
                x1 = bx;
                x2 = bx + C * (cx - bx);
            } else {
                x2 = bx;
                x1 = bx - C * (bx - ax);
            }
            
            //The initial function evaluations. Note that we never need to
            //evaluate the function at the original endpoints
            double f1 = listener.evaluate(x1);
            double f2 = listener.evaluate(x2);
            
            while (Math.abs(x3 - x0) > tolerance * (Math.abs(x1) + 
                    Math.abs(x2))) {
                if (f2 < f1) {
                    //One possible outcome, its housekeeping and a new function
                    //evaluation
                    v1[0] = x0;
                    v2[0] = x1;
                    v3[0] = x2;
                    shft3(v1, v2, v3, R * x2 + C * x3);
                    x0 = v1[0];
                    x1 = v2[0];
                    x2 = v3[0];
                    
                    v1[0] = f1;
                    v2[0] = f2;
                    shft2(v1, v2, listener.evaluate(x2));
                    f1 = v1[0];
                    f2 = v2[0];
                } else {
                    //The other outcome, and its new function evaluation
                    v1[0] = x3;
                    v2[0] = x2;
                    v3[0] = x1;
                    shft3(v1, v2, v3, R * x1 + C * x0);
                    x3 = v1[0];
                    x2 = v2[0];
                    x1 = v3[0];
                    
                    v1[0] = f2;
                    v2[0] = f1;
                    shft2(v1, v2, listener.evaluate(x1));
                    f2 = v1[0];
                    f1 = v2[0];
                }
                //Back to see if we are done.
            }
            
            //We are done. Output the best of the current values
            if (f1 < f2) {
                xmin = x1;
                fmin = f1;
            } else {
                xmin = x2;
                fmin = f2;
            }
        } catch (Throwable t) {
            locked = false;
            throw new OptimizationException(t);
        }
        
        resultAvailable = true;
        locked = false;
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a minimum or not.
     * The instance is ready when both the listener and the bracket are 
     * available.
     * @return True if this instance is ready, false otherwise.
     */    
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isBracketAvailable();
    }
    
    /**
     * Internal method to set algorithm tolerance. This method does not check
     * whether this instance is locked or not.
     * The algorithm will iterate until the result converges below this value of
     * accuracy or until the maximum number of iterations is achieved (and in
     * such case, convergence will be assumed to have failed).
     * @param tolerance Tolerance or accuracy to be obtained in estimated 
     * minimum.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */    
    private void internalSetTolerance(double tolerance) 
            throws IllegalArgumentException {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }    
}
