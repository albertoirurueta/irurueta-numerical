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
package com.irurueta.numerical.roots;

import com.irurueta.numerical.*;

/**
 * Finds a single dimensional function's root within a bracket of values using
 * Newton-Raphson's method.
 * This class is based on the implementation found in Numerical Recipes 3rd ed.
 * Secion 9.4. page 456.
 */
@SuppressWarnings("WeakerAccess")
public class NewtonRaphsonSingleRootEstimator 
    extends DerivativeSingleRootEstimator {
    
    /**
     * Maximum number of iterations.
     */
    public static final int JMAX = 20;
    
    /**
     * Constant defining default accuracy of the estimated root.
     */    
    public static final double DEFAULT_TOLERANCE = 1e-6;
    
    /**
     * Constant defining minimum allowed tolerance.
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
    public NewtonRaphsonSingleRootEstimator() {
        super();
        tolerance = DEFAULT_TOLERANCE;
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @param tolerance Tolerance to be achieved in the estimated root.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */        
    public NewtonRaphsonSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener, 
            double minEvalPoint, double maxEvalPoint, double tolerance)
            throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }
        
    /**
     * Constructor.
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     * @param derivativeListener Listener to evaluate the function's derivative.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @param tolerance Tolerance to be achieved in the estimated root.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */            
    public NewtonRaphsonSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener,
            SingleDimensionFunctionEvaluatorListener derivativeListener,
            double minEvalPoint, double maxEvalPoint, double tolerance)
            throws InvalidBracketRangeException {
        super(listener, derivativeListener, minEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }
        
    /**
     * Returns tolerance value.
     * Tolerance is the accuracy to be achieved when estimating a root.
     * If a root is found by this class, it is ensured to have an accuracy below
     * the tolerance value.
     * @return Tolerance value.
     */
    public double getTolerance() {
        return tolerance;
    }

    /**
     * Sets tolerance value.
     * Tolerance is the accuracy to be achieved when estimating a root.
     * If a root is found by this class, it is ensured to have an accuracy below
     * provided tolerance value.
     * @param tolerance Tolerance value.
     * @throws LockedException Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided tolerance value is
     * negative.
     */    
    public void setTolerance(double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }
    
    /**
     * Estimates a local root for a given single dimension function being 
     * evaluated by provided listener.
     * @throws LockedException Exception raised if this instance is already 
     * locked.
     * @throws NotReadyException Exception raised if either a listener has not
     * yet been provided or a bracket has not been provided or computed.
     * @throws RootEstimationException Raised if the root estimation failed for
     * some other reason (usually inability to evaluate the function, 
     * numerical instability or convergence problems, or no roots are found).
     */        
    @Override
    @SuppressWarnings("Duplicates")
    public void estimate() throws LockedException, NotReadyException,
            RootEstimationException {
        
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        rootAvailable = false;
        
        double x1 = minEvalPoint;
        double x2 = maxEvalPoint;
        double xacc = tolerance;
        
        double rtn = 0.5 * (x1 + x2);
        double f;
        double df;
        for (int j = 0; j < JMAX; j++) {
            try {
                f = listener.evaluate(rtn);
                df = derivativeListener.evaluate(rtn);
            } catch (EvaluationException e) {
                throw new RootEstimationException(e);
            }
            
            double dx = f / df;
            rtn -= dx;
            if ((x1 - rtn) * (rtn - x2) < 0.0) {
                //jumped out of brackets
                locked = false;
                throw new RootEstimationException();
            }
            if (Math.abs(dx) < xacc) {
                //root found
                root = rtn;
                rootAvailable = true;
                locked = false;
                return;
            }
        }
        //maximum number of iterations exceeded
        locked = false;
        throw new RootEstimationException();
    }

    /**
     * Internal method to set tolerance value.
     * Tolerance is the accuracy to be achieved when estimating a root.
     * If a root is found by this class, it is ensured to have an accuracy below
     * provided tolerance value.
     * This method does not check whether this instance is locked or not.
     * @param tolerance Tolerance value.
     * @throws IllegalArgumentException Raised if provided tolerance value is
     * negative.
     */
    private void internalSetTolerance(double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }
}
