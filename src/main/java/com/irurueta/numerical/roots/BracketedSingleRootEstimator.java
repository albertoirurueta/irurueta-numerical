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
 * Computes a root for a single dimension function inside a given bracket of 
 * values, in other words, root will only be searched within provided minimum
 * and maximum evaluation points.
 */
@SuppressWarnings("WeakerAccess")
public abstract class BracketedSingleRootEstimator extends SingleRootEstimator {
    
    /**
     * Number tries to automatically compute a bracket of values for a given
     * function.
     */
    public static final int NTRY = 50;
    
    /**
     * Factor to use to modify initial values when searching for a bracket.
     */
    public static final double FACTOR = 1.6;
    
    /**
     * Default minimum evaluation point.
     */
    public static final double DEFAULT_MIN_EVAL_POINT = -Double.MAX_VALUE;
    
    /**
     * Default maximum evaluation point.
     */
    public static final double DEFAULT_MAX_EVAL_POINT = Double.MAX_VALUE;
    
    /**
     * Constant defining the value by which the largest bracket evaluation value
     * is increased respect the minimum.
     */
    public static final double BRACKET_EPS = 1e-8;
    
    /**
     * Minimum evaluation point.
     */
    protected double minEvalPoint;
    
    /**
     * Maximum evaluation point.
     */
    protected double maxEvalPoint;
    
    /**
     * Boolean indicating whether a bracket has been computed and is available.
     */
    protected boolean bracketAvailable;
    
    /**
     * Constructor.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     */
    public BracketedSingleRootEstimator(double minEvalPoint, 
            double maxEvalPoint) throws InvalidBracketRangeException {
        super();
        internalSetBracket(minEvalPoint, maxEvalPoint);
    }
    
    /**
     * Constructor.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched. The largest value inside the bracket will be
     * Double.MAX_VALUE.
     * @throws InvalidBracketRangeException Raised if minEvalPoint equals the
     * the maximum value a double can contain, which is Double.MAX_VALUE.
     */
    public BracketedSingleRootEstimator(double minEvalPoint) 
            throws InvalidBracketRangeException {
        super();
        internalSetBracket(minEvalPoint, DEFAULT_MAX_EVAL_POINT);
    }
    
    /**
     * Empty constructor.
     */
    public BracketedSingleRootEstimator() {
        super();
        minEvalPoint = DEFAULT_MIN_EVAL_POINT;
        maxEvalPoint = DEFAULT_MAX_EVAL_POINT;
        bracketAvailable = true;
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     */
    public BracketedSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener, 
            double minEvalPoint, double maxEvalPoint) 
            throws InvalidBracketRangeException {
        super(listener);
        internalSetBracket(minEvalPoint, maxEvalPoint);
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched. The largest value inside the bracket will be
     * Double.MAX_VALUE.
     * @throws InvalidBracketRangeException Raised if minEvalPoint equals the
     * the maximum value a double can contain, which is Double.MAX_VALUE.
     */
    public BracketedSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener, 
            double minEvalPoint) throws InvalidBracketRangeException {
        super(listener);
        internalSetBracket(minEvalPoint, DEFAULT_MAX_EVAL_POINT);
    }
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a single dimension function f(x)
     * to find its roots.
     */
    public BracketedSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener) {
        super(listener);
        minEvalPoint = DEFAULT_MIN_EVAL_POINT;
        maxEvalPoint = DEFAULT_MAX_EVAL_POINT;
        bracketAvailable = true;        
    }
    
    /**
     * Sets the bracket of values (i.e. range of values) where the root will be
     * searched.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @throws LockedException Raised if this instance is locked while doing
     * some computations.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     */
    public void setBracket(double minEvalPoint, double maxEvalPoint)
            throws LockedException, InvalidBracketRangeException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetBracket(minEvalPoint, maxEvalPoint);
    }
    
    /**
     * Internal method to set the bracket of values (i.e. range of values) where
     * the root will be searched.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     * root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     * root will be searched.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     */
    private void internalSetBracket(double minEvalPoint, double maxEvalPoint)
            throws InvalidBracketRangeException {
        if (minEvalPoint >= maxEvalPoint) {
            throw new InvalidBracketRangeException();
        }
        
        this.minEvalPoint = minEvalPoint;
        this.maxEvalPoint = maxEvalPoint;
        bracketAvailable = true;
    }
    
    /**
     * Returns boolean indicating whether bracket has been set or not.
     * @return True if bracket has been set, false otherwise.
     */
    public boolean isBracketAvailable() {
        return bracketAvailable;
    }
    
    /**
     * Returns smallest value inside the bracket of values where the root will
     * be searched.
     * @return Smallest value inside the bracket.
     * @throws NotAvailableException Raised if bracket has not been set.
     */
    public double getMinEvaluationPoint() throws NotAvailableException {
        if (!isBracketAvailable()) {
            throw new NotAvailableException();
        }
        return minEvalPoint;
    }
    
    /**
     * Returns largest value inside the bracket of values where the root will
     * be searched.
     * @return Largest values inside the bracket.
     * @throws NotAvailableException Raised if bracket has not been set.
     */
    public double getMaxEvaluationPoint() throws NotAvailableException {
        if (!isBracketAvailable()) {
            throw new NotAvailableException();
        }
        return maxEvalPoint;
    }
    
    /**
     * Starting from provided minimum and maximum values, this method expands
     * the range (i.e. bracket of values) until a zero crossing is found where 
     * a root is present or until the bracket becomes unacceptably large, where
     * an exception will be raised.
     * Notice that this method searches for zero crossings, hence, functions
     * such as Math.pow(x - root, 2.0), will raise a RootEstimationException
     * because the only root present in them does not produce a zero crossing.
     * @param minEvalPoint Smallest initial value to estimate a new larger
     * bracket.
     * @param maxEvalPoint Largest initial value to estimate a new larger 
     * bracket.
     * @throws LockedException Raised if this instance is locked while doing
     * some computations.
     * @throws NotReadyException Raised if this instance is not ready (e.g. a
     * listener has not yet been provided, etc.)
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint
     * @throws RootEstimationException Raised if a bracket couldn't be found
     * inside the provided limits because no zero crossing was present or
     * convergence was not achieved.
     */
    public void computeBracket(double minEvalPoint, double maxEvalPoint)
            throws LockedException, NotReadyException, 
            InvalidBracketRangeException, RootEstimationException {
        
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        if (minEvalPoint >= maxEvalPoint) {
            throw new InvalidBracketRangeException();
        }
        
        locked = true;
        
        //expand initial bracket until function contains a sign change
        double x1 = minEvalPoint;
        double x2 = maxEvalPoint;
        double f1;
        double f2;
        try {
            f1 = listener.evaluate(x1);
            f2 = listener.evaluate(x2);

            boolean found = false;
            for (int j = 0; j < NTRY; j++) {
                if (f1 * f2 < 0.0) {
                    found = true;
                    break;
                }
                if (Math.abs(f1) < Math.abs(f2)) {
                    x1 += FACTOR * (x1 - x2);
                    f1 = listener.evaluate(x1);
                } else {
                    x2 += FACTOR * (x2 - x1);
                    f2 = listener.evaluate(x2);
                }
            }

            if (!found) {
                throw new RootEstimationException();
            }

            this.minEvalPoint = x1;
            this.maxEvalPoint = x2;
            bracketAvailable = true;

        } catch (EvaluationException e) {
            throw new RootEstimationException(e);
        } finally {
            locked = false;
        }
    }
    
    /**
     * Starting from provided point, this method expands the range (i.e. 
     * bracket of values) until a zero crossing is found where a root is present
     * or until the bracket becomes unacceptably large, where an exception will 
     * be raised.
     * Notice that this method searches for zero crossings, hence, functions
     * such as Math.pow(x - root, 2.0), will raise a RootEstimationException
     * because the only root present in them does not produce a zero crossing.
     * @param point Initial value to start the bracket computation. Bracket 
     * range is expanded starting at the point that was provided.
     * @throws LockedException Raised if this instance is locked while doing
     * some computations.
     * @throws NotReadyException Raised if this instance is not ready (e.g. a
     * listener has not yet been provided, etc.).
     * @throws InvalidBracketRangeException Raised if point is close to 
     * Double.MAX_VALUE.
     * @throws RootEstimationException Raised if a bracket couldn't be found
     * inside the provided limits because no zero crossing was present or
     * convergence was not achieved.
     */    
    public void computeBracket(double point) throws LockedException,
            NotReadyException, InvalidBracketRangeException, 
            RootEstimationException {
        computeBracket(point, FACTOR * point + BRACKET_EPS);
    }

    /**
     * Starting at zero, this method expands the range (i.e. bracket of values) 
     * until a zero crossing is found where a root is present or until the 
     * bracket becomes unacceptably large, where an exception will be raised.
     * Notice that this method searches for zero crossings, hence, functions
     * such as Math.pow(x - root, 2.0), will raise a RootEstimationException
     * because the only root present in them does not produce a zero crossing.
     * @throws LockedException Raised if this instance is locked while doing
     * some computations.
     * @throws NotReadyException Raised if this instance is not ready (e.g. a
     * listener has not yet been provided, etc.).
     * @throws RootEstimationException Raised if a bracket couldn't be found
     * inside the provided limits because no zero crossing was present or
     * convergence was not achieved.
     */        
    public void computeBracket() throws LockedException, NotReadyException,
            RootEstimationException {
        try {
            computeBracket(0.0, BRACKET_EPS);
        } catch (InvalidBracketRangeException ignore) {
            //never happens
        }
    }
    
    /**
     * Internal method to swap two values. Value inside a[0] will be swapped 
     * with value provided in b[0].
     * @param a Value to be swapped.
     * @param b Value to be swapped.
     */
    protected void swap(double[] a, double[] b) {
        double tmp = a[0];
        a[0] = b[0];
        b[0] = tmp;
    }
    
    /**
     * Internal method to determine whether a and b have the same sign.
     * @param a Value to be compared.
     * @param b Value to be compared.
     * @return Returns a if a and b have the same sign or -a otherwise.
     */
    protected double sign(double a, double b) {
        if (b >= 0.0) {
            return a >= 0.0 ? a : -a;
        } else {
            return a >= 0.0 ? -a : a;
        }
    }
}
