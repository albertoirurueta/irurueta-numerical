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

import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * This class estimates the root of a single dimension continuous function using
 * Brent's method.
 * The implementation of this class is based on Numerical Recipes 3rd ed. 
 * Section 9.3. Page 454.
 */
@SuppressWarnings("WeakerAccess")
public class BrentSingleRootEstimator extends BracketedSingleRootEstimator {
    
    /**
     * Constant defining maximum number of iterations.
     */
    public static final int ITMAX = 100;
    
    /**
     * Constant defining a small value which is considered as machine precision.
     */
    public static final double EPS = 1e-10;
    
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
    public BrentSingleRootEstimator() {
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
    public BrentSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener,
            double minEvalPoint, double maxEvalPoint, double tolerance)
            throws InvalidBracketRangeException, IllegalArgumentException {
        super(listener, minEvalPoint, maxEvalPoint);
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
     * Internal method to set tolerance value.
     * Tolerance is the accuracy to be achieved when estimating a root.
     * If a root is found by this class, it is ensured to have an accuracy below
     * provided tolerance value.
     * This method does not check whether this instance is locked or not.
     * @param tolerance Tolerance value.
     * @throws IllegalArgumentException Raised if provided tolerance value is
     * negative.
     */
    private void internalSetTolerance(double tolerance) 
            throws IllegalArgumentException {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
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
    public void setTolerance(double tolerance) throws LockedException,
            IllegalArgumentException {
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
    public void estimate() throws LockedException, NotReadyException,
            RootEstimationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        
        locked = true;
        
        double x1 = minEvalPoint;
        double x2 = maxEvalPoint;
        double tol = tolerance;
        double a = x1, b = x2, c = x2, d = 0.0, e = 0.0, fc, p, q,
                r, s, tol1, xm;
        double fa, fb;
        try {
            fa = listener.evaluate(a);
            fb = listener.evaluate(b);
        } catch (Throwable t) {
            throw new RootEstimationException(t);
        }
        
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
            //root must be bracketed
            locked = false;
            throw new RootEstimationException();
        }
        fc = fb;
        for (int iter = 0; iter < ITMAX; iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a;
                fc = fa;
                e = d = b - a;
            }
            if (Math.abs(fc) < Math.abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol1 = 2.0 * EPS * Math.abs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);
            if (Math.abs(xm) <= tol1 || fb == 0.0) {
                //root found
                root = b;
                rootAvailable = true;
                locked = false;
                return;
            }
            if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
                s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0) {
                    q = -q;
                }
                p = Math.abs(p);
                double min1 = 3.0 * xm * q - Math.abs(tol1 * q);
                double min2 = Math.abs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            if (Math.abs(d) > tol1) {
                b += d;
            } else {
                b += sign(tol1, xm);
            }
            try {
                fb = listener.evaluate(b);
            } catch (Throwable t) {
                throw new RootEstimationException(t);
            }
        }
        //maximum number of iterations exceeded
        locked = false;
        throw new RootEstimationException();
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to start 
     * estimating a root.
     * This class will be ready once a listener is provided and a bracket is
     * either provided or computed.
     * @return True if this instance is ready, false otherwise.
     */
    @Override
    public boolean isReady(){
        return isListenerAvailable() && isBracketAvailable();
    }
}
