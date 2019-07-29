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
 * This class searches for a single REAL root on a single dimension function
 * (i.e. f(x) ). The root to be found must be inside a given or computed bracket
 * of values.
 * 
 * This class uses the bisection method, which is one of the slowest but most
 * reliable existing methods (i.e. the method doubles the accuracy on each
 * iteration, so in 64 iterations the highest accuracy provided by a double
 * value can be obtained).
 * 
 * In comparison to other methods, the Bisection one can find roots even in
 * situations where isn't a zero crossing. In other words, bracket estimation
 * might fail for this class, but even then a root might be obtained.
 * 
 * In order to find a root around a given range of values, a bracket of values
 * can be provided. Otherwise, a bracket can be computed by attempting to find
 * a zero-crossing while expanding an initial range of values.
 * 
 * Bracket computation estimates a larger range of values, which can later be
 * refined in order to estimate a given root.
 * 
 * This class is based on the implementation of Numerical Recipes 3rd ed. 
 * Section 9.1.1 page 447.
 */
@SuppressWarnings("WeakerAccess")
public class BisectionSingleRootEstimator extends BracketedSingleRootEstimator {
    
    /**
     * Constant defining maximum number of iterations to estimate a root.
     */
    public static final int JMAX = 50;
    
    /**
     * Constant defining default tolerance to find a root.
     */
    public static final double DEFAULT_TOLERANCE = 1e-6;
    
    /**
     * Minimum allowed tolerance that can be set.
     */
    public static final double MIN_TOLERANCE = 0.0;
    
    /**
     * Tolerance to find a root. Whenever the variation of the estimated root is
     * smaller than the provided tolerance, then the algorithm is assumed to be 
     * converged.
     */
    private double tolerance;
     
    /**
     * Empty constructor.
     */
    public BisectionSingleRootEstimator() {
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
     * @param tolerance Tolerance to find a root. During the estimation of a 
     * root, if the variation of the estimated root is smaller than the provided
     * tolerance, then the algorithm is assumed to be converged, and the root
     * is ensured to have an accuracy that equals tolerance.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     * maxEvalPoint.
     * @throws IllegalArgumentException Raised if provided tolerance is negative.
     */
    public BisectionSingleRootEstimator(
            SingleDimensionFunctionEvaluatorListener listener,
            double minEvalPoint, double maxEvalPoint, double tolerance)
            throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }
    
    /**
     * Returns tolerance to find a root. Whenever the variation of the estimated 
     * root is smaller than returned tolerance, then the algorithm is assumed to
     * be converged, and the estimated root is ensured to have an accuracy that 
     * equals the returned tolerance.
     * @return Tolerance to find a root.
     */
    public double getTolerance() {
        return tolerance;
    }
    
    /**
     * Sets tolerance to find a root. Whenever the variation of the estimated
     * root is smaller than provided tolerance, then the algorithm is assumed to
     * be converged, and the estimated root is ensured to have an accuracy that
     * equals provided tolerance.
     * @param tolerance Tolerance to find a root.
     * @throws LockedException Raised if this instance is locked while doing
     * some computations.
     * @throws IllegalArgumentException Raised if provided tolerance is negative.
     */
    public void setTolerance(double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }
    
    /**
     * Internal method to set tolerance to find a root. This method does not
     * check whether this instance is locked.
     * Whenever the variation of the estimated root is smaller than provided
     * tolerance, then the algorithm is assumed to be converged, and the 
     * estimated root is ensured to have an accuracy that equals provided 
     * tolerance.
     * @param tolerance Tolerance to find a root.
     * @throws IllegalArgumentException Raised if provided tolerance is negative.
     */
    private void internalSetTolerance(double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }
    
    /**
     * Estimates a single root of the provided single dimension function 
     * contained within a given bracket of values.
     * @throws LockedException Exception raised if this instance is already 
     * locked.
     * @throws NotReadyException Exception raised if not enough parameters have
     * been provided in order to start the estimation.
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
        rootAvailable = false;
        
        try {
            double x1 = minEvalPoint;
            double x2 = maxEvalPoint;
            double xacc = tolerance;
            double dx;
            double xmid;
            double rtb;
            double f = listener.evaluate(x1);
            double fmid = listener.evaluate(x2);
        
            if (f * fmid >= 0.0) {
                //check that bracket contains a sign change in function
                locked = false;
                throw new RootEstimationException();
            }
            if (f < 0.0) {
                dx = x2 - x1;
                rtb = x1;
            } else {
                dx = x1 - x2;
                rtb = x2;
            }
            for (int j = 0; j < JMAX; j++) {
                dx *= 0.5;
                xmid = rtb + dx;
                fmid = listener.evaluate(xmid);
                if (fmid <= 0.0) {
                    rtb = xmid;
                }
                if (Math.abs(dx) < xacc || fmid == 0.0) {
                    //result obtained
                    root = rtb;
                    rootAvailable = true;
                    locked = false;
                    return;
                }
            }
        } catch (EvaluationException e) {
            throw new RootEstimationException(e);
        } finally {
            locked = false;
        }
        //too many iterations and error exceeds desired tolerance
        throw new RootEstimationException();
    }
    
    /**
     * Returns boolean indicating whether this instance is ready to compute a
     * root.
     * This instance is ready as soon as a listener is provided and a bracket
     * is provided or computed.
     * @return True if instance is ready, false otherwise
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isBracketAvailable();
    }
}
