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

import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Computes a root for a single dimension function inside a given bracket of
 * values, in other words, root will only be searched within provided minimum
 * and maximum evaluation points.
 * This class searches for REAL roots only!
 * This implementation is based on Numerical Recipes 3rd ed. Section 9.4, page
 * 456.
 */
public class SafeNewtonRaphsonSingleRootEstimator extends DerivativeSingleRootEstimator {

    /**
     * Maximum number of iterations.
     */
    public static final int MAXIT = 100;

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
    public SafeNewtonRaphsonSingleRootEstimator() {
        super();
        tolerance = DEFAULT_TOLERANCE;
    }

    /**
     * Constructor.
     *
     * @param listener     Listener to evaluate a single dimension function f(x)
     *                     to find its roots.
     * @param minEvalPoint Smallest value inside the bracket of values where the
     *                     root will be searched.
     * @param maxEvalPoint Largest value inside the bracket of values where the
     *                     root will be searched.
     * @param tolerance    Tolerance to be achieved in the estimated root.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     *                                      maxEvalPoint.
     * @throws IllegalArgumentException     Raised if tolerance is negative.
     */
    public SafeNewtonRaphsonSingleRootEstimator(
            final SingleDimensionFunctionEvaluatorListener listener, final double minEvalPoint,
            final double maxEvalPoint, final double tolerance) throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }

    /**
     * Constructor.
     *
     * @param listener           Listener to evaluate a single dimension function f(x)
     *                           to find its roots.
     * @param derivativeListener Listener to evaluate the function's derivative
     * @param minEvalPoint       Smallest value inside the bracket of values where the
     *                           root will be searched.
     * @param maxEvalPoint       Largest value inside the bracket of values where the
     *                           root will be searched.
     * @param tolerance          Tolerance to be achieved in the estimated root.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     *                                      maxEvalPoint.
     * @throws IllegalArgumentException     Raised if tolerance is negative.
     */
    public SafeNewtonRaphsonSingleRootEstimator(
            final SingleDimensionFunctionEvaluatorListener listener,
            final SingleDimensionFunctionEvaluatorListener derivativeListener, final double minEvalPoint,
            final double maxEvalPoint, final double tolerance) throws InvalidBracketRangeException {
        super(listener, derivativeListener, minEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }

    /**
     * Returns tolerance value.
     * Tolerance is the accuracy to be achieved when estimating a root.
     * If a root is found by this class, it is ensured to have an accuracy below
     * the tolerance value.
     *
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
     *
     * @param tolerance Tolerance value.
     * @throws IllegalArgumentException Raised if provided tolerance value is
     *                                  negative.
     */
    private void internalSetTolerance(final double tolerance) {
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
     *
     * @param tolerance Tolerance value.
     * @throws LockedException          Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided tolerance value is
     *                                  negative.
     */
    public void setTolerance(final double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }

    /**
     * Estimates a local root for a given single dimension function being
     * evaluated by provided listener.
     *
     * @throws LockedException         Exception raised if this instance is already
     *                                 locked.
     * @throws NotReadyException       Exception raised if either a listener has not
     *                                 yet been provided or a bracket has not been provided or computed.
     * @throws RootEstimationException Raised if the root estimation failed for
     *                                 some other reason (usually inability to evaluate the function,
     *                                 numerical instability or convergence problems, or no roots are found).
     */
    @Override
    @SuppressWarnings("Duplicates")
    public void estimate() throws LockedException, NotReadyException, RootEstimationException {

        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        locked = true;
        rootAvailable = false;

        final var x1 = minEvalPoint;
        final var x2 = maxEvalPoint;
        final var xacc = tolerance;

        double xh;
        double xl;
        final double fl;
        final double fh;
        try {
            fl = listener.evaluate(x1);
            fh = listener.evaluate(x2);
        } catch (final EvaluationException e) {
            throw new RootEstimationException(e);
        }

        if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
            // root must be bracketed
            locked = false;
            throw new RootEstimationException();
        }
        if (fl == 0.0) {
            root = x1;
            rootAvailable = true;
            locked = false;
            return;
        }
        if (fh == 0.0) {
            root = x2;
            rootAvailable = true;
            locked = false;
            return;
        }
        if (fl < 0.0) {
            xl = x1;
            xh = x2;
        } else {
            xh = x1;
            xl = x2;
        }
        var rts = 0.5 * (x1 + x2);
        var dxold = Math.abs(x2 - x1);
        var dx = dxold;
        double f;
        double df;
        try {
            f = listener.evaluate(rts);
            df = derivativeListener.evaluate(rts);
        } catch (final EvaluationException e) {
            throw new RootEstimationException(e);
        }

        for (var j = 0; j < MAXIT; j++) {
            if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (Math.abs(2.0 * f) > Math.abs(dxold * df))) {
                dxold = dx;
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
                if (xl == rts) {
                    root = rts;
                    rootAvailable = true;
                    locked = false;
                    return;
                }
            } else {
                dxold = dx;
                dx = f / df;
                final var temp = rts;
                rts -= dx;
                if (temp == rts) {
                    root = rts;
                    rootAvailable = true;
                    locked = false;
                    return;
                }
            }
            if (Math.abs(dx) < xacc) {
                root = rts;
                rootAvailable = true;
                locked = false;
                return;
            }

            try {
                f = listener.evaluate(rts);
                df = derivativeListener.evaluate(rts);
            } catch (final EvaluationException e) {
                throw new RootEstimationException(e);
            }

            if (f < 0.0) {
                xl = rts;
            } else {
                xh = rts;
            }
        }
        // maximum number of iterations exceeded
        locked = false;
        throw new RootEstimationException();
    }
}
