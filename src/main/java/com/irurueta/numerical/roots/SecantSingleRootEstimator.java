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
 * This implementation is based on Numerical Recipes 3rd ed. Section 9.2, page
 * 449
 */
public class SecantSingleRootEstimator extends BracketedSingleRootEstimator {

    /**
     * Maximum number of iterations.
     */
    public static final int MAXIT = 30;

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
    public SecantSingleRootEstimator() {
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
    public SecantSingleRootEstimator(
            final SingleDimensionFunctionEvaluatorListener listener, final double minEvalPoint,
            final double maxEvalPoint, final double tolerance) throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
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
    public void estimate() throws LockedException, NotReadyException, RootEstimationException {
        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        locked = true;
        rootAvailable = false;

        try {
            final var x1 = minEvalPoint;
            final var x2 = maxEvalPoint;
            final var xacc = tolerance;
            double xl;
            double rts;
            var fl = listener.evaluate(x1);
            var f = listener.evaluate(x2);
            final var v1 = new double[1];
            final var v2 = new double[1];

            if (Math.abs(fl) < Math.abs(f)) {
                rts = x1;
                xl = x2;
                v1[0] = fl;
                v2[0] = f;
                swap(v1, v2);
                fl = v1[0];
                f = v2[0];
            } else {
                xl = x1;
                rts = x2;
            }
            for (int j = 0; j < MAXIT; j++) {
                final var dx = (xl - rts) * f / (f - fl);
                xl = rts;
                fl = f;
                rts += dx;
                f = listener.evaluate(rts);
                if (Math.abs(dx) < xacc || f == 0.0) {
                    // Result obtained
                    root = rts;
                    rootAvailable = true;
                    locked = false;
                    return;
                }
            }
        } catch (final EvaluationException e) {
            throw new RootEstimationException(e);
        } finally {
            locked = false;
        }
        // too many iterations and error exceeds desired tolerance
        throw new RootEstimationException();
    }

    /**
     * Returns boolean indicating whether this instance is ready to start
     * estimating a root.
     * This class will be ready once a listener is provided and a bracket is
     * either provided or computed.
     *
     * @return True if this instance is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isBracketAvailable();
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
}
