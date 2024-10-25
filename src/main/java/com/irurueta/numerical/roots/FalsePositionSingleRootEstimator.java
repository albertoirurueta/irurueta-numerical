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
 * This class is based on the false position algorithm.
 * This implementation is based on Numerical Recipes 3rd ed. Section 9.2
 * page 451.
 */
public class FalsePositionSingleRootEstimator extends BracketedSingleRootEstimator {

    /**
     * Maximum allowed number of iterations.
     */
    public static final int MAXIT = 30;

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
    public FalsePositionSingleRootEstimator() {
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
     * @param tolerance    Tolerance to find a root. During the estimation of a
     *                     root, if the variation of the estimated root is smaller than the provided
     *                     tolerance, then the algorithm is assumed to be converged, and the root
     *                     is ensured to have an accuracy that equals tolerance.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     *                                      maxEvalPoint.
     * @throws IllegalArgumentException     Raised if provided tolerance is negative.
     */
    public FalsePositionSingleRootEstimator(
            final SingleDimensionFunctionEvaluatorListener listener, final double minEvalPoint,
            final double maxEvalPoint, final double tolerance) throws InvalidBracketRangeException {
        super(listener, minEvalPoint, maxEvalPoint);
        internalSetTolerance(tolerance);
    }

    /**
     * Returns tolerance to find a root. Whenever the variation of the estimated
     * root is smaller than returned tolerance, then the algorithm is assumed to
     * be converged, and the estimated root is ensured to have an accuracy that
     * equals the returned tolerance.
     *
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
     *
     * @param tolerance Tolerance to find a root.
     * @throws LockedException          Raised if this instance is locked while doing
     *                                  some computations.
     * @throws IllegalArgumentException Raised if provided tolerance is negative.
     */
    public void setTolerance(final double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }

    /**
     * Estimates a single root of the provided single dimension function
     * contained within a given bracket of values.
     *
     * @throws LockedException         Exception raised if this instance is already
     *                                 locked.
     * @throws NotReadyException       Exception raised if not enough parameters have
     *                                 been provided in order to start the estimation.
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
            final var v1 = new double[1];
            final var v2 = new double[1];

            final var x1 = minEvalPoint;
            final var x2 = maxEvalPoint;
            final var xacc = tolerance;
            double xl;
            double xh;
            double del;
            var fl = listener.evaluate(x1);
            var fh = listener.evaluate(x2);
            if (fl * fh > 0.0) {
                // root must be bracketed
                locked = false;
                throw new RootEstimationException();
            }
            if (fl < 0.0) {
                xl = x1;
                xh = x2;
            } else {
                xl = x2;
                xh = x1;

                v1[0] = fl;
                v2[0] = fh;
                swap(v1, v2);
                fl = v1[0];
                fh = v2[0];
            }
            var dx = xh - xl;
            for (var j = 0; j < MAXIT; j++) {
                final var rtf = xl + dx * fl / (fl - fh);
                final var f = listener.evaluate(rtf);
                if (f < 0.0) {
                    del = xl - rtf;
                    xl = rtf;
                    fl = f;
                } else {
                    del = xh - rtf;
                    xh = rtf;
                    fh = f;
                }
                dx = xh - xl;
                if (Math.abs(del) < xacc || f == 0.0) {
                    // result obtained
                    root = rtf;
                    rootAvailable = true;
                    locked = false;
                    return;
                }
            }
        } catch (final EvaluationException e) {
            throw new RootEstimationException(e);
        }
        // too many iterations and error exceeds desired tolerance
        locked = false;
        throw new RootEstimationException();
    }

    /**
     * Returns boolean indicating whether this instance is ready to compute a
     * root.
     * This instance is ready as soon as a listener is provided and a bracket
     * is provided or computed.
     *
     * @return True if instance is ready, false otherwise
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isBracketAvailable();
    }

    /**
     * Internal method to set tolerance to find a root. This method does not
     * check whether this instance is locked.
     * Whenever the variation of the estimated root is smaller than provided
     * tolerance, then the algorithm is assumed to be converged, and the
     * estimated root is ensured to have an accuracy that equals provided
     * tolerance.
     *
     * @param tolerance Tolerance to find a root.
     * @throws IllegalArgumentException Raised if provided tolerance is negative.
     */
    private void internalSetTolerance(final double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }
}
