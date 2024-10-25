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
 * This implementation is based on Numerical Recipes 3rd ed. Section 9.2.1, page
 * 452.
 */
public class RidderSingleRootEstimator extends BracketedSingleRootEstimator {

    /**
     * Maximum number of iterations.
     */
    public static final int MAXIT = 60;

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
    public RidderSingleRootEstimator() {
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
    public RidderSingleRootEstimator(
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

        final var x1 = minEvalPoint;
        final var x2 = maxEvalPoint;
        final var xacc = tolerance;
        double fl;
        double fh;
        try {
            fl = listener.evaluate(x1);
            fh = listener.evaluate(x2);
        } catch (final EvaluationException e) {
            throw new RootEstimationException(e);
        }

        double ans;
        var found = false;
        if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
            var xl = x1;
            var xh = x2;
            ans = -9.99e99;
            try {
                for (var j = 0; j < MAXIT; j++) {
                    final var xm = 0.5 * (xl + xh);
                    final var fm = listener.evaluate(xm);
                    final var s = Math.sqrt(fm * fm - fl * fh);
                    if (s == 0.0) {
                        found = true;
                        break;
                    }
                    final var xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
                    if (Math.abs(xnew - ans) <= xacc) {
                        // result found
                        found = true;
                        break;
                    }
                    ans = xnew;
                    final var fnew = listener.evaluate(ans);
                    if (sign(fm, fnew) != fm) {
                        xl = xm;
                        fl = fm;
                        xh = ans;
                        fh = fnew;
                    } else if (sign(fl, fnew) != fl) {
                        xh = ans;
                        fh = fnew;
                    } else if (sign(fh, fnew) != fh) {
                        xl = ans;
                        fl = fnew;
                    } else {
                        // never get here
                        locked = false;
                        throw new RootEstimationException();
                    }
                    if (Math.abs(xh - xl) <= xacc) {
                        // result found
                        found = true;
                        break;
                    }

                }
            } catch (final EvaluationException e) {
                throw new RootEstimationException(e);
            }
            if (!found) {
                // too many iterations and error exceeds desired tolerance
                locked = false;
                throw new RootEstimationException();
            }
        } else {
            if (fl == 0.0) {
                // result found
                ans = x1;
            } else if (fh == 0.0) {
                // result found
                ans = x2;
            } else {
                locked = false;
                throw new RootEstimationException();
            }
        }

        // result found
        root = ans;
        rootAvailable = true;
        locked = false;
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
}
