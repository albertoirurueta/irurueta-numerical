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

import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Class to compute local minimum on single dimension functions using a
 * modification of Brent's algorithm that takes into account the function's
 * derivative.
 * This class will search for a local minimum within a bracket of values.
 * A bracket is a set of points: "a" a minimum evaluation point,
 * "b" a middle evaluation point and "c" a maximum evaluation where a &lt;= b
 * &lt;= c, and where f(b) &lt;= f(a) and f(b) &lt;= f(c).
 * This class is based on the implementation of Numerical Recipes 3rd ed.
 * Section 10.4. Page 500.
 */
public class DerivativeBrentSingleOptimizer extends BracketedSingleOptimizer {

    /**
     * Maximum number of iterations to perform. If convergence is not found
     * within this number of iterations, the minimum search will be considered
     * as failed.
     */
    public static final int ITMAX = 100;

    /**
     * Constant defining machine precision.
     */
    public static final double ZEPS = 1e-8;

    /**
     * Default tolerance. Estimated result will be found with an accuracy below
     * or equal to provided tolerance value.
     */
    public static final double DEFAULT_TOLERANCE = 3e-8;

    /**
     * Minimum allowed tolerance value.
     */
    public static final double MIN_TOLERANCE = 0.0;

    /**
     * Listener to evaluate the functions derivative. If the function's
     * derivative is not know (e.g. does not have a closed expression), then
     * a DerivativeEstimator might be used inside the listener implementation.
     */
    private SingleDimensionFunctionEvaluatorListener derivativeListener;

    /**
     * Tolerance. Estimated result will be found with an accuracy below or equal
     * to provided tolerance value.
     */
    private double tolerance;

    /**
     * Empty constructor.
     */
    protected DerivativeBrentSingleOptimizer() {
        super();
        tolerance = DEFAULT_TOLERANCE;
    }

    /**
     * Constructor. Creates an instance with provided bracket of values and a
     * listener to get single dimension function evaluations.
     *
     * @param listener           Listener to evaluate a function.
     * @param derivativeListener Listener to get function derivative.
     * @param minEvalPoint       Minimum bracket evaluation point.
     * @param middleEvalPoint    Middle bracket evaluation point.
     * @param maxEvalPoint       Maximum bracket evaluation point.
     * @param tolerance          tolerance to find result with. Estimated result will be
     *                           found with an accuracy below or equal to provided tolerance value.
     * @throws InvalidBracketRangeException Raised if the following condition is
     *                                      not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     * @throws IllegalArgumentException     Raised if tolerance is negative.
     */
    protected DerivativeBrentSingleOptimizer(
            final SingleDimensionFunctionEvaluatorListener listener,
            final SingleDimensionFunctionEvaluatorListener derivativeListener,
            final double minEvalPoint, final double middleEvalPoint, final double maxEvalPoint,
            final double tolerance) throws InvalidBracketRangeException {
        super(listener, minEvalPoint, middleEvalPoint, maxEvalPoint);
        this.derivativeListener = derivativeListener;
        internalSetTolerance(tolerance);
    }

    /**
     * Returns derivative listener to get function derivative.
     *
     * @return Derivative listener.
     * @throws NotAvailableException Raised if derivative listener is not
     *                               available for retrieval.
     */
    public SingleDimensionFunctionEvaluatorListener getDerivativeListener()
            throws NotAvailableException {
        if (!isDerivativeListenerAvailable()) {
            throw new NotAvailableException();
        }
        return derivativeListener;
    }

    /**
     * Sets derivative listener that gets function derivative.
     *
     * @param derivativeListener Sets derivative listener.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setDerivativeListener(
            final SingleDimensionFunctionEvaluatorListener derivativeListener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        this.derivativeListener = derivativeListener;
    }

    /**
     * Returns boolean indicating whether derivative listener has been provided
     * and is available for retrieval.
     *
     * @return Boolean indicating whether derivative listener is available.
     */
    public boolean isDerivativeListenerAvailable() {
        return derivativeListener != null;
    }

    /**
     * Returns tolerance value. Estimated result will be found with an accuracy
     * below or equal to provided tolerance value.
     *
     * @return Tolerance value.
     */
    public double getTolerance() {
        return tolerance;
    }

    /**
     * Sets tolerance value. Estimated result will be found with an accuracy
     * below or equal to provided tolerance value.
     *
     * @param tolerance Tolerance value.
     * @throws LockedException          Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public void setTolerance(final double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }

    /**
     * This function estimates a function minimum within provided or computed
     * bracket of values.
     * Given a function f that computes a function and also its derivative
     * function df, and given a bracketing triplet of abscissas ax, bx, cx (such
     * that bx is between ax and cx, and f(bx) is less than both f(ax) and
     * f(cx), this routine isolates the minimum to a fractional precision of
     * about tolerance using a modification of Brent's method that uses
     * derivatives. The abscissa of the minimum is returned as xmin and the
     * minimum function value is returned as fmin.
     *
     * @throws LockedException       Raised if this instance is locked, because
     *                               estimation is being computed.
     * @throws NotReadyException     Raised if this instance is not ready because
     *                               either a listener or a bracket has not yet been provided or computed.
     * @throws OptimizationException Raised if the algorithm failed because of
     *                               lack of convergence or because function couldn't be evaluated.
     */
    @SuppressWarnings("DuplicatedCode")
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

        final double[] v1 = new double[1];
        final double[] v2 = new double[2];
        final double[] v3 = new double[3];

        try {
            // Will be used as flags for whether proposed steps are accpetable or
            // not
            boolean ok1;
            boolean ok2;
            double a;
            double b;
            double d = 0.0;
            double d1;
            double d2;
            double du;
            double dv;
            double dw;
            double dx;
            double e = 0.0;
            double fu;
            double fv;
            double fw;
            double fx;
            double olde;
            double tol1;
            double tol2;
            double u;
            double u1;
            double u2;
            double v;
            double w;
            double x;
            double xm;

            // Comments following will point out only differences from the Brent
            // single optimizer. Read that routine first.
            a = Math.min(ax, cx);
            b = Math.max(ax, cx);
            x = w = v = bx;
            fw = fv = fx = listener.evaluate(x);
            dw = dv = dx = derivativeListener.evaluate(x);

            // All out housekeeping chores are doubled by the necessity of moving
            // around derivative values as well as function values
            for (int iter = 0; iter < ITMAX; iter++) {
                xm = 0.5 * (a + b);
                tol1 = tolerance * Math.abs(x) + ZEPS;
                tol2 = 2.0 * tol1;
                if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                    fmin = fx;
                    xmin = x;

                    resultAvailable = true;
                    locked = false;
                    return;
                }

                double tmp = dx >= 0.0 ? a - x : b - x;
                if (Math.abs(e) > tol1) {
                    // Initialize these d's to an out-of-bracket value
                    d1 = 2.0 * (b - a);
                    d2 = d1;
                    // Secant method with one point
                    if (dw != dx) {
                        d1 = (w - x) * dx / (dx - dw);
                    }
                    // And the other
                    if (dv != dx) {
                        d2 = (v - x) * dx / (dx - dv);
                    }
                    // Which of these two estimates of d shall we take? We will
                    // insist that they be within the bracket, and on the side
                    // pointed to by the derivative at x
                    u1 = x + d1;
                    u2 = x + d2;
                    ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
                    ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
                    // Movement on the step before last
                    olde = e;
                    e = d;
                    if (ok1 || ok2) {
                        // Take only an acceptable d, and if both are acceptable,
                        // then take the smallest one.
                        if (ok1 && ok2) {
                            d = Math.abs(d1) < Math.abs(d2) ? d1 : d2;
                        } else if (ok1) {
                            d = d1;
                        } else {
                            d = d2;
                        }

                        if (Math.abs(d) <= Math.abs(0.5 * olde)) {
                            u = x + d;
                            if (u - a < tol2 || b - u < tol2)
                                d = sign(tol1, xm - x);
                        } else {
                            // Bisect, not golden section.
                            e = tmp;
                            d = 0.5 * (e);
                            // Decide which segment by the sign of the derivative
                        }
                    } else {
                        e = tmp;
                        d = 0.5 * e;
                    }
                } else {
                    e = tmp;
                    d = 0.5 * e;
                }

                if (Math.abs(d) >= tol1) {
                    u = x + d;
                    fu = listener.evaluate(u);
                } else {
                    u = x + sign(tol1, d);
                    fu = listener.evaluate(u);
                    if (fu > fx) {
                        // If the minimum step in the downhill direction takes us
                        // uphill, then we are done
                        fmin = fx;
                        xmin = x;

                        resultAvailable = true;
                        locked = false;
                        return;
                    }
                }

                // Now all the housekeeping, sigh
                du = derivativeListener.evaluate(u);
                if (fu <= fx) {
                    if (u >= x) {
                        a = x;
                    } else {
                        b = x;
                    }
                    v1[0] = v;
                    v2[0] = fv;
                    v3[0] = dv;
                    mov3(v1, v2, v3, w, fw, dw);
                    v = v1[0];
                    fv = v2[0];
                    dv = v3[0];


                    v1[0] = w;
                    v2[0] = fw;
                    v3[0] = dw;
                    mov3(v1, v2, v3, x, fx, dx);
                    w = v1[0];
                    fw = v2[0];
                    dw = v3[0];

                    v1[0] = x;
                    v2[0] = fx;
                    v3[0] = dx;
                    mov3(v1, v2, v3, u, fu, du);
                    x = v1[0];
                    fx = v2[0];
                    dx = v3[0];
                } else {
                    if (u < x) {
                        a = u;
                    } else {
                        b = u;
                    }
                    if (fu <= fw || w == x) {
                        v1[0] = v;
                        v2[0] = fv;
                        v3[0] = dv;
                        mov3(v1, v2, v3, w, fw, dw);
                        v = v1[0];
                        fv = v2[0];
                        dv = v3[0];

                        v1[0] = w;
                        v2[0] = fw;
                        v3[0] = dw;
                        mov3(v1, v2, v3, u, fu, du);
                        w = v1[0];
                        fw = v2[0];
                        dw = v3[0];
                    } else if (fu < fv || v == x || v == w) {
                        v1[0] = v;
                        v2[0] = fv;
                        v3[0] = dv;
                        mov3(v1, v2, v3, u, fu, du);
                        v = v1[0];
                        fv = v2[0];
                        dv = v3[0];
                    }
                }

                if (iterationCompletedListener != null) {
                    iterationCompletedListener.onIterationCompleted(this, iter, ITMAX);
                }
            }

        } catch (final EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }

        // Too many iterations in Derivative Brent
        throw new OptimizationException();
    }

    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a local minimum.
     * This instance will be ready once a listener, derivative listener and
     * bracket are available.
     *
     * @return True if ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isDerivativeListenerAvailable() &&
                isBracketAvailable();
    }

    /**
     * Internal method to set tolerance. Estimated result will be found with an
     * accuracy below or equal to provided tolerance value.
     *
     * @param tolerance Tolerance value.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    private void internalSetTolerance(final double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }
}
