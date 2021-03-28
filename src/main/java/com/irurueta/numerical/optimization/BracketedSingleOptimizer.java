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
 * This class searches for brackets of values containing a minimum in a single
 * dimension function.
 * A bracket is a set of points: "a" a minimum evaluation point,
 * "b" a middle evaluation point and "c" a maximum evaluation where a &lt;= b
 * &lt;= c, and where f(b) &lt;= f(a) and f(b) &lt;= f(c).
 * This class uses a downhill algorithm that is better suited to continuous
 * functions. Other functions might not obtain reliable results when using this
 * algorithm to obtain a bracket of points.
 * Some subclasses of this class will implement algorithms to refine the
 * solution obtained in a bracket in order to find an accurate estimation of a
 * minimum.
 * Some algorithms might not need to previously compute a bracket and will
 * simply search for a minimum in all the range of possible values, whereas
 * other algorithms will require first the computation of a bracket.
 * In either case, computing a bracket prior estimating a minimum will always
 * ensure that a more reliable solution will be found.
 * Besides, bracket computation is required when a function contains several
 * minima and search of an accurate minimum estimation is desired to be
 * restricted to a certain range of values.
 */
public abstract class BracketedSingleOptimizer extends SingleOptimizer {

    /**
     * The default ratio by which intervals are magnified and.
     */
    public static final double GOLD = 1.618034;

    /**
     * The maximum magnification allowed for a parabolic-fit step.
     */
    public static final double GLIMIT = 100.0;

    /**
     * Small value representing machine precision.
     */
    public static final double TINY = 1e-20;

    /**
     * Default minimum evaluation point where the bracket is supposed to start
     * By default, if no bracket is computed, the whole range of values is used
     * for minimum estimation.
     */
    public static final double DEFAULT_MIN_EVAL_POINT = -Double.MAX_VALUE;

    /**
     * Default middle evaluation point where the bracket is supposed to start
     * By default, if no bracket is computed, the whole range of values is used
     * for minimum estimation.
     */
    public static final double DEFAULT_MIDDLE_EVAL_POINT = 0.0;

    /**
     * Default maximum evaluation point where the bracket is supposed to start
     * By default, if no bracket is computed, the whole range of values is used
     * for minimum estimation.
     */
    public static final double DEFAULT_MAX_EVAL_POINT = Double.MAX_VALUE;

    /**
     * Minimum evaluation point inside the bracket.
     */
    protected double ax;

    /**
     * Middle evaluation point inside the bracket.
     */
    protected double bx;

    /**
     * Maximum evaluation point inside the bracket.
     */
    protected double cx;

    /**
     * Boolean indicating whether a bracket has been provided or computed.
     */
    private boolean bracketAvailable;

    /**
     * Function evaluation value at minimum evaluation point inside the bracket.
     */
    private double fa;

    /**
     * Function evaluation value at middle evaluation point inside the bracket.
     */
    private double fb;

    /**
     * Function evaluation value at maximum evaluation point inside the bracket.
     */
    private double fc;

    /**
     * Boolean indicating whether function evaluation at bracket limits and
     * middle point are available or not.
     */
    private boolean bracketEvaluationAvailable;

    /**
     * Constructor. Creates an instance with provided bracket of values.
     *
     * @param minEvalPoint    Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint    Maximum bracket evaluation point.
     * @throws InvalidBracketRangeException Raised if the following condition is
     *                                      not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     */
    protected BracketedSingleOptimizer(final double minEvalPoint,
                                       final double middleEvalPoint,
                                       final double maxEvalPoint)
            throws InvalidBracketRangeException {
        internalSetBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);
    }

    /**
     * Empty Constructor. Creates an instance using default bracket values.
     */
    protected BracketedSingleOptimizer() {
        ax = DEFAULT_MIN_EVAL_POINT;
        bx = DEFAULT_MIDDLE_EVAL_POINT;
        cx = DEFAULT_MAX_EVAL_POINT;
        bracketAvailable = true;
    }

    /**
     * Constructor. Creates an instance with provided bracket of values and a
     * listener to get single dimension function evaluations.
     *
     * @param listener        Listener to evaluate a function.
     * @param minEvalPoint    Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint    Maximum bracket evaluation point.
     * @throws InvalidBracketRangeException Raised if the following condition is
     *                                      not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     */
    protected BracketedSingleOptimizer(
            final SingleDimensionFunctionEvaluatorListener listener,
            final double minEvalPoint, final double middleEvalPoint, final double maxEvalPoint)
            throws InvalidBracketRangeException {
        super(listener);
        internalSetBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);
    }

    /**
     * Sets a bracket of values to later search for a minimum. A local minimum
     * will only be search within the minimum and maximum evaluation points of
     * a given bracket.
     * If bracket is not provided, it can also be computed from a default or
     * coarse set of points in order to obtain a more refined bracket so that
     * a minimum search can be estimated more precisely.
     *
     * @param minEvalPoint    Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint    Maximum bracket evaluation point.
     * @throws InvalidBracketRangeException Raised if the following condition is
     *                                      not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     * @throws LockedException              Raised if this instance is locked. This instance
     *                                      will be locked while doing some operations. Attempting to change any
     *                                      parameter while being locked will raise this exception.
     */
    public void setBracket(final double minEvalPoint, final double middleEvalPoint,
                           final double maxEvalPoint) throws LockedException,
            InvalidBracketRangeException {

        if (isLocked()) {
            throw new LockedException();
        }
        internalSetBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);
    }

    /**
     * Returns boolean indicating whether a bracket has been provided or
     * computed and is available for retrieval.
     *
     * @return true if a bracket has been provided, false otherwise.
     */
    public boolean isBracketAvailable() {
        return bracketAvailable;
    }

    /**
     * Returns minimum evaluation point where the bracket starts
     *
     * @return Minimum evaluation point.
     * @throws NotAvailableException Raised if not provided or computed.
     */
    public double getMinEvaluationPoint() throws NotAvailableException {
        if (!isBracketAvailable()) {
            throw new NotAvailableException();
        }

        return ax;
    }

    /**
     * Returns middle evaluation point within the bracket.
     *
     * @return Middle evaluation point.
     * @throws NotAvailableException Raised if not provided or computed.
     */
    public double getMiddleEvaluationPoint() throws NotAvailableException {
        if (!isBracketAvailable()) {
            throw new NotAvailableException();
        }

        return bx;
    }

    /**
     * Returns maximum evaluation point whether the bracket finishes.
     *
     * @return Maximum evaluation point.
     * @throws NotAvailableException Raised if not provided or computed.
     */
    public double getMaxEvaluationPoint() throws NotAvailableException {
        if (!isBracketAvailable()) {
            throw new NotAvailableException();
        }

        return cx;
    }

    /**
     * Returns single dimension function evaluation at provided or computed
     * minimum evaluation point where the bracket starts.
     *
     * @return Function evaluation at bracket's minimum evaluation point.
     * @throws NotAvailableException Raised if bracket evaluations are not
     *                               available.
     */
    public double getEvaluationAtMin() throws NotAvailableException {
        if (!areBracketEvaluationsAvailable()) {
            throw new NotAvailableException();
        }

        return fa;
    }

    /**
     * Returns single dimension function evaluation at provided or computed
     * middle evaluation point within the bracket.
     *
     * @return Function evaluation at bracket's middle evaluation point.
     * @throws NotAvailableException Raised if bracket evaluations are not
     *                               available.
     */
    public double getEvaluationAtMiddle() throws NotAvailableException {
        if (!areBracketEvaluationsAvailable()) {
            throw new NotAvailableException();
        }

        return fb;
    }

    /**
     * Returns single dimension function evaluation at provided or computed
     * maximum evaluation point where the bracket finishes.
     *
     * @return Function evaluation at bracket's maximum evaluation point.
     * @throws NotAvailableException Raised if bracket evaluations are not
     *                               available.
     */
    public double getEvaluationAtMax() throws NotAvailableException {
        if (!areBracketEvaluationsAvailable()) {
            throw new NotAvailableException();
        }

        return fc;
    }

    /**
     * Computes a bracket of values using provided values as a starting point.
     * Given a function f, and given distinct initial points ax and bx, this
     * routine searches in the downhill direction (defined by the function as
     * evaluated at the initial points) and returns.
     * ax (minimum evaluation point), bx (middle evaluation point), cx (maximum
     * evaluation point) that bracket a minimum of the function. Also returned
     * are the function values at the three points fa, fb, and fc, which are the
     * function evaluations at minimum, middle and maximum bracket points.
     *
     * @param minEvalPoint    Initial minimum evaluation point of bracket.
     * @param middleEvalPoint Initial middle evaluation point of bracket.
     * @throws LockedException              Raised if this instance is locked. This instance
     *                                      will be locked while doing some operations. Attempting to change any
     *                                      parameter while being locked will raise this exception.
     * @throws NotReadyException            Raised if this instance is not ready because a
     *                                      listener has not yet been provided.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt;
     *                                      middleEvalPoint.
     * @throws OptimizationException        Raised if a bracket couldn't be found .
     *                                      because convergence was not achieved or function evaluation failed.
     */
    @SuppressWarnings("Duplicates")
    public void computeBracket(final double minEvalPoint, final double middleEvalPoint)
            throws LockedException, NotReadyException,
            InvalidBracketRangeException, OptimizationException {

        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }
        if (minEvalPoint > middleEvalPoint) {
            throw new InvalidBracketRangeException();
        }

        locked = true;

        final double[] a = new double[1];
        final double[] b = new double[1];
        final double[] c = new double[1];

        try {
            ax = minEvalPoint;
            bx = middleEvalPoint;
            double fu;
            fa = listener.evaluate(ax);
            fb = listener.evaluate(bx);

            //switch roles of a and b so that we can go downhill in the
            //direction from a to b
            if (fb > fa) {
                a[0] = ax;
                b[0] = bx;
                swap(a, b);
                ax = a[0];
                bx = b[0];

                a[0] = fa;
                b[0] = fb;
                swap(a, b);
                a[0] = fa;
                b[0] = fb;
            }

            //First guess for c
            cx = bx + GOLD * (bx - ax);
            fc = listener.evaluate(cx);

            //Keep returning here until we bracket.
            while (fb > fc) {
                //Compute u by parabolic extrapolation from a, b, c. TINY is
                //used to prevent any possible division by zero.
                final double r = (bx - ax) * (fb - fc);
                final double q = (bx - cx) * (fb - fa);
                double u = bx - ((bx - cx) * q - (bx - ax) * r) /
                        (2.0 * sign(Math.max(Math.abs(q - r), TINY), q - r));
                final double ulim = bx + GLIMIT * (cx - bx);

                //We won't go farther than this. Test various possibilities:
                if ((bx - u) * (u - cx) > 0.0) {
                    //Parabolic u is between b and c: try it.
                    fu = listener.evaluate(u);
                    if (fu < fc) {
                        //Got a minimum between b and c.
                        ax = bx;
                        bx = u;
                        fa = fb;
                        fb = fu;
                        break;

                    } else if (fu > fb) {
                        //Got a minimum between a and u
                        cx = u;
                        fc = fu;
                        break;
                    }

                    //Parabolic fit was no use. Use default magnification.
                    u = cx + GOLD * (cx - bx);
                    fu = listener.evaluate(u);

                } else if ((cx - u) * (u - ulim) > 0.0) {
                    //Parabolic fit is between c and its allowed limit
                    fu = listener.evaluate(u);

                    if (fu < fc) {
                        a[0] = bx;
                        b[0] = cx;
                        c[0] = u;
                        shft3(a, b, c, u + GOLD * (u - cx));
                        bx = a[0];
                        cx = b[0];
                        u = c[0];

                        a[0] = fb;
                        b[0] = fc;
                        c[0] = fu;
                        shft3(a, b, c, listener.evaluate(u));
                        fb = a[0];
                        fc = b[0];
                        fu = c[0];
                    }
                } else if ((u - ulim) * (ulim - cx) >= 0.0) {
                    //Limit parabolic u to maximum allowed value.
                    u = ulim;
                    fu = listener.evaluate(u);
                } else {
                    //Reject parabolic u, use default magnification
                    u = cx + GOLD * (cx - bx);
                    fu = listener.evaluate(u);
                }
                //Eliminate oldest point and continue
                a[0] = ax;
                b[0] = bx;
                c[0] = cx;
                shft3(a, b, c, u);
                ax = a[0];
                bx = b[0];
                cx = c[0];

                a[0] = fa;
                b[0] = fb;
                c[0] = fc;
                shft3(a, b, c, fu);
                fa = a[0];
                fb = b[0];
                fc = c[0];
            }
        } catch (final EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }

        bracketAvailable = true;
        bracketEvaluationAvailable = true;
    }

    /**
     * Computes a bracket of values using provided value as a starting point,
     * and assuming that bracket finishes at Double.MAX_VALUE.
     * Given a function f, and given distinct initial points ax and bx = 0.0,
     * this routine searches in the downhill direction (defined by the function
     * as evaluated at the initial points) and returns
     * ax (minimum evaluation point), bx (middle evaluation point), cx (maximum
     * evaluation point) that bracket a minimum of the function. Also returned
     * are the function values at the three points fa, fb, and fc, which are the
     * function evaluations at minimum, middle and maximum bracket points.
     *
     * @param minEvalPoint Initial minimum evaluation point of bracket.
     * @throws LockedException              Raised if this instance is locked. This instance
     *                                      will be locked while doing some operations. Attempting to change any
     *                                      parameter while being locked will raise this exception.
     * @throws NotReadyException            Raised if this instance is not ready because a
     *                                      listener has not yet been provided.
     * @throws InvalidBracketRangeException Raised if minEvalPoint &lt; 0.0.
     * @throws OptimizationException        Raised if a bracket couldn't be found
     *                                      because convergence was not achieved or function evaluation failed.
     */
    public void computeBracket(final double minEvalPoint) throws LockedException,
            NotReadyException, OptimizationException,
            InvalidBracketRangeException {
        computeBracket(minEvalPoint, DEFAULT_MIDDLE_EVAL_POINT);
    }

    /**
     * Computes a bracket of values using the whole range of possible values as
     * an initial guess.
     * Given a function f, and given distinct initial points ax =
     * -Double.MAX_VALUE and bx = 0.0, this
     * routine searches in the downhill direction (defined by the function as
     * evaluated at the initial points) and returns
     * ax (minimum evaluation point), bx (middle evaluation point), cx (maximum
     * evaluation point) that bracket a minimum of the function. Also returned
     * are the function values at the three points fa, fb, and fc, which are the
     * function evaluations at minimum, middle and maximum bracket points
     *
     * @throws LockedException       Raised if this instance is locked. This instance
     *                               will be locked while doing some operations. Attempting to change any
     *                               parameter while being locked will raise this exception.
     * @throws NotReadyException     Raised if this instance is not ready because a
     *                               listener has not yet been provided.
     * @throws OptimizationException Raised if a bracket couldn't be found
     *                               because convergence was not achieved or function evaluation failed.
     */
    public void computeBracket() throws LockedException, NotReadyException,
            OptimizationException {
        try {
            computeBracket(DEFAULT_MIN_EVAL_POINT, DEFAULT_MIDDLE_EVAL_POINT);
        } catch (InvalidBracketRangeException ignore) {
            //never happens
        }
    }

    /**
     * Computes function evaluations at provided or estimated bracket locations.
     * After calling this method bracket evaluations will be available.
     *
     * @throws LockedException       Raised if this instance is locked. This instance
     *                               will be locked while doing some operations. Attempting to change any
     *                               parameter while being locked will raise this exception.
     * @throws NotReadyException     Raised if this instance is not ready because a
     *                               listener has not yet been provided.
     * @throws OptimizationException Raised if function evaluation failed.
     */
    public void evaluateBracket() throws LockedException, NotReadyException,
            OptimizationException {

        if (isLocked()) {
            throw new LockedException();
        }
        if (!isReady()) {
            throw new NotReadyException();
        }

        locked = true;

        try {
            fa = listener.evaluate(ax);
            fb = listener.evaluate(bx);
            fc = listener.evaluate(cx);
        } catch (final EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }

        bracketEvaluationAvailable = true;
    }

    /**
     * Returns boolean indicating whether bracket evaluations are available for
     * retrieval.
     *
     * @return True if bracket evaluations are available, false otherwise.
     */
    public boolean areBracketEvaluationsAvailable() {
        return bracketEvaluationAvailable;
    }

    /**
     * Internal method to determine whether a and b have the same sign.
     *
     * @param a Value to be compared.
     * @param b Value to be compared.
     * @return Returns a if a and b have the same sign or -a otherwise.
     */
    protected double sign(final double a, final double b) {
        if (b >= 0.0) {
            return a >= 0.0 ? a : -a;
        } else {
            return a >= 0.0 ? -a : a;
        }
    }

    /**
     * Pushes b value into a, and c value into b. a and b are in/out parameters.
     * Results will be available at a[0] and b[0] after executing this method.
     *
     * @param a a value to be lost.
     * @param b a value to be shifted into a.
     * @param c a value to be shifted into b.
     */
    protected void shft2(final double[] a, final double[] b, final double c) {
        a[0] = b[0];
        b[0] = c;
    }

    /**
     * Pushes b value into a, and c value into b and d value into c. a, b and c
     * are in/out parameters.
     * Results will be available at a[0], b[0] and c[0] after executing this
     * method.
     *
     * @param a a value to be lost.
     * @param b a value to be shifted into a.
     * @param c a value to be shifted into b.
     * @param d a value to be shifted into c.
     */
    protected void shft3(final double[] a, final double[] b, final double[] c, final double d) {
        a[0] = b[0];
        b[0] = c[0];
        c[0] = d;
    }

    /**
     * Moves d, e and f into a[0], b[0] and c[0]. Previously existing values
     * into a, b, c will be lost after executing this method.
     *
     * @param a a value to be set.
     * @param b a value to be set.
     * @param c a value to be set.
     * @param d a value to be copied.
     * @param e a value to be copied.
     * @param f a value to be copied.
     */
    protected void mov3(final double[] a, final double[] b, final double[] c, final double d,
                        final double e, final double f) {
        a[0] = d;
        b[0] = e;
        c[0] = f;
    }

    /**
     * Internal method to swap two values. Value inside a[0] will be swapped
     * with value provided in b[0].
     *
     * @param a Value to be swapped.
     * @param b Value to be swapped.
     */
    private void swap(final double[] a, final double[] b) {
        double tmp = a[0];
        a[0] = b[0];
        b[0] = tmp;
    }

    /**
     * Internal method to set a bracket of values. This method does not check
     * whether this instance is locked.
     *
     * @param minEvalPoint    Minimum bracket evaluation point.
     * @param middleEvalPoint Middle bracket evaluation point.
     * @param maxEvalPoint    Maximum bracket evaluation point.
     * @throws InvalidBracketRangeException Raised if the following condition is
     *                                      not met: minEvalPoint &lt;= middleEvalPoint &lt;= maxEvalPoint.
     */
    private void internalSetBracket(final double minEvalPoint, final double middleEvalPoint,
                                    final double maxEvalPoint) throws InvalidBracketRangeException {

        if ((minEvalPoint > middleEvalPoint) ||
                (middleEvalPoint > maxEvalPoint)) { //which also means || (minEvalPoint > maxEvalPoint))
            throw new InvalidBracketRangeException();
        }

        ax = minEvalPoint;
        bx = middleEvalPoint;
        cx = maxEvalPoint;

        bracketAvailable = true;
    }
}
