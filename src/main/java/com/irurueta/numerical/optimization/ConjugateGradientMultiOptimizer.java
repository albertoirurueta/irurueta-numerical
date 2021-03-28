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
import com.irurueta.numerical.GradientFunctionEvaluatorListener;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;

/**
 * This class searches for a multi dimension function local minimum.
 * The local minimum is searched by starting the algorithm at a start point
 * and a given direction which should be close and point to the local minimum to
 * be found to achieve the best accuracy with the lowest number of iterations.
 * NOTE: this algorithm might not have proper convergence in some situations,
 * but it is ensured to provide faster convergence than other algorithms such
 * as Brent because gradient information is also provided.
 * The implementation of this class is based on Numerical Recipes 3rd ed.
 * Section 10.8 page 515.
 */
@SuppressWarnings({"WeakerAccess", "DuplicatedCode"})
public class ConjugateGradientMultiOptimizer extends LineMultiOptimizer {
    /**
     * Constant defining default tolerance or accuracy to be achieved on the
     * minimum being estimated by this class.
     */
    public static final double DEFAULT_TOLERANCE = 3e-8;

    /**
     * Minimum allowed tolerance value.
     */
    public static final double MIN_TOLERANCE = 0.0;

    /**
     * Maximum allowed iterations.
     */
    public static final int ITMAX = 200;

    /**
     * Constant defining a value to be considered as machine precision.
     */
    public static final double EPS = 1e-18;

    /**
     * Convergence criterion for the zero gradient test.
     */
    public static final double GTOL = 1e-8;

    /**
     * Defines whether Polak-Ribiere is used if true, otherwise Fletcher-Reeves
     * will be used.
     */
    public static final boolean DEFAULT_USE_POLAK_RIBIERE = true;

    /**
     * The fractional tolerance in the function value such that failure to
     * decrease by more than this amount on one iteration signals doneness.
     */
    private double tolerance;

    /**
     * Member contains number of iterations that were needed to estimate a
     * minimum.
     */
    private int iter;

    /**
     * Value of the function at the minimum.
     */
    private double fret;

    /**
     * Listener to obtain gradient values for the multi dimension function being
     * evaluated.
     * If the gradient is unknown (e.g. doesn't have a closed expression), the
     * provided listener could use a GradientEstimator to obtain one.
     */
    private GradientFunctionEvaluatorListener gradientListener;

    /**
     * Boolean indicating whether Polak-Ribiere method is used if true,
     * otherwise Fletcher-Reeves will be used.
     */
    private boolean usePolakRibiere;

    /**
     * Empty constructor.
     */
    public ConjugateGradientMultiOptimizer() {
        super();
        tolerance = DEFAULT_TOLERANCE;
        usePolakRibiere = DEFAULT_USE_POLAK_RIBIERE;
        gradientListener = null;
        iter = 0;
    }

    /**
     * Constructor.
     *
     * @param listener         Listener to evaluate a multidimension function.
     * @param gradientListener Listener to obtain gradient value for the
     *                         multidimension function being evaluated.
     * @param point            Start point where algorithm will be started. Start point
     *                         should be close to the local minimum to be found. Provided array must
     *                         have a length equal to the number of dimensions of the function being
     *                         evaluated, otherwise and exception will be raised when searching for the
     *                         minimum.
     * @param direction        Direction to start looking for a minimum. Provided array
     *                         must have the same length as the number of dimensions of the function
     *                         being evaluated. Provided direction is considered as a vector pointing
     *                         to the minimum to be found.
     * @param tolerance        Tolerance or accuracy to be expected on estimated local
     *                         minimum.
     * @param usePolakRibiere  True if Polak-Ribiere method is used, otherwise
     *                         Fletcher-Reeves will be used.
     * @throws IllegalArgumentException Raised if provided point and direction
     *                                  don't have the same length or if provided tolerance is negative.
     */
    public ConjugateGradientMultiOptimizer(
            final MultiDimensionFunctionEvaluatorListener listener,
            final GradientFunctionEvaluatorListener gradientListener, final double[] point,
            final double[] direction, final double tolerance, final boolean usePolakRibiere) {
        super(listener, point, direction);
        internalSetTolerance(tolerance);
        this.usePolakRibiere = usePolakRibiere;
        this.gradientListener = gradientListener;
        iter = 0;
    }

    /**
     * Constructor.
     *
     * @param listener         Listener to evaluate a multidimensional function.
     * @param gradientListener Listener to obtain gradient value for the
     *                         multidimensional function being evaluated.
     * @param point            Start point where algorithm will be started. Start point
     *                         should be close to the local minimum to be found. Provided array must
     *                         have a length equal to the number of dimensions of the function being
     *                         evaluated, otherwise and exception will be raised when searching for the
     *                         minimum.
     * @param tolerance        Tolerance or accuracy to be expected on estimated local
     *                         minimum.
     * @param usePolakRibiere  True if Polak-Ribiere method is used, otherwise
     *                         Fletcher-Reeves will be used.
     * @throws IllegalArgumentException Raised if tolerance is negative.
     */
    public ConjugateGradientMultiOptimizer(
            final MultiDimensionFunctionEvaluatorListener listener,
            final GradientFunctionEvaluatorListener gradientListener, final double[] point,
            final double tolerance, final boolean usePolakRibiere) {
        super(listener);
        internalSetStartPoint(point);
        internalSetTolerance(tolerance);
        this.usePolakRibiere = usePolakRibiere;
        this.gradientListener = gradientListener;
        iter = 0;
    }

    /**
     * This function estimates a function minimum.
     * Implementations of this class will usually search a local minimum close
     * to a start point and will start looking into provided start direction.
     * Minimization of a function f. Input consists of an initial starting point
     * p. The initial matrix xi, whose columns contain the initial set of
     * directions, is set to the identity. Returned is the best point found, at
     * which point fret is the minimum function value and iter is the number of
     * iterations taken.
     *
     * @throws LockedException       Raised if this instance is locked, because
     *                               estimation is being computed.
     * @throws NotReadyException     Raised if this instance is not ready, because
     *                               a listener, a gradient listener and a start point haven't been provided.
     * @throws OptimizationException Raised if the algorithm failed because of
     *                               lack of convergence or because function couldn't be evaluated.
     */
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

        final int n = p.length;

        //set vector of directions
        if (!isDirectionAvailable()) {
            xi = new double[n];
        } else {
            if (xi.length != n) {
                xi = new double[n];
            }
        }

        boolean validResult = false;
        try {
            double gg;
            double dgg;

            final double[] g = new double[n];
            final double[] h = new double[n];

            double fp = listener.evaluate(p);
            gradientListener.evaluateGradient(p, xi);

            for (int j = 0; j < n; j++) {
                g[j] = -xi[j];
                h[j] = g[j];
                xi[j] = h[j];
            }
            for (int its = 0; its < ITMAX; its++) {
                iter = its;
                fret = linmin();
                if (2.0 * Math.abs(fret - fp) <= tolerance * (Math.abs(fret) +
                        Math.abs(fp) + EPS)) {
                    //minimum found
                    validResult = true;

                    if (iterationCompletedListener != null) {
                        iterationCompletedListener.onIterationCompleted(this, its, ITMAX);
                    }
                    break;
                }

                fp = fret;

                gradientListener.evaluateGradient(p, xi);

                double test = 0.0;
                final double den = Math.max(Math.abs(fp), 1.0);
                for (int j = 0; j < n; j++) {
                    final double temp = Math.abs(xi[j]) *
                            Math.max(Math.abs(p[j]), 1.0) / den;

                    if (temp > test) {
                        test = temp;
                    }
                }
                if (test < GTOL) {
                    // minimum found
                    validResult = true;

                    if (iterationCompletedListener != null) {
                        iterationCompletedListener.onIterationCompleted(this, its, ITMAX);
                    }
                    break;
                }

                dgg = gg = 0.0;
                for (int j = 0; j < n; j++) {
                    gg += g[j] * g[j];

                    if (isPolakRibiereEnabled()) {
                        // This statement for Polak-Ribiere
                        dgg += (xi[j] + g[j]) * xi[j];
                    } else {
                        // This statement for Fletcher-Reeves
                        dgg += xi[j] * xi[j];
                    }
                }

                if (gg == 0.0) {
                    // minimum found
                    validResult = true;

                    if (iterationCompletedListener != null) {
                        iterationCompletedListener.onIterationCompleted(this, its, ITMAX);
                    }
                    break;
                }

                double gam = dgg / gg;
                for (int j = 0; j < n; j++) {
                    g[j] = -xi[j];
                    h[j] = g[j] + gam * h[j];
                    xi[j] = h[j];
                }

                if (iterationCompletedListener != null) {
                    iterationCompletedListener.onIterationCompleted(this, its, ITMAX);
                }
            }

            if (!validResult) {
                //too many iterations
                locked = false;
                throw new OptimizationException();
            }
        } catch (final EvaluationException e) {
            throw new OptimizationException(e);
        } finally {
            locked = false;
        }

        // set result
        xmin = p;
        resultAvailable = true;
        fmin = fret;
    }

    /**
     * Returns boolean indicating whether this instance is ready to start the
     * estimation of a local minimum.
     * An instance is ready once a listener, a gradient listener and a start
     * point are provided.
     *
     * @return True if this instance is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isGradientListenerAvailable() &&
                isStartPointAvailable();
    }

    /**
     * Returns tolerance or accuracy to be expected on estimated local minimum.
     *
     * @return Tolerance or accuracy to be expected on estimated local minimum.
     */
    public double getTolerance() {
        return tolerance;
    }

    /**
     * Sets tolerance or accuracy to be expected on estimated local minimum.
     *
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     *                  minimum.
     * @throws LockedException          Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided tolerance is
     *                                  negative.
     */
    public void setTolerance(final double tolerance) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetTolerance(tolerance);
    }

    /**
     * Returns gradient listener in charge of obtaining gradient values for the
     * function to be evaluated.
     *
     * @return Gradient listener.
     * @throws NotAvailableException Raised if gradient listener has not yet
     *                               been provided.
     */
    public GradientFunctionEvaluatorListener getGradientListener()
            throws NotAvailableException {
        if (!isGradientListenerAvailable()) {
            throw new NotAvailableException();
        }
        return gradientListener;
    }

    /**
     * Sets gradient listener in charge of obtaining gradient values for the
     * function to be evaluated.
     *
     * @param gradientListener Gradient listener.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setGradientListener(
            final GradientFunctionEvaluatorListener gradientListener)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }

        this.gradientListener = gradientListener;
    }

    /**
     * Returns boolean indicating whether a gradient listener has already been
     * provided and is available for retrieval.
     *
     * @return True if available, false otherwise.
     */
    public boolean isGradientListenerAvailable() {
        return gradientListener != null;
    }

    /**
     * Returns boolean indicating whether Polak-Ribiere method is used or
     * Fletcher-Reeves is used instead.
     *
     * @return If true, Polak-Ribiere method is used, otherwise Fletcher-Reeves
     * is used.
     */
    public boolean isPolakRibiereEnabled() {
        return usePolakRibiere;
    }

    /**
     * Sets boolean indicating whether Polak-Ribiere method or Fletcher-Reeves
     * method is used.
     * If provided value is true, Polak-Ribiere method will be used, otherwise
     * Flecther-Reeves method will be used.
     *
     * @param useIt Boolean to determine method.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setUsePolakRibiere(final boolean useIt) throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        usePolakRibiere = useIt;
    }

    /**
     * Sets start point where local minimum is searched nearby.
     *
     * @param point Start point to search for a local minimum.
     * @throws LockedException Raised if this instance is locked.
     */
    public void setStartPoint(final double[] point)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetStartPoint(point);
    }

    /**
     * Return number of iterations that were needed to estimate a minimum.
     *
     * @return number of iterations that were needed.
     */
    public int getIterations() {
        return iter;
    }

    /**
     * Internal method to set tolerance or accuracy to be expected on estimated
     * local minimum.
     * This method does not check whether this instance is locked.
     *
     * @param tolerance Tolerance or accuracy to be expected on estimated local
     *                  minimum.
     * @throws IllegalArgumentException Raised if provided tolerance is
     *                                  negative.
     */
    private void internalSetTolerance(final double tolerance) {
        if (tolerance < MIN_TOLERANCE) {
            throw new IllegalArgumentException();
        }
        this.tolerance = tolerance;
    }

    /**
     * Internal method to set start point where local minimum is searched
     * nearby.
     * This method does not check whether this instance is locked.
     *
     * @param point Start point to search for a local minimum.
     */
    private void internalSetStartPoint(final double[] point) {
        p = point;
    }
}
