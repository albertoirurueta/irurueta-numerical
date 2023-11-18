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

import com.irurueta.numerical.DirectionalEvaluator;
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NumericalException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Abstract class to search for a local minimum on a multidimensional function
 * along a given line of input parameters.
 * Line minimization implementations exist such as Powell's method or conjugate
 * gradients. Among those, usually conjugate gradients provide faster
 * convergence because search directions during the estimation don't need to be
 * perpendicular.
 */
public abstract class LineMultiOptimizer extends MultiOptimizer {
    /**
     * n-dimensional point containing a minimum in a given line.
     */
    protected double[] p;

    /**
     * Direction to make the search.
     */
    protected double[] xi;

    /**
     * Number of dimensions on function being evaluated.
     */
    private int n;

    /**
     * Class in charge of evaluating a function through a given line.
     */
    private DirectionalEvaluator evaluator;

    /**
     * Internal optimizer to find a minimum of a function along a line of
     * input values. Hence, input is converted to a single dimension using a
     * DirectionalEvaluator.
     */
    private BrentSingleOptimizer brent;

    /**
     * Empty constructor.
     */
    protected LineMultiOptimizer() {
        super();
        p = xi = null;
        n = 0;
        evaluator = null;
        brent = null;
    }

    /**
     * Constructor.
     *
     * @param listener Listener to evaluate a multi-dimension function.
     */
    protected LineMultiOptimizer(
            final MultiDimensionFunctionEvaluatorListener listener) {
        super(listener);
        p = xi = null;
        n = 0;
        evaluator = null;
        brent = null;
    }

    /**
     * Constructor.
     *
     * @param listener  Listener to evaluate a multi-dimension function.
     * @param point     Start point where algorithm will be started. Start point
     *                  should be close to the local minimum to be found. Provided array must
     *                  have a length equal to the number of dimensions of the function being
     *                  evaluated, otherwise and exception will be raised when searching for the
     *                  minimum.
     * @param direction Direction to start looking for a minimum. Provided array
     *                  must have the same length as the number of dimensions of the function
     *                  being evaluated. Provided direction is considered as a vector pointing
     *                  to the minimum to be found.
     * @throws IllegalArgumentException Raised if provided point and direction
     *                                  don't have the same length.
     */
    protected LineMultiOptimizer(final MultiDimensionFunctionEvaluatorListener listener,
                                 final double[] point, final double[] direction) {
        super(listener);
        internalSetStartPointAndDirection(point, direction);
        n = 0;
        evaluator = null;
        brent = null;
    }

    /**
     * Returns boolean indicating whether start point has already been provided
     * and is ready for retrieval.
     *
     * @return True if available, false otherwise.
     */
    public boolean isStartPointAvailable() {
        return p != null;
    }

    /**
     * Returns start point where algorithm will be started. Start point should
     * be close to the local minimum to be found.
     *
     * @return Start point where algorithm will be started.
     * @throws NotAvailableException Raised if start point has not yet been
     *                               provided and is not available.
     */
    public double[] getStartPoint() throws NotAvailableException {
        if (!isStartPointAvailable()) {
            throw new NotAvailableException();
        }
        return p;
    }

    /**
     * Returns boolean indicating whether direction has already been provided
     * and is ready for retrieval.
     *
     * @return True if available, false otherwise.
     */
    public boolean isDirectionAvailable() {
        return xi != null;
    }

    /**
     * Returns direction to start looking for a minimum. Provided array must
     * have the same length as the number of dimensions of the function being
     * evaluated. Provided direction is considered as a vector pointing to the
     * minimum to be found.
     *
     * @return Direction to start looking for a minimum.
     * @throws NotAvailableException Raised if direction has not yet been
     *                               provided and is not available.
     */
    public double[] getDirection() throws NotAvailableException {
        if (!isDirectionAvailable()) {
            throw new NotAvailableException();
        }
        return xi;
    }

    /**
     * Internal method to set start point and direction to start the search for
     * a local minimum.
     *
     * @param point     Start point where algorithm will be started. Start point
     *                  should be close to the local minimum to be found. Provided array must
     *                  have a length equal to the number of dimensions of the function being
     *                  evaluated, otherwise and exception will be raised when searching for the
     *                  minimum.
     * @param direction Direction to start looking for a minimum. Provided array
     *                  must have the same length as the number of dimensions of the function
     *                  being evaluated. Provided direction is considered as a vector pointing
     *                  to the minimum to be found.
     * @throws LockedException          Raised if this instance is locked.
     * @throws IllegalArgumentException Raised if provided point and direction
     *                                  don't have the same length.
     */
    public void setStartPointAndDirection(final double[] point, final double[] direction)
            throws LockedException {
        if (isLocked()) {
            throw new LockedException();
        }
        internalSetStartPointAndDirection(point, direction);
    }

    /**
     * Returns boolean indicating whether this instance is considered to be
     * ready to start the estimation of a minimum.
     * This instance is considered to be ready once a listener, start point and
     * direction are provided.
     *
     * @return True if this instance is ready, false otherwise.
     */
    @Override
    public boolean isReady() {
        return isListenerAvailable() && isStartPointAvailable() &&
                isDirectionAvailable();
    }

    /**
     * Searches for a minimum along a given line of input values.
     * The line being searched is obtained by using a start point and direction.
     *
     * @return Returns function evaluation at minimum that has been found.
     */
    @SuppressWarnings("Duplicates")
    protected double linmin() {
        final double ax;
        final double xx;
        final double linxmin;
        n = p.length;

        if (evaluator == null) {
            // attempt to reuse evaluator
            evaluator = new DirectionalEvaluator(listener, p, xi);
        }
        if (evaluator.getListener() != listener) {
            // update listener
            evaluator.setListener(listener);
        }
        if (evaluator.getPoint() != p || evaluator.getDirection() != xi) {
            evaluator.setPointAndDirection(p, xi);
        }

        ax = 0.0;
        xx = 1.0;

        try {
            if (brent == null) {
                // attempt to reuse brent single optimizer
                brent = new BrentSingleOptimizer(
                        new SingleDimensionFunctionEvaluatorListener() {

                            @Override
                            public double evaluate(double point) throws EvaluationException {
                                return evaluator.evaluateAt(point);
                            }
                        }, BrentSingleOptimizer.DEFAULT_MIN_EVAL_POINT,
                        BrentSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT,
                        BrentSingleOptimizer.DEFAULT_MAX_EVAL_POINT,
                        BrentSingleOptimizer.DEFAULT_TOLERANCE);
            }

            brent.computeBracket(ax, xx);
            brent.minimize();
            linxmin = brent.getResult();

            for (int j = 0; j < n; j++) {
                xi[j] *= linxmin;
                p[j] += xi[j];
            }

            return brent.getEvaluationAtResult();
        } catch (final NumericalException e) {
            // if minimization fails, try to evaluate at best point found so far
            try {
                return listener.evaluate(p);
            } catch (EvaluationException e2) {
                // if minimization fails here we assume that obtained result is
                // the worst possible one
                return Double.MAX_VALUE;
            }
        }
    }

    /**
     * Internal method to set start point and direction to start the search for
     * a local minimum.
     * This method does not check whether this instance is locked.
     *
     * @param point     Start point where algorithm will be started. Start point
     *                  should be close to the local minimum to be found. Provided array must
     *                  have a length equal to the number of dimensions of the function being
     *                  evaluated, otherwise and exception will be raised when searching for the
     *                  minimum.
     * @param direction Direction to start looking for a minimum. Provided array
     *                  must have the same length as the number of dimensions of the function
     *                  being evaluated. Provided direction is considered as a vector pointing
     *                  to the minimum to be found.
     * @throws IllegalArgumentException Raised if provided point and direction
     *                                  don't have the same length.
     */
    private void internalSetStartPointAndDirection(final double[] point,
                                                   final double[] direction) {
        if (point.length != direction.length) {
            throw new IllegalArgumentException();
        }
        p = point;
        xi = direction;
    }
}
