/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;

/**
 * Class to estimate the Jacobian of a multi variate and multidimensional
 * function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the Jacobian at such point.
 * <p>
 * For each input parameter {@code x}, the finite difference step is computed as
 * {@code h = max(relativeStep * |x|, minStep)}. If the resulting step is exactly zero
 * (e.g. when {@code x == 0.0} and {@code minStep == 0.0}), {@code relativeStep} is used as
 * step. The step is finally adjusted to the value that is exactly representable
 * as the difference between the perturbed and the original parameter.
 * <p>
 * By default {@code relativeStep} is {@link #EPS} and {@code minStep} is {@code 0.0},
 * which preserves the behaviour of a purely relative step.
 * Setting an absolute minimum step assumes that parameters are scaled to order 1,
 * since the same absolute step is used regardless of the magnitude of the parameter's
 * natural scale.
 * <p>
 * Instances of this class reuse internal work arrays between calls to avoid
 * allocations, hence they are not thread-safe nor re-entrant: a single instance
 * must not be used concurrently from several threads, nor from within the
 * evaluated function itself.
 */
public class JacobianEstimator {

    /**
     * Constant considered as machine precision.
     * This is the default relative step used to estimate the Jacobian.
     */
    public static final double EPS = 1e-8;

    /**
     * Default minimum absolute step.
     */
    private static final double DEFAULT_MIN_STEP = 0.0;

    /**
     * Listener to evaluate a multivariate function.
     */
    private final MultiVariateFunctionEvaluatorListener listener;

    /**
     * Relative step respect to the magnitude of each parameter.
     */
    private double relativeStep = EPS;

    /**
     * Minimum absolute step.
     */
    private double minStep = DEFAULT_MIN_STEP;

    /**
     * Internal array to hold input parameter's values.
     */
    private double[] xh;

    /**
     * Internal array to hold function values at the evaluated point.
     */
    private double[] fold;

    /**
     * Internal array to hold function values at the perturbed point.
     */
    private double[] fh;

    /**
     * Constructor.
     * Uses {@link #EPS} as relative step and {@code 0.0} as minimum step.
     *
     * @param listener listener to evaluate a multivariate function.
     */
    public JacobianEstimator(final MultiVariateFunctionEvaluatorListener listener) {
        this(listener, EPS, DEFAULT_MIN_STEP);
    }

    /**
     * Constructor.
     *
     * @param listener     listener to evaluate a multivariate function.
     * @param relativeStep relative step respect to the magnitude of each parameter. Must be
     *                     strictly positive.
     * @param minStep      minimum absolute step. Must be positive or zero.
     * @throws IllegalArgumentException if {@code relativeStep} is not strictly positive (or is
     *                                  NaN or infinite), or {@code minStep} is negative (or is NaN or infinite).
     */
    public JacobianEstimator(
            final MultiVariateFunctionEvaluatorListener listener, final double relativeStep, final double minStep) {
        this.listener = listener;
        this.relativeStep = checkRelativeStep(relativeStep);
        this.minStep = checkMinStep(minStep);
    }

    /**
     * Gets relative step respect to the magnitude of each parameter.
     * Default value is {@link #EPS}.
     *
     * @return relative step.
     */
    public double getRelativeStep() {
        return relativeStep;
    }

    /**
     * Sets relative step respect to the magnitude of each parameter.
     *
     * @param relativeStep relative step. Must be strictly positive.
     * @throws IllegalArgumentException if provided value is not strictly positive, or is NaN or infinite.
     */
    public void setRelativeStep(final double relativeStep) {
        this.relativeStep = checkRelativeStep(relativeStep);
    }

    /**
     * Gets minimum absolute step.
     * Default value is {@code 0.0}.
     * An absolute minimum step assumes that parameters are scaled to order 1.
     *
     * @return minimum absolute step.
     */
    public double getMinStep() {
        return minStep;
    }

    /**
     * Sets minimum absolute step.
     * An absolute minimum step assumes that parameters are scaled to order 1.
     *
     * @param minStep minimum absolute step. Must be positive or zero.
     * @throws IllegalArgumentException if provided value is negative, NaN or infinite.
     */
    public void setMinStep(final double minStep) {
        this.minStep = checkMinStep(minStep);
    }

    /**
     * Returns the Jacobian of a multivariate function at provided point.
     *
     * @param point input point.
     * @return jacobian.
     * @throws EvaluationException raised if function cannot be evaluated.
     */
    public Matrix jacobian(final double[] point) throws EvaluationException {
        try {
            final var result = new Matrix(listener.getNumberOfVariables(), point.length);
            jacobian(point, result);
            return result;
        } catch (final WrongSizeException e) {
            throw new EvaluationException(e);
        }
    }

    /**
     * Sets estimated jacobian in provided result matrix of a multivariate
     * function at provided point.
     * This method is preferred respect to jacobian(double[]) because result
     * matrix can be reused and hence is more memory efficient.
     * The step used for each parameter {@code x} is
     * {@code max(relativeStep * |x|, minStep)}, or {@code relativeStep} if such value is zero.
     *
     * @param point  input point.
     * @param result output matrix containing estimated jacobian. This parameter
     *               must be a matrix having the same number of columns as the length of the
     *               point and the same number of rows as the number of variables returned by
     *               function evaluations.
     * @throws EvaluationException      raised if function cannot be evaluated.
     * @throws IllegalArgumentException if size of result is not valid.
     */
    public void jacobian(final double[] point, final Matrix result) throws EvaluationException {
        final var numdims = point.length;
        final var numvars = listener.getNumberOfVariables();
        if (result.getColumns() != numdims) {
            throw new IllegalArgumentException();
        }
        if (result.getRows() != numvars) {
            throw new IllegalArgumentException();
        }

        if (fold == null || fold.length != numvars) {
            fold = new double[numvars];
            fh = new double[numvars];
        }
        if (xh == null || xh.length != numdims) {
            xh = new double[numdims];
        }
        System.arraycopy(point, 0, xh, 0, numdims);

        double temp;
        double h;
        listener.evaluate(point, fold);
        for (var j = 0; j < numdims; j++) {
            temp = point[j];
            h = computeStep(temp);
            xh[j] = temp + h;
            h = xh[j] - temp;
            listener.evaluate(xh, fh);
            xh[j] = temp;
            for (var i = 0; i < numvars; i++) {
                result.setElementAt(i, j, (fh[i] - fold[i]) / h);
            }
        }
    }

    /**
     * Computes the step used to perturb an input parameter.
     * The step is {@code max(relativeStep * |x|, minStep)}, and if such value is zero, then
     * {@code relativeStep} is used instead.
     * Notice that the actually applied step might differ slightly because of finite machine precision.
     *
     * @param x value of the input parameter to be perturbed.
     * @return step to be used.
     */
    protected double computeStep(final double x) {
        final var h = Math.max(relativeStep * Math.abs(x), minStep);
        // Trick to reduce finite-precision error
        return h == 0.0 ? relativeStep : h;
    }

    /**
     * Validates a relative step.
     *
     * @param relativeStep value to be checked.
     * @return provided value if valid.
     * @throws IllegalArgumentException if value is not strictly positive, or is NaN or infinite.
     */
    private static double checkRelativeStep(final double relativeStep) {
        if (!(relativeStep > 0.0) || Double.isInfinite(relativeStep)) {
            throw new IllegalArgumentException("relativeStep must be strictly positive and finite");
        }
        return relativeStep;
    }

    /**
     * Validates a minimum step.
     *
     * @param minStep value to be checked.
     * @return provided value if valid.
     * @throws IllegalArgumentException if value is negative, NaN or infinite.
     */
    private static double checkMinStep(final double minStep) {
        if (!(minStep >= 0.0) || Double.isInfinite(minStep)) {
            throw new IllegalArgumentException("minStep must be positive or zero and finite");
        }
        return minStep;
    }
}
