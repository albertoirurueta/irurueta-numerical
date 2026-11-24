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
package com.irurueta.numerical;

/**
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtained.
 * This implementation is faster although less accurate than
 * SymmetricDerivativeEstimator.
 * The step used to perturb the evaluation point is computed as
 * {@code max(relativeStep * |x|, minStep)}. By default the relative step is
 * {@link #EPS} and the minimum step is zero, which preserves the historical
 * behaviour. Setting a minimum step prevents catastrophic loss of precision
 * when evaluating derivatives at points close to zero. Notice that an
 * absolute minimum step assumes that parameters are scaled to order 1.
 */
public class DerivativeEstimator {

    /**
     * Constant defining machine precision for this algorithm.
     */
    public static final double EPS = 1e-8;

    /**
     * Listener to evaluate a single dimension function.
     */
    protected final SingleDimensionFunctionEvaluatorListener listener;

    /**
     * Relative step, which is multiplied by the absolute value of the
     * evaluation point to obtain the step used to perturb such point.
     */
    private double relativeStep;

    /**
     * Minimum absolute step used to perturb the evaluation point.
     */
    private double minStep;

    /**
     * Constructor using default relative step ({@link #EPS}) and no minimum
     * step.
     *
     * @param listener listener to evaluate a single dimension function
     */
    public DerivativeEstimator(final SingleDimensionFunctionEvaluatorListener listener) {
        this(listener, EPS, 0.0);
    }

    /**
     * Constructor.
     *
     * @param listener     listener to evaluate a single dimension function
     * @param relativeStep relative step, which is multiplied by the absolute
     *                     value of the evaluation point. Must be strictly
     *                     positive.
     * @param minStep      minimum absolute step. Must be zero or positive.
     *                     Notice that an absolute minimum step assumes that
     *                     parameters are scaled to order 1.
     * @throws IllegalArgumentException if relative step is not strictly
     *                                  positive, or minimum step is negative
     *                                  or NaN or infinite
     */
    public DerivativeEstimator(
            final SingleDimensionFunctionEvaluatorListener listener, final double relativeStep,
            final double minStep) {
        this.listener = listener;
        this.relativeStep = checkRelativeStep(relativeStep);
        this.minStep = checkMinStep(minStep);
    }

    /**
     * Gets relative step, which is multiplied by the absolute value of the
     * evaluation point to obtain the perturbation step.
     *
     * @return relative step
     */
    public double getRelativeStep() {
        return relativeStep;
    }

    /**
     * Sets relative step, which is multiplied by the absolute value of the
     * evaluation point to obtain the perturbation step.
     *
     * @param relativeStep relative step. Must be strictly positive.
     * @throws IllegalArgumentException if provided value is not strictly
     *                                  positive, or is NaN or infinite
     */
    public void setRelativeStep(final double relativeStep) {
        this.relativeStep = checkRelativeStep(relativeStep);
    }

    /**
     * Gets minimum absolute step used to perturb the evaluation point.
     *
     * @return minimum step
     */
    public double getMinStep() {
        return minStep;
    }

    /**
     * Sets minimum absolute step used to perturb the evaluation point.
     * Notice that an absolute minimum step assumes that parameters are scaled
     * to order 1.
     *
     * @param minStep minimum step. Must be zero or positive.
     * @throws IllegalArgumentException if provided value is negative, NaN or infinite
     */
    public void setMinStep(final double minStep) {
        this.minStep = checkMinStep(minStep);
    }

    /**
     * Computes the function derivative at provided point x.
     *
     * @param x Point where derivative is estimated
     * @return Derivative of function at provided point
     * @throws EvaluationException Raised if function cannot be properly
     *                             evaluated
     */
    public double derivative(final double x) throws EvaluationException {
        final var fold = listener.evaluate(x);

        final var xh = x + computeStep(x);
        final var h = xh - x;

        final var fh = listener.evaluate(xh);
        return (fh - fold) / h;
    }

    /**
     * Computes the step used to perturb the evaluation point, which is
     * {@code max(relativeStep * |x|, minStep)}. If such value is zero, the
     * relative step is used instead to reduce finite-precision error.
     *
     * @param x point where derivative is estimated
     * @return step to perturb the evaluation point
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
            throw new IllegalArgumentException("relative step must be strictly positive and finite");
        }
        return relativeStep;
    }

    /**
     * Validates a minimum step.
     *
     * @param minStep value to be checked.
     * @return provided value if valid.
     * @throws IllegalArgumentException if value is negative, or is NaN or infinite.
     */
    private static double checkMinStep(final double minStep) {
        if (!(minStep >= 0.0) || Double.isInfinite(minStep)) {
            throw new IllegalArgumentException("minimum step must be zero or positive and finite");
        }
        return minStep;
    }
}
