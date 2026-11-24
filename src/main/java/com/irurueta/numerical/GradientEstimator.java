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
 * Class to estimate the gradient of a multidimensional function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the gradient at such point.
 * The step used to perturb each input parameter x is computed as
 * {@code max(relativeStep * |x|, minStep)}, falling back to {@code relativeStep} if the result is zero.
 * By default, {@code relativeStep = EPS} and {@code minStep = 0.0}, which keeps the classical behavior
 * ({@code EPS * |x|}). For small non-zero parameters such step collapses and precision is lost, in which case
 * a positive minimum step can be configured (e.g. {@code relativeStep = 1e-7} and {@code minStep = 1e-7} produce
 * steps of {@code 1e-7 * max(1, |x|)}).
 * Notice that the default {@code relativeStep} is tuned for first order (forward difference) estimation, whose
 * optimal step scales with the square root of machine precision. Second order estimators (symmetric or
 * Savitzky-Golay) usually have an optimal relative step that scales with the cubic root of machine precision
 * (around 1e-5), so when using them larger relative steps might be more appropriate.
 */
public class GradientEstimator {

    /**
     * Constant considered as machine precision.
     */
    public static final double EPS = 1e-8;

    /**
     * Listener to evaluate a multidimensional function.
     */
    protected final MultiDimensionFunctionEvaluatorListener listener;

    /**
     * Internal array to hold input parameter's values.
     */
    private double[] xh;

    /**
     * Relative step respect to the absolute value of each input parameter.
     */
    private double relativeStep = EPS;

    /**
     * Minimum step used to perturb each input parameter.
     */
    private double minStep;

    /**
     * Constructor.
     *
     * @param listener Listener to evaluate a multidimensional function.
     */
    public GradientEstimator(final MultiDimensionFunctionEvaluatorListener listener) {
        this.listener = listener;
    }

    /**
     * Constructor.
     *
     * @param listener     Listener to evaluate a multidimensional function.
     * @param relativeStep Relative step respect to the absolute value of each input parameter. Must be
     *                     strictly greater than zero.
     * @param minStep      Minimum step used to perturb each input parameter. Must be zero or positive.
     * @throws IllegalArgumentException Raised if relative step is not strictly positive or if minimum step is
     *                                  negative or NaN.
     */
    public GradientEstimator(
            final MultiDimensionFunctionEvaluatorListener listener, final double relativeStep,
            final double minStep) {
        this.listener = listener;
        setRelativeStep(relativeStep);
        setMinStep(minStep);
    }

    /**
     * Gets relative step respect to the absolute value of each input parameter.
     * Default value is {@link #EPS}.
     *
     * @return relative step.
     */
    public double getRelativeStep() {
        return relativeStep;
    }

    /**
     * Sets relative step respect to the absolute value of each input parameter.
     *
     * @param relativeStep relative step. Must be strictly greater than zero.
     * @throws IllegalArgumentException Raised if provided value is zero, negative or NaN.
     */
    public void setRelativeStep(final double relativeStep) {
        if (!(relativeStep > 0.0)) {
            throw new IllegalArgumentException("relativeStep must be strictly positive");
        }
        this.relativeStep = relativeStep;
    }

    /**
     * Gets minimum step used to perturb each input parameter.
     * Default value is zero, which means that only the relative step is used.
     *
     * @return minimum step.
     */
    public double getMinStep() {
        return minStep;
    }

    /**
     * Sets minimum step used to perturb each input parameter.
     *
     * @param minStep minimum step. Must be zero or positive.
     * @throws IllegalArgumentException Raised if provided value is negative or NaN.
     */
    public void setMinStep(final double minStep) {
        if (!(minStep >= 0.0)) {
            throw new IllegalArgumentException("minStep must be zero or positive");
        }
        this.minStep = minStep;
    }

    /**
     * Returns the gradient of a multidimensional function at provided point.
     *
     * @param point Input point.
     * @return Gradient.
     * @throws EvaluationException Raised if function cannot be evaluated.
     */
    public double[] gradient(final double[] point) throws EvaluationException {
        final var result = new double[point.length];
        gradient(point, result);
        return result;
    }

    /**
     * Sets estimated gradient in provided result array of a multidimensional
     * function at provided point.
     * This method is preferred respect to gradient(double[]) because result
     * array can be reused and hence is more memory efficient.
     *
     * @param point  Input point.
     * @param result Output parameter containing estimated array. This parameter
     *               must be an array of length equal to point.
     * @throws EvaluationException      Raised if function cannot be evaluated.
     * @throws IllegalArgumentException Raised if length of result and point are
     *                                  not equal.
     */
    public void gradient(final double[] point, final double[] result) throws EvaluationException {
        final var length = point.length;
        if (result.length != length) {
            throw new IllegalArgumentException();
        }

        if (xh == null || xh.length != length) {
            xh = new double[length];
        }
        System.arraycopy(point, 0, xh, 0, length);

        double temp;
        double h;
        double fh;
        final var fold = listener.evaluate(point);
        for (var j = 0; j < length; j++) {
            temp = point[j];
            h = computeStep(temp);
            xh[j] = temp + h;
            h = xh[j] - temp;
            fh = listener.evaluate(xh);
            xh[j] = temp;
            result[j] = (fh - fold) / h;
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
}
