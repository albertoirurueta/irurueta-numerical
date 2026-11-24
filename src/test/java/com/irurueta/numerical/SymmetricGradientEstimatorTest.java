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

import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SymmetricGradientEstimatorTest implements MultiDimensionFunctionEvaluatorListener {

    private static final int MIN_DIMS = 1;
    private static final int MAX_DIMS = 3;

    private static final double MIN_EVAL_POINT = -10.0;
    private static final double MAX_EVAL_POINT = 10.0;

    private static final double MIN_OFFSET = -10.0;
    private static final double MAX_OFFSET = 10.0;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;

    private static final double ABSOLUTE_ERROR = 1e-1;

    private static final int TIMES = 100;

    /**
     * Relative step and minimum step used to obtain steps equal to 1e-7 * max(1, |x|).
     */
    private static final double CUSTOM_STEP = 1e-7;

    /**
     * Machine precision constant used by default as relative step.
     */
    private static final double DEFAULT_EPS = GradientEstimator.EPS;

    /**
     * Position, among the points evaluated by the estimator, of the first perturbed one.
     */
    private static final int PERTURBED_INDEX = 0;

    /**
     * Relative tolerance when comparing expected steps with the steps actually applied.
     */
    private static final double RELATIVE_STEP_ERROR = 1e-6;

    /**
     * Small non-zero parameter values where default steps collapse.
     */
    private static final double[] SMALL_VALUES = {1e-5, 1e-8, 1e-10};

    /**
     * Parameter values used to check computed steps.
     */
    private static final double[] STEP_VALUES = {0.0, 1e-9, 1e-3, 1.0, 50.0, -50.0};

    /**
     * Maximum error expected when minimum step is configured.
     */
    private static final double SMALL_VALUE_ERROR = 1e-6;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    @Test
    void testConstructor() {

        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        minimum = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var point = new double[ndims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

        final var estimator = new SymmetricGradientEstimator(this);
        assertNotNull(estimator);
    }

    @Test
    void testGradient() throws EvaluationException {

        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            final var point = new double[ndims];
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

            final var estimator = new SymmetricGradientEstimator(this);

            final var gradient1 = estimator.gradient(point);
            final var gradient2 = new double[ndims];
            estimator.gradient(point, gradient2);

            // check correctness
            final var gradient3 = gradient(point);

            assertEquals(gradient1.length, ndims);
            assertEquals(gradient2.length, ndims);
            assertEquals(gradient3.length, ndims);
            for (var i = 0; i < ndims; i++) {
                assertEquals(gradient1[i], gradient3[i], ABSOLUTE_ERROR);
                assertEquals(gradient2[i], gradient3[i], ABSOLUTE_ERROR);
            }
        }
    }

        /**
     * Checks default settings.
     */
    @Test
    void constructor_whenDefault_hasDefaultSettings() {
        final var estimator = new SymmetricGradientEstimator(this);

        assertEquals(DEFAULT_EPS, estimator.getRelativeStep(), 0.0);
        assertEquals(0.0, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that settings are stored by the constructor.
     */
    @Test
    void constructor_whenCustomSettings_setsValues() {
        final var estimator = new SymmetricGradientEstimator(this, CUSTOM_STEP, 2.0 * CUSTOM_STEP);

        assertEquals(CUSTOM_STEP, estimator.getRelativeStep(), 0.0);
        assertEquals(2.0 * CUSTOM_STEP, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that an invalid relative step is rejected by the constructor.
     */
    @Test
    void constructor_whenInvalidRelativeStep_throwsIllegalArgumentException() {
        assertThrows(IllegalArgumentException.class, () -> new SymmetricGradientEstimator(this, 0.0, 0.0));
        assertThrows(IllegalArgumentException.class, () -> new SymmetricGradientEstimator(this, -CUSTOM_STEP, 0.0));
        assertThrows(IllegalArgumentException.class, () -> new SymmetricGradientEstimator(this, Double.NaN, 0.0));
    }

    /**
     * Checks that an invalid minimum step is rejected by the constructor.
     */
    @Test
    void constructor_whenInvalidMinStep_throwsIllegalArgumentException() {
        assertThrows(IllegalArgumentException.class,
                () -> new SymmetricGradientEstimator(this, CUSTOM_STEP, -CUSTOM_STEP));
        assertThrows(IllegalArgumentException.class,
                () -> new SymmetricGradientEstimator(this, CUSTOM_STEP, Double.NaN));
    }

    /**
     * Checks that relative step can be set and invalid values are rejected leaving the value unchanged.
     */
    @Test
    void setRelativeStep_whenValidOrInvalid_setsOrThrows() {
        final var estimator = new SymmetricGradientEstimator(this);

        estimator.setRelativeStep(CUSTOM_STEP);
        assertEquals(CUSTOM_STEP, estimator.getRelativeStep(), 0.0);

        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(0.0));
        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(-CUSTOM_STEP));
        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(Double.NaN));
        assertEquals(CUSTOM_STEP, estimator.getRelativeStep(), 0.0);
    }

    /**
     * Checks that minimum step can be set and invalid values are rejected leaving the value unchanged.
     */
    @Test
    void setMinStep_whenValidOrInvalid_setsOrThrows() {
        final var estimator = new SymmetricGradientEstimator(this);

        estimator.setMinStep(CUSTOM_STEP);
        assertEquals(CUSTOM_STEP, estimator.getMinStep(), 0.0);
        estimator.setMinStep(0.0);
        assertEquals(0.0, estimator.getMinStep(), 0.0);

        assertThrows(IllegalArgumentException.class, () -> estimator.setMinStep(-CUSTOM_STEP));
        assertThrows(IllegalArgumentException.class, () -> estimator.setMinStep(Double.NaN));
        assertEquals(0.0, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that with relative and minimum steps equal to 1e-7 the step is 1e-7 * max(1, |x|).
     *
     * @throws EvaluationException if function cannot be evaluated.
     */
    @Test
    void gradient_whenCustomSettings_usesRelativeStepWithMinimum() throws EvaluationException {
        for (final var x : STEP_VALUES) {
            final var points = evaluatedPoints(
                    listener -> new SymmetricGradientEstimator(listener, CUSTOM_STEP, CUSTOM_STEP), x);

            final var expected = CUSTOM_STEP * Math.max(1.0, Math.abs(x));
            final var applied = points.get(PERTURBED_INDEX) - x;
            assertEquals(expected, applied, expected * RELATIVE_STEP_ERROR);
        }
    }

    /**
     * Checks that with default settings the step is EPS * |x|, and EPS when x is zero.
     *
     * @throws EvaluationException if function cannot be evaluated.
     */
    @Test
    void gradient_whenDefaultSettings_usesLegacyStep() throws EvaluationException {
        for (final var x : STEP_VALUES) {
            final var points = evaluatedPoints(
                    listener -> new SymmetricGradientEstimator(listener), x);

            final var expected = x == 0.0 ? DEFAULT_EPS : DEFAULT_EPS * Math.abs(x);
            final var applied = points.get(PERTURBED_INDEX) - x;
            assertEquals(expected, applied, expected * RELATIVE_STEP_ERROR);
        }
    }

    /**
     * Checks that explicitly providing default settings gives identical results to the default constructor.
     *
     * @throws EvaluationException if function cannot be evaluated.
     */
    @Test
    void gradient_whenExplicitDefaultSettings_resultsAreIdentical() throws EvaluationException {
        final MultiDimensionFunctionEvaluatorListener listener = SymmetricGradientEstimatorTest::exp;
        final var defaultEstimator = new SymmetricGradientEstimator(listener);
        final var explicitEstimator = new SymmetricGradientEstimator(listener, DEFAULT_EPS, 0.0);

        for (final var x : STEP_VALUES) {
            final var point = new double[]{x};
            assertArrayEquals(defaultEstimator.gradient(point), explicitEstimator.gradient(point), 0.0);
        }
    }

    /**
     * Checks that a zero valued parameter is still perturbed, both with default and custom settings.
     *
     * @throws EvaluationException if function cannot be evaluated.
     */
    @Test
    void gradient_whenZeroParameter_isStillPerturbed() throws EvaluationException {
        final var defaultPoints = evaluatedPoints(listener -> new SymmetricGradientEstimator(listener), 0.0);
        assertNotEquals(0.0, defaultPoints.get(PERTURBED_INDEX), 0.0);

        final var customPoints = evaluatedPoints(
                listener -> new SymmetricGradientEstimator(listener, CUSTOM_STEP, 0.0), 0.0);
        assertNotEquals(0.0, customPoints.get(PERTURBED_INDEX), 0.0);
        assertEquals(CUSTOM_STEP, customPoints.get(PERTURBED_INDEX), CUSTOM_STEP * RELATIVE_STEP_ERROR);

        final var minPoints = evaluatedPoints(
                listener -> new SymmetricGradientEstimator(listener, DEFAULT_EPS, CUSTOM_STEP), 0.0);
        assertEquals(CUSTOM_STEP, minPoints.get(PERTURBED_INDEX), CUSTOM_STEP * RELATIVE_STEP_ERROR);
    }

    /**
     * Checks that for small non-zero parameters a minimum step is more accurate than default settings.
     *
     * @throws EvaluationException if function cannot be evaluated.
     */
    @Test
    void gradient_whenSmallParameterAndMinStep_isMoreAccurateThanDefault() throws EvaluationException {
        final MultiDimensionFunctionEvaluatorListener listener = SymmetricGradientEstimatorTest::exp;
        for (final var x : SMALL_VALUES) {
            final var defaultError = expError(new SymmetricGradientEstimator(listener), x);
            final var minStepError = expError(new SymmetricGradientEstimator(listener, CUSTOM_STEP, CUSTOM_STEP), x);

            assertTrue(minStepError < SMALL_VALUE_ERROR, "x=" + x + " error=" + minStepError);
            assertFalse(defaultError <= minStepError, "x=" + x + " default=" + defaultError);
        }
    }

    @Override
    public double evaluate(final double[] point) {
        final var dims = Math.min(Math.min(point.length, minimum.length), width.length);

        var value = 1.0;
        for (var i = 0; i < dims; i++) {
            value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
        }

        value += offset;
        return value;
    }

    private double[] gradient(final double[] params) {
        final var dims = Math.min(Math.min(params.length, minimum.length), width.length);

        final var gradient = new double[dims];

        for (var j = 0; j < dims; j++) {
            var value = 1.0;
            for (var i = 0; i < dims; i++) {
                if (i != j) {
                    value *= Math.pow(params[i] - minimum[i], 2.0) / width[i];
                } else {
                    value *= 2.0 * (params[i] - minimum[i]) / width[i];
                }
            }

            gradient[j] = value;
        }

        return gradient;
    }

    /**
     * Evaluates exp(x) on the first coordinate of provided point.
     *
     * @param point point to evaluate.
     * @return evaluated value.
     */
    private static double exp(final double[] point) {
        return Math.exp(point[0]);
    }

    /**
     * Estimates the gradient of exp(x) at provided value using provided estimator factory and returns
     * the absolute error respect to the analytical gradient.
     *
     * @param estimator estimator to be used.
     * @param x         point where gradient is estimated.
     * @return absolute error.
     * @throws EvaluationException if function cannot be evaluated.
     */
    private static double expError(final GradientEstimator estimator, final double x) throws EvaluationException {
        final var gradient = estimator.gradient(new double[]{x});
        return Math.abs(gradient[0] - Math.exp(x));
    }

    /**
     * Gets the first coordinate of the points where a function is evaluated when estimating the gradient
     * at provided value.
     *
     * @param factory factory creating the estimator from a recording listener.
     * @param x       point where gradient is estimated.
     * @return list containing the first coordinate of each evaluated point, in evaluation order.
     * @throws EvaluationException if function cannot be evaluated.
     */
    private static List<Double> evaluatedPoints(
            final Function<MultiDimensionFunctionEvaluatorListener, GradientEstimator> factory,
            final double x) throws EvaluationException {
        final var points = new ArrayList<Double>();
        final MultiDimensionFunctionEvaluatorListener recorder = point -> {
            points.add(point[0]);
            return exp(point);
        };
        factory.apply(recorder).gradient(new double[]{x});
        return points;
    }
}
