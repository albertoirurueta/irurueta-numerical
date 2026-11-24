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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DerivativeEstimatorTest implements SingleDimensionFunctionEvaluatorListener {

    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;

    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;

    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;

    public static final double ABSOLUTE_ERROR = 1e-2;

    /**
     * Relative step used in tests.
     */
    private static final double RELATIVE_STEP = 1e-7;

    /**
     * Minimum step used in tests.
     */
    private static final double MIN_STEP = 1e-7;

    /**
     * Absolute error allowed when evaluating at zero.
     */
    private static final double ZERO_ERROR = 1e-5;

    /**
     * Relative error allowed for the estimation using a minimum step at very small points.
     */
    private static final double CUSTOM_RELATIVE_ERROR = 1e-5;

    private double minimum;
    private double width;
    private double offset;

    @Test
    void testConstructor() {
        final var estimator = new DerivativeEstimator(this);
        assertNotNull(estimator);
    }

    @Test
    void testDerivative() throws EvaluationException {

        final var randomizer = new UniformRandomizer();
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var x = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        final var estimator = new DerivativeEstimator(this);

        // estimate derivative
        final var estDerivative = estimator.derivative(x);

        // real derivative
        final var realDerivative = derivative(x);

        // compare both results
        assertEquals(estDerivative, realDerivative, ABSOLUTE_ERROR);
    }

    /**
     * Checks default relative and minimum steps.
     */
    @Test
    void constructor_whenDefault_usesDefaultSteps() {
        final var estimator = new DerivativeEstimator(this);

        assertEquals(DerivativeEstimator.EPS, estimator.getRelativeStep(), 0.0);
        assertEquals(0.0, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks constructor with custom steps.
     */
    @Test
    void constructor_whenCustomSteps_setsSteps() {
        final var estimator = new DerivativeEstimator(this, RELATIVE_STEP, MIN_STEP);

        assertEquals(RELATIVE_STEP, estimator.getRelativeStep(), 0.0);
        assertEquals(MIN_STEP, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks constructor validation of invalid steps.
     */
    @Test
    void constructor_whenInvalidSteps_throwsIllegalArgumentException() {
        assertThrows(IllegalArgumentException.class, () -> new DerivativeEstimator(this, 0.0, MIN_STEP));
        assertThrows(IllegalArgumentException.class, () -> new DerivativeEstimator(this, -RELATIVE_STEP, MIN_STEP));
        assertThrows(IllegalArgumentException.class, () -> new DerivativeEstimator(this, Double.NaN, MIN_STEP));
        assertThrows(IllegalArgumentException.class, () -> new DerivativeEstimator(this, RELATIVE_STEP, -MIN_STEP));
        assertThrows(IllegalArgumentException.class, () -> new DerivativeEstimator(this, RELATIVE_STEP, Double.NaN));
    }

    /**
     * Checks setters with valid and invalid values.
     */
    @Test
    void setters_whenValidAndInvalidValues_setOrThrow() {
        final var estimator = new DerivativeEstimator(this);

        estimator.setRelativeStep(RELATIVE_STEP);
        estimator.setMinStep(MIN_STEP);
        assertEquals(RELATIVE_STEP, estimator.getRelativeStep(), 0.0);
        assertEquals(MIN_STEP, estimator.getMinStep(), 0.0);

        estimator.setMinStep(0.0);
        assertEquals(0.0, estimator.getMinStep(), 0.0);

        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(0.0));
        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(-1.0));
        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(Double.NaN));
        assertThrows(IllegalArgumentException.class, () -> estimator.setMinStep(-1.0));
        assertThrows(IllegalArgumentException.class, () -> estimator.setMinStep(Double.NaN));

        assertEquals(RELATIVE_STEP, estimator.getRelativeStep(), 0.0);
        assertEquals(0.0, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that step is {@code 1e-7 * max(1, |x|)} when both relative and minimum steps are 1e-7.
     */
    @Test
    void computeStep_whenRelativeAndMinStepAreEqual_returnsScaledStep() {
        final var estimator = new DerivativeEstimator(this, RELATIVE_STEP, MIN_STEP);

        for (final var x : new double[]{0.0, 1e-10, 1e-5, 0.5, 1.0, 10.0, -1e3}) {
            assertEquals(1e-7 * Math.max(1.0, Math.abs(x)), estimator.computeStep(x), 1e-20);
        }
    }

    /**
     * Checks that default step is relative to x, falling back to the relative step at zero.
     */
    @Test
    void computeStep_whenDefaults_returnsRelativeStepOrFallbackAtZero() {
        final var estimator = new DerivativeEstimator(this);

        assertEquals(DerivativeEstimator.EPS * 10.0, estimator.computeStep(-10.0), 0.0);
        assertEquals(DerivativeEstimator.EPS, estimator.computeStep(0.0), 0.0);
    }

    /**
     * Checks that x equal to zero is still perturbed and derivative is accurate.
     *
     * @throws EvaluationException if function evaluation fails.
     */
    @Test
    void derivative_whenXIsZero_isPerturbed() throws EvaluationException {
        setUpFunction();

        final var defaults = new DerivativeEstimator(this);
        final var custom = new DerivativeEstimator(this, RELATIVE_STEP, MIN_STEP);

        assertEquals(derivative(0.0), defaults.derivative(0.0), ZERO_ERROR);
        assertEquals(derivative(0.0), custom.derivative(0.0), ZERO_ERROR);
    }

    /**
     * Checks that explicitly using default steps gives exactly the same result as the default constructor.
     *
     * @throws EvaluationException if function evaluation fails.
     */
    @Test
    void derivative_whenExplicitDefaultSteps_matchesDefaultConstructor() throws EvaluationException {
        setUpFunction();

        final var defaults = new DerivativeEstimator(this);
        final var explicit = new DerivativeEstimator(this, DerivativeEstimator.EPS, 0.0);

        for (final var x : new double[]{-3.0, 0.0, 0.25, 7.0}) {
            assertEquals(defaults.derivative(x), explicit.derivative(x), 0.0);
        }
    }

    /**
     * Checks that a minimum step improves accuracy at points very close to zero.
     *
     * @throws EvaluationException if function evaluation fails.
     */
    @Test
    void derivative_whenXIsVerySmall_minStepIsMoreAccurateThanDefaults() throws EvaluationException {
        setUpFunction();

        final var defaults = new DerivativeEstimator(this);
        final var custom = new DerivativeEstimator(this, RELATIVE_STEP, MIN_STEP);

        for (final var x : new double[]{1e-5, 1e-8, 1e-10}) {
            final var real = derivative(x);
            final var defaultError = Math.abs(defaults.derivative(x) - real) / Math.abs(real);
            final var customError = Math.abs(custom.derivative(x) - real) / Math.abs(real);

            assertTrue(customError < CUSTOM_RELATIVE_ERROR, "x = " + x + ", error = " + customError);
            assertTrue(customError < defaultError, "x = " + x);
        }
    }

    /**
     * Sets up the function to be evaluated.
     */
    private void setUpFunction() {
        minimum = 0.5;
        width = 1.5;
        offset = 1.0;
    }

    @Override
    public double evaluate(final double point) {
        return (point - minimum) * (point - minimum) / width + offset;
    }

    public double derivative(final double x) {
        return 2.0 * (x - minimum) / width;
    }
}
