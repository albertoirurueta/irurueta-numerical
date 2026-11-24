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
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class JacobianEstimatorTest implements MultiVariateFunctionEvaluatorListener {

    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 3;

    private static final int MIN_VARS = 1;
    private static final int MAX_VARS = 4;

    private static final double MIN_EVAL_POINT = -10.0;
    private static final double MAX_EVAL_POINT = 10.0;

    private static final double MIN_OFFSET = 0.0;
    private static final double MAX_OFFSET = 1.0;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 10.0;

    private static final double ABSOLUTE_ERROR = 1e-1;

    private static final int TIMES = 100;

    /**
     * Default relative step expected for a new estimator.
     */
    private static final double DEFAULT_RELATIVE_STEP = 1e-8;

    /**
     * Default minimum step expected for a new estimator.
     */
    private static final double DEFAULT_MIN_STEP = 0.0;

    /**
     * Custom step used both as relative and as minimum step.
     */
    private static final double CUSTOM_STEP = 1e-7;

    /**
     * Tolerance when comparing a step with its nominal value, accounting for rounding of the perturbed parameter.
     */
    private static final double STEP_TOLERANCE = 1e-13;

    /**
     * Constant offset of the linear function used by the small parameter case.
     */
    private static final double LINEAR_OFFSET = 1.0;

    /**
     * Slope of the linear function used by the small parameter case.
     */
    private static final double LINEAR_SLOPE = 3.0;

    /**
     * Small parameter values used to compare the accuracy of different steps.
     */
    private static final double[] SMALL_PARAMETERS = {1e-5, 1e-7, 1e-10};

    /**
     * Maximum error expected from the estimation using a minimum step of {@link #CUSTOM_STEP}.
     */
    private static final double SMALL_PARAMETER_ERROR = 1e-6;

    /**
     * Number of dimensions used by tests recording perturbations.
     */
    private static final int RECORDING_DIMS = 3;

    private int ndims;
    private int nvars;
    private Matrix minimums;
    private Matrix widths;
    private double[] offsets;

    /**
     * When not null, every point received by {@link #evaluate(double[], double[])} is recorded here.
     */
    private List<double[]> recordedPoints;

    /**
     * When true, {@link #evaluate(double[], double[])} evaluates a linear function of the first parameter.
     */
    private boolean linearMode;

    @Test
    void testConstructor() throws WrongSizeException {

        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        nvars = randomizer.nextInt(MIN_VARS, MAX_VARS);

        minimums = Matrix.createWithUniformRandomValues(nvars, ndims, MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var point = new double[ndims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offsets = new double[nvars];
        randomizer.fill(offsets, MIN_OFFSET, MAX_OFFSET);

        widths = Matrix.createWithUniformRandomValues(nvars, ndims, MIN_WIDTH, MAX_WIDTH);

        final var estimator = new JacobianEstimator(this);
        assertNotNull(estimator);
    }

    @Test
    void testJacobian() throws EvaluationException, WrongSizeException {

        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
            nvars = randomizer.nextInt(MIN_VARS, MAX_VARS);

            minimums = Matrix.createWithUniformRandomValues(nvars, ndims, MIN_EVAL_POINT, MAX_EVAL_POINT);
            final var point = new double[ndims];
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
            offsets = new double[nvars];
            randomizer.fill(offsets, MIN_OFFSET, MAX_OFFSET);

            widths = Matrix.createWithUniformRandomValues(nvars, ndims, MIN_WIDTH, MAX_WIDTH);

            final var estimator = new JacobianEstimator(this);

            final var jacobian1 = estimator.jacobian(point);
            final var jacobian2 = new Matrix(nvars, ndims);
            estimator.jacobian(point, jacobian2);

            // check correctness
            final var jacobian3 = jacobian(point);
            assertNotNull(jacobian3);
            for (var j = 0; j < nvars; j++) {
                for (var i = 0; i < ndims; i++) {
                    assertEquals(jacobian1.getElementAt(j, i), jacobian3.getElementAt(j, i), ABSOLUTE_ERROR);
                    assertEquals(jacobian2.getElementAt(j, i), jacobian3.getElementAt(j, i), ABSOLUTE_ERROR);
                }
            }
        }
    }

    /**
     * Checks that a new estimator uses the default steps.
     */
    @Test
    void constructor_whenDefault_usesDefaultSteps() {
        final var estimator = new JacobianEstimator(this);

        assertEquals(DEFAULT_RELATIVE_STEP, JacobianEstimator.EPS, 0.0);
        assertEquals(JacobianEstimator.EPS, estimator.getRelativeStep(), 0.0);
        assertEquals(DEFAULT_MIN_STEP, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that the constructor stores valid steps.
     */
    @Test
    void constructor_whenValidSteps_storesSteps() {
        final var estimator = new JacobianEstimator(this, CUSTOM_STEP, 2.0 * CUSTOM_STEP);

        assertEquals(CUSTOM_STEP, estimator.getRelativeStep(), 0.0);
        assertEquals(2.0 * CUSTOM_STEP, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that the constructor rejects a non strictly positive relative step.
     *
     * @param relativeStep invalid relative step.
     */
    @ParameterizedTest
    @ValueSource(doubles = {0.0, -1e-8, Double.NaN, Double.POSITIVE_INFINITY})
    void constructor_whenInvalidRelativeStep_throwsIllegalArgumentException(final double relativeStep) {
        assertThrows(IllegalArgumentException.class, () -> new JacobianEstimator(this, relativeStep, 0.0));
    }

    /**
     * Checks that the constructor rejects a negative minimum step.
     *
     * @param minStep invalid minimum step.
     */
    @ParameterizedTest
    @ValueSource(doubles = {-1e-8, Double.NaN, Double.POSITIVE_INFINITY})
    void constructor_whenInvalidMinStep_throwsIllegalArgumentException(final double minStep) {
        assertThrows(IllegalArgumentException.class, () -> new JacobianEstimator(this, CUSTOM_STEP, minStep));
    }

    /**
     * Checks that setters store valid values.
     */
    @Test
    void setSteps_whenValid_updatesSteps() {
        final var estimator = new JacobianEstimator(this);

        estimator.setRelativeStep(CUSTOM_STEP);
        estimator.setMinStep(0.0);
        assertEquals(CUSTOM_STEP, estimator.getRelativeStep(), 0.0);
        assertEquals(0.0, estimator.getMinStep(), 0.0);

        estimator.setMinStep(CUSTOM_STEP);
        assertEquals(CUSTOM_STEP, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that the relative step setter rejects non strictly positive values and keeps previous value.
     *
     * @param relativeStep invalid relative step.
     */
    @ParameterizedTest
    @ValueSource(doubles = {0.0, -1e-8, Double.NaN, Double.POSITIVE_INFINITY})
    void setRelativeStep_whenInvalid_throwsIllegalArgumentException(final double relativeStep) {
        final var estimator = new JacobianEstimator(this);

        assertThrows(IllegalArgumentException.class, () -> estimator.setRelativeStep(relativeStep));
        assertEquals(JacobianEstimator.EPS, estimator.getRelativeStep(), 0.0);
    }

    /**
     * Checks that the minimum step setter rejects negative values and keeps previous value.
     *
     * @param minStep invalid minimum step.
     */
    @ParameterizedTest
    @ValueSource(doubles = {-1e-8, Double.NaN, Double.POSITIVE_INFINITY})
    void setMinStep_whenInvalid_throwsIllegalArgumentException(final double minStep) {
        final var estimator = new JacobianEstimator(this);

        assertThrows(IllegalArgumentException.class, () -> estimator.setMinStep(minStep));
        assertEquals(DEFAULT_MIN_STEP, estimator.getMinStep(), 0.0);
    }

    /**
     * Checks that with default steps results are bit-for-bit equal to the original step rule.
     */
    @Test
    void jacobian_whenDefaultSteps_matchesOriginalRuleExactly() throws EvaluationException, WrongSizeException {
        final var randomizer = new UniformRandomizer();
        for (var t = 0; t < TIMES; t++) {
            setUpFunction(RECORDING_DIMS, MAX_VARS);
            final var point = new double[RECORDING_DIMS];
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
            point[0] = 0.0;

            final var estimator = new JacobianEstimator(this);
            final var result = estimator.jacobian(point);
            final var expected = originalRuleJacobian(point);

            for (var i = 0; i < nvars; i++) {
                for (var j = 0; j < RECORDING_DIMS; j++) {
                    assertEquals(expected.getElementAt(i, j), result.getElementAt(i, j), 0.0);
                }
            }
        }
    }

    /**
     * Checks that a parameter exactly equal to zero still gets a non-zero step with default values.
     */
    @Test
    void jacobian_whenParameterIsZeroAndDefaultSteps_usesNonZeroStep() throws EvaluationException {
        setUpFunction(RECORDING_DIMS, 1);
        final var point = new double[RECORDING_DIMS];
        final var estimator = new JacobianEstimator(this);

        recordedPoints = new ArrayList<>();
        estimator.jacobian(point);

        assertEquals(RECORDING_DIMS + 1, recordedPoints.size());
        for (var j = 0; j < RECORDING_DIMS; j++) {
            final var step = recordedPoints.get(j + 1)[j] - point[j];
            assertNotEquals(0.0, step);
            assertEquals(JacobianEstimator.EPS, step, 0.0);
        }
    }

    /**
     * Checks that using the same relative and minimum step yields a step of {@code step * max(1, |x|)}.
     */
    @Test
    void jacobian_whenRelativeAndMinStepAreEqual_usesStepTimesMaxOfOneAndAbsoluteValue() throws EvaluationException {
        setUpFunction(RECORDING_DIMS, 1);
        final var point = new double[]{0.5, -3.0, 20.0};
        final var estimator = new JacobianEstimator(this, CUSTOM_STEP, CUSTOM_STEP);

        recordedPoints = new ArrayList<>();
        estimator.jacobian(point);

        assertEquals(RECORDING_DIMS + 1, recordedPoints.size());
        for (var j = 0; j < RECORDING_DIMS; j++) {
            final var step = recordedPoints.get(j + 1)[j] - point[j];
            final var expected = CUSTOM_STEP * Math.max(1.0, Math.abs(point[j]));
            assertEquals(expected, step, STEP_TOLERANCE);
        }
    }

    /**
     * Checks that for small parameters a minimum step is closer to the analytical Jacobian than defaults.
     */
    @Test
    void jacobian_whenSmallParametersAndMinStep_isMoreAccurateThanDefaults() throws EvaluationException {
        linearMode = true;
        nvars = 1;
        final var defaults = new JacobianEstimator(this);
        final var withMinStep = new JacobianEstimator(this, CUSTOM_STEP, CUSTOM_STEP);

        for (final var value : SMALL_PARAMETERS) {
            final var point = new double[]{value};

            final var defaultsError = Math.abs(defaults.jacobian(point).getElementAt(0, 0) - LINEAR_SLOPE);
            final var minStepError = Math.abs(withMinStep.jacobian(point).getElementAt(0, 0) - LINEAR_SLOPE);

            assertTrue(minStepError < defaultsError);
            assertTrue(minStepError < SMALL_PARAMETER_ERROR);
        }
    }

    /**
     * Checks that the estimator can be reused when the number of function variables changes.
     */
    @Test
    void jacobian_whenNumberOfVariablesChanges_resizesInternalBuffers() throws EvaluationException {
        final var point = new double[RECORDING_DIMS];
        final var estimator = new JacobianEstimator(this);

        setUpFunction(RECORDING_DIMS, MIN_VARS);
        assertEquals(MIN_VARS, estimator.jacobian(point).getRows());

        setUpFunction(RECORDING_DIMS, MAX_VARS);
        assertEquals(MAX_VARS, estimator.jacobian(point).getRows());
    }

    @Override
    public void evaluate(final double[] point, final double[] result) {
        if (recordedPoints != null) {
            recordedPoints.add(point.clone());
        }
        if (linearMode) {
            result[0] = LINEAR_OFFSET + LINEAR_SLOPE * point[0];
            return;
        }

        final var dims = Math.min(Math.min(point.length, minimums.getColumns()), widths.getColumns());
        final var vars = result.length;

        for (var j = 0; j < vars; j++) {
            var value = 1.0;
            for (var i = 0; i < dims; i++) {
                value *= Math.pow(point[i] - minimums.getElementAt(j, i), 2.0) / widths.getElementAt(j, i);
            }

            value += offsets[j];

            result[j] = value;
        }
    }

    @Override
    public int getNumberOfVariables() {
        return nvars;
    }

    /**
     * Sets up the quadratic test function with deterministic parameters.
     *
     * @param dims number of input dimensions.
     * @param vars number of function variables.
     */
    private void setUpFunction(final int dims, final int vars) {
        ndims = dims;
        nvars = vars;
        linearMode = false;
        recordedPoints = null;
        try {
            minimums = new Matrix(vars, dims);
            widths = Matrix.createWithUniformRandomValues(vars, dims, MIN_WIDTH, MAX_WIDTH);
        } catch (final WrongSizeException e) {
            throw new IllegalStateException(e);
        }
        offsets = new double[vars];
    }

    /**
     * Estimates the Jacobian using the original step rule: {@code EPS * |x|}, or {@code EPS} if zero.
     *
     * @param point input point.
     * @return estimated Jacobian.
     * @throws WrongSizeException if matrix cannot be created.
     */
    private Matrix originalRuleJacobian(final double[] point) throws WrongSizeException {
        final var numdims = point.length;
        final var fold = new double[nvars];
        final var fh = new double[nvars];
        final var xh = point.clone();
        final var result = new Matrix(nvars, numdims);

        evaluate(point, fold);
        for (var j = 0; j < numdims; j++) {
            final var temp = point[j];
            var h = DEFAULT_RELATIVE_STEP * Math.abs(temp);
            if (h == 0.0) {
                h = DEFAULT_RELATIVE_STEP;
            }
            xh[j] = temp + h;
            h = xh[j] - temp;
            evaluate(xh, fh);
            xh[j] = temp;
            for (var i = 0; i < nvars; i++) {
                result.setElementAt(i, j, (fh[i] - fold[i]) / h);
            }
        }
        return result;
    }

    private Matrix jacobian(final double[] params) {
        final var dims = Math.min(Math.min(params.length, minimums.getColumns()), widths.getColumns());
        final var vars = nvars;

        try {
            final var jacobian = new Matrix(vars, dims);
            for (var k = 0; k < vars; k++) {
                for (var j = 0; j < dims; j++) {
                    var value = 1.0;
                    for (var i = 0; i < dims; i++) {
                        if (i != j) {
                            value *= Math.pow(params[i] - minimums.getElementAt(k, i), 2.0)
                                    / widths.getElementAt(k, i);
                        } else {
                            value *= 2.0 * (params[i] - minimums.getElementAt(k, i)) / widths.getElementAt(k, i);
                        }
                    }

                    jacobian.setElementAt(k, j, value);
                }
            }

            return jacobian;
        } catch (final WrongSizeException e) {
            return null;
        }
    }
}
