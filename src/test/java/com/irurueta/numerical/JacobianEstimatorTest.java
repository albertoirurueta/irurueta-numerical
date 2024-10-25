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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

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

    private int ndims;
    private int nvars;
    private Matrix minimums;
    private Matrix widths;
    private double[] offsets;

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

    @Override
    public void evaluate(final double[] point, final double[] result) {
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
