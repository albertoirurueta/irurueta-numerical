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

class SavitzkyGolayGradientEstimatorTest implements MultiDimensionFunctionEvaluatorListener {

    private static final int MIN_DIMS = 1;
    private static final int MAX_DIMS = 3;

    private static final double MIN_EVAL_POINT = -10.0;
    private static final double MAX_EVAL_POINT = 10.0;

    private static final double MIN_OFFSET = -10.0;
    private static final double MAX_OFFSET = 10.0;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;

    private static final double ABSOLUTE_ERROR = 1e-2;

    private static final int TIMES = 100;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    @Test
    void testConstructor() {

        final var randomizer = new UniformRandomizer();
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        minimum = new double[ndims];
        final var point = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

        final var estimator = new SavitzkyGolayGradientEstimator(this);
        assertNotNull(estimator);
    }

    @Test
    void testGradient() throws EvaluationException {

        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            final var point = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

            final var estimator = new SavitzkyGolayGradientEstimator(this);

            final var gradient1 = estimator.gradient(point);

            final var gradient2 = new double[ndims];
            estimator.gradient(point, gradient2);

            final var gradient3 = gradient(point);

            // check correctness
            assertEquals(gradient1.length, ndims);
            assertEquals(gradient2.length, ndims);
            assertEquals(gradient3.length, ndims);
            for (var i = 0; i < ndims; i++) {
                assertEquals(gradient1[i], gradient3[i], 5 * ABSOLUTE_ERROR);
                assertEquals(gradient2[i], gradient3[i], 5 * ABSOLUTE_ERROR);
                assertEquals(gradient1[i], gradient2[i], 0.0);
            }
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
}
