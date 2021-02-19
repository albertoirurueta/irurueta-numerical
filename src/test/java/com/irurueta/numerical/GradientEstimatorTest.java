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
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class GradientEstimatorTest
        implements MultiDimensionFunctionEvaluatorListener {

    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 3;

    private static final double MIN_EVAL_POINT = -10.0;
    private static final double MAX_EVAL_POINT = 10.0;

    private static final double MIN_OFFSET = 0.0;
    private static final double MAX_OFFSET = 1.0;

    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 10.0;

    private static final double ABSOLUTE_ERROR = 1e-1;

    private static final int TIMES = 100;

    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;

    @Test
    public void testConstructor() {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        minimum = new double[ndims];
        final double[] point = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

        final GradientEstimator estimator = new GradientEstimator(this);
        assertNotNull(estimator);
    }

    @Test
    public void testGradient() throws EvaluationException {

        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

            minimum = new double[ndims];
            final double[] point = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);

            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);

            final GradientEstimator estimator = new GradientEstimator(this);

            final double[] gradient1 = estimator.gradient(point);
            final double[] gradient2 = new double[ndims];
            estimator.gradient(point, gradient2);

            // check correctness
            final double[] gradient3 = gradient(point);
            for (int i = 0; i < ndims; i++) {
                assertEquals(gradient1[i], gradient3[i], ABSOLUTE_ERROR);
                assertEquals(gradient2[i], gradient3[i], ABSOLUTE_ERROR);
            }
        }
    }

    @Override
    public double evaluate(final double[] point) throws EvaluationException {
        final int dims = Math.min(Math.min(point.length, minimum.length),
                width.length);

        double value = 1.0;

        for (int i = 0; i < dims; i++) {
            value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
        }

        value += offset;

        return value;
    }

    public double[] gradient(final double[] params) {
        final int dims = Math.min(Math.min(params.length, minimum.length),
                width.length);

        final double[] gradient = new double[dims];

        double value;
        for (int j = 0; j < dims; j++) {
            value = 1.0;
            for (int i = 0; i < dims; i++) {
                if (i != j) {
                    value *= Math.pow(params[i] - minimum[i], 2.0) /
                            width[i];
                } else {
                    value *= 2.0 * (params[i] - minimum[i]) / width[i];
                }
            }

            gradient[j] = value;
        }

        return gradient;
    }
}
