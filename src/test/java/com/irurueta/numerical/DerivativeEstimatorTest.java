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

class DerivativeEstimatorTest implements SingleDimensionFunctionEvaluatorListener {

    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;

    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;

    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;

    public static final double ABSOLUTE_ERROR = 1e-2;

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

    @Override
    public double evaluate(final double point) {
        return (point - minimum) * (point - minimum) / width + offset;
    }

    public double derivative(final double x) {
        return 2.0 * (x - minimum) / width;
    }
}
