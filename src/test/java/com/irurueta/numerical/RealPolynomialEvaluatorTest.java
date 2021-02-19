/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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
import static org.junit.Assert.assertSame;

public class RealPolynomialEvaluatorTest {

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;
    private static final int MIN_LENGTH = 1;
    private static final int MAX_LENGTH = 5;

    @Test
    public void testConstructor() {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final double[] polyParams = new double[length];

        final RealPolynomialEvaluator evaluator = new RealPolynomialEvaluator(
                polyParams);

        // check correctness
        assertSame(evaluator.getPolyParams(), polyParams);
    }

    @Test
    public void testEvaluate() {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final double[] polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final RealPolynomialEvaluator evaluator = new RealPolynomialEvaluator(
                polyParams);

        assertEquals(evaluator.evaluate(x),
                PolynomialEvaluator.evaluate(polyParams, x), 0.0);
    }
}
