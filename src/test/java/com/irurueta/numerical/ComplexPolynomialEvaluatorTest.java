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

import com.irurueta.algebra.Complex;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

public class ComplexPolynomialEvaluatorTest {

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;
    private static final int MIN_LENGTH = 1;
    private static final int MAX_LENGTH = 5;

    @Test
    public void testConstructor() {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final Complex[] polyParams = new Complex[length];

        final ComplexPolynomialEvaluator evaluator =
                new ComplexPolynomialEvaluator(polyParams);

        // check correctness
        assertSame(evaluator.getPolyParams(), polyParams);
    }

    @Test
    public void testEvaluate() {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final Complex[] polyParams = new Complex[length];
        for (int i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }

        final Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE));

        final ComplexPolynomialEvaluator evaluator =
                new ComplexPolynomialEvaluator(polyParams);

        assertEquals(evaluator.evaluate(x),
                PolynomialEvaluator.evaluate(polyParams, x));
    }
}
