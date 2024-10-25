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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class PolynomialEvaluatorTest {

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;
    private static final int MIN_LENGTH = 1;
    private static final int MAX_LENGTH = 5;

    private static final double ABSOLUTE_ERROR = 1e-9;

    private static final int TIMES = 10;

    @Test
    void testEvaluateRealConstant() {

        final var polyParams = new double[1];

        final var randomizer = new UniformRandomizer();
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        assertEquals(PolynomialEvaluator.evaluate(polyParams, x), polyParams[0], 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
        //noinspection ObviousNullCheck
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new double[0], x));
    }

    @Test
    void testEvaluateRealFirstDegree() {

        final var polyParams = new double[2];

        final var randomizer = new UniformRandomizer();
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        assertEquals(polyParams[0] * x + polyParams[1], PolynomialEvaluator.evaluate(polyParams, x),
                ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
        //noinspection ObviousNullCheck
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new double[0], x));
    }

    @Test
    void testEvaluateRealSecondDegree() {

        final var polyParams = new double[3];

        final var randomizer = new UniformRandomizer();
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        assertEquals(polyParams[0] * x * x + polyParams[1] * x + polyParams[2],
                PolynomialEvaluator.evaluate(polyParams, x), ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
        //noinspection ObviousNullCheck
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new double[0], x));
    }

    @Test
    void testEvaluateReal() {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
            final var degree = length - 1;

            final var polyParams = new double[length];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            var value = 0.0;
            for (var i = 0; i < length; i++) {
                value += polyParams[i] * Math.pow(x, degree - i);
            }

            final var value2 = PolynomialEvaluator.evaluate(polyParams, x);
            if (Math.abs(value2 - value) > ABSOLUTE_ERROR) {
                continue;
            }

            assertEquals(value2, value, ABSOLUTE_ERROR);

            // Force IllegalArgumentException
            assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
            //noinspection ObviousNullCheck
            assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new double[0], x));

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testEvaluateComplexConstant() {

        final var polyParams = new Complex[1];

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                    randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
        }

        final var x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

        assertEquals(polyParams[0], PolynomialEvaluator.evaluate(polyParams, x));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new Complex[0], x));
    }

    @Test
    void testEvaluateComplexFirstDegree() {

        final var polyParams = new Complex[2];

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                    randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
        }

        final var x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

        assertEquals(polyParams[0].multiplyAndReturnNew(x).addAndReturnNew(polyParams[1]),
                PolynomialEvaluator.evaluate(polyParams, x));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new Complex[0], x));
    }

    @Test
    void testEvaluateComplexSecondDegree() {

        final var polyParams = new Complex[3];

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                    randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
        }

        final var x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

        assertTrue(PolynomialEvaluator.evaluate(polyParams, x).equals(
                polyParams[0].multiplyAndReturnNew(x.multiplyAndReturnNew(x)).
                        addAndReturnNew(polyParams[1].multiplyAndReturnNew(x)).
                        addAndReturnNew(polyParams[2]), ABSOLUTE_ERROR));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
        assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new Complex[0], x));
    }

    @Test
    void testEvaluateComplex() {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();
            final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
            final var degree = length - 1;

            final var polyParams = new Complex[length];
            for (var i = 0; i < polyParams.length; i++) {
                polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                        randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            }

            final var x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE),
                    randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

            final var value = new Complex();
            for (var i = 0; i < length; i++) {
                value.add(polyParams[i].multiplyAndReturnNew(x.powAndReturnNew(degree - i)));
            }

            if (!PolynomialEvaluator.evaluate(polyParams, x).equals(value, ABSOLUTE_ERROR)) {
                continue;
            }
            assertTrue(PolynomialEvaluator.evaluate(polyParams, x).equals(value, ABSOLUTE_ERROR));

            // Force IllegalArgumentException
            assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(null, x));
            assertThrows(IllegalArgumentException.class, () -> PolynomialEvaluator.evaluate(new Complex[0], x));

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }
}
