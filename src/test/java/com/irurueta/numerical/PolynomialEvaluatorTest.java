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

import static org.junit.Assert.*;

public class PolynomialEvaluatorTest {

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;
    private static final int MIN_LENGTH = 1;
    private static final int MAX_LENGTH = 5;

    private static final double ABSOLUTE_ERROR = 1e-9;

    private static final int TIMES = 10;

    @Test
    public void testEvaluateRealConstant() {

        final double[] polyParams = new double[1];

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0], 0.0);

        // Force IllegalArgumentException
        double evaluate = 0.0;
        try {
            evaluate = PolynomialEvaluator.evaluate(null, x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            evaluate = PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertEquals(0.0, evaluate, 0.0);
    }

    @Test
    public void testEvaluateRealFirstDegree() {

        final double[] polyParams = new double[2];

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0] * x + polyParams[1], ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        double evaluate = 0.0;
        try {
            evaluate = PolynomialEvaluator.evaluate(null, x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            evaluate = PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertEquals(0.0, evaluate, 0.0);
    }

    @Test
    public void testEvaluateRealSecondDegree() {

        final double[] polyParams = new double[3];

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0] * x * x + polyParams[1] * x + polyParams[2],
                ABSOLUTE_ERROR);

        // Force IllegalArgumentException
        double evaluate = 0.0;
        try {
            evaluate = PolynomialEvaluator.evaluate(null, x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            evaluate = PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertEquals(0.0, evaluate, 0.0);
    }

    @Test
    public void testEvaluateReal() {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            final int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
            final int degree = length - 1;

            final double[] polyParams = new double[length];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            double value = 0.0;
            for (int i = 0; i < length; i++) {
                value += polyParams[i] * Math.pow(x, degree - i);
            }

            final double value2 = PolynomialEvaluator.evaluate(polyParams, x);
            if (Math.abs(value2 - value) > ABSOLUTE_ERROR) {
                continue;
            }

            assertEquals(value2, value, ABSOLUTE_ERROR);

            // Force IllegalArgumentException
            double evaluate = 0.0;
            try {
                evaluate = PolynomialEvaluator.evaluate(null, x);
                fail("IllegalArgumentException expected but not thrown");
            } catch (final IllegalArgumentException ignore) {
            }
            try {
                evaluate = PolynomialEvaluator.evaluate(new double[0], x);
                fail("IllegalArgumentException expected but not thrown");
            } catch (final IllegalArgumentException ignore) {
            }
            assertEquals(0.0, evaluate, 0.0);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testEvaluateComplexConstant() {

        final Complex[] polyParams = new Complex[1];

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for (int i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }

        final Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE));

        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0]);

        // Force IllegalArgumentException
        Complex evaluate = null;
        try {
            evaluate = PolynomialEvaluator.evaluate(null, x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            evaluate = PolynomialEvaluator.evaluate(new Complex[0], x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(evaluate);
    }

    @Test
    public void testEvaluateComplexFirstDegree() {

        final Complex[] polyParams = new Complex[2];

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for (int i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }

        final Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE));

        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0].multiplyAndReturnNew(x).
                        addAndReturnNew(polyParams[1]));

        // Force IllegalArgumentException
        Complex evaluate = null;
        try {
            evaluate = PolynomialEvaluator.evaluate(null, x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            evaluate = PolynomialEvaluator.evaluate(new Complex[0], x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(evaluate);
    }

    @Test
    public void testEvaluateComplexSecondDegree() {

        final Complex[] polyParams = new Complex[3];

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for (int i = 0; i < polyParams.length; i++) {
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }

        final Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE));

        assertTrue(PolynomialEvaluator.evaluate(polyParams, x).equals(
                polyParams[0].multiplyAndReturnNew(x.multiplyAndReturnNew(x)).
                        addAndReturnNew(polyParams[1].multiplyAndReturnNew(x)).
                        addAndReturnNew(polyParams[2]), ABSOLUTE_ERROR));

        // Force IllegalArgumentException
        Complex evaluate = null;
        try {
            evaluate = PolynomialEvaluator.evaluate(null, x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            evaluate = PolynomialEvaluator.evaluate(new Complex[0], x);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(evaluate);
    }

    @Test
    public void testEvaluateComplex() {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            final int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
            final int degree = length - 1;

            final Complex[] polyParams = new Complex[length];
            for (int i = 0; i < polyParams.length; i++) {
                polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));
            }

            final Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));

            final Complex value = new Complex();
            for (int i = 0; i < length; i++) {
                value.add(polyParams[i].multiplyAndReturnNew(
                        x.powAndReturnNew(degree - i)));
            }

            if (!PolynomialEvaluator.evaluate(polyParams, x).equals(value, ABSOLUTE_ERROR)) {
                continue;
            }
            assertTrue(PolynomialEvaluator.evaluate(polyParams, x).
                    equals(value, ABSOLUTE_ERROR));

            // Force IllegalArgumentException
            Complex evaluate = null;
            try {
                evaluate = PolynomialEvaluator.evaluate(null, x);
                fail("IllegalArgumentException expected but not thrown");
            } catch (final IllegalArgumentException ignore) {
            }
            try {
                evaluate = PolynomialEvaluator.evaluate(new Complex[0], x);
                fail("IllegalArgumentException expected but not thrown");
            } catch (final IllegalArgumentException ignore) {
            }
            assertNull(evaluate);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }
}
