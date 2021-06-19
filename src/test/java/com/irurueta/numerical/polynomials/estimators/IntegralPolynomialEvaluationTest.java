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
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.SerializationHelper;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.io.IOException;
import java.util.Random;

import static org.junit.Assert.*;

public class IntegralPolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final int MIN_ORDER = 1;
    private static final int MAX_ORDER = 5;

    @Test
    public void testConstructor() {
        IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();

        // check default values
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        assertEquals(eval.getX(), 0.0, 0.0);
        assertNull(eval.getConstants());
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);

        // test constructor with values
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double[] constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval = new IntegralPolynomialEvaluation(x, evaluation, constants,
                order);

        // check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertSame(eval.getConstants(), constants);
        assertEquals(eval.getIntegralOrder(), order);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);


        eval = new IntegralPolynomialEvaluation(x, evaluation, constants);

        // check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertSame(eval.getConstants(), constants);
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);


        eval = new IntegralPolynomialEvaluation(x, evaluation, order);

        // check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertNull(eval.getConstants());
        assertEquals(eval.getIntegralOrder(), order);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);


        eval = new IntegralPolynomialEvaluation(x, evaluation);

        // check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertNull(eval.getConstants());
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);
    }

    @Test
    public void testGetSetX() {
        final IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();

        // check default value
        assertEquals(eval.getX(), 0.0, 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setX(x);

        // check correctness
        assertEquals(eval.getX(), x, 0.0);
    }

    @Test
    public void testGetSetConstants() {
        final IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();

        // check default value
        assertNull(eval.getConstants());

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double[] constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setConstant(constants);

        // check correctness
        assertSame(eval.getConstants(), constants);
    }

    @Test
    public void testGetSetIntegralOrder() {
        final IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();

        // check default value
        assertEquals(eval.getIntegralOrder(), 1);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setIntegralOrder(order);

        // check correctness
        assertEquals(eval.getIntegralOrder(), order);

        // Force IllegalArgumentException
        try {
            eval.setIntegralOrder(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double[] constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        final IntegralPolynomialEvaluation eval1 = new IntegralPolynomialEvaluation(x, evaluation, constants,
                order);

        // check default values
        assertEquals(eval1.getEvaluation(), evaluation, 0.0);
        assertEquals(eval1.getX(), x, 0.0);
        assertSame(eval1.getConstants(), constants);
        assertEquals(eval1.getIntegralOrder(), order);
        assertEquals(eval1.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);

        // serialize and deserialize
        final byte[] bytes = SerializationHelper.serialize(eval1);
        final IntegralPolynomialEvaluation eval2 = SerializationHelper.deserialize(bytes);

        // check correctness
        assertEquals(eval2.getEvaluation(), evaluation, 0.0);
        assertEquals(eval2.getX(), x, 0.0);
        assertArrayEquals(eval2.getConstants(), constants, 0.0);
        assertEquals(eval2.getIntegralOrder(), order);
        assertEquals(eval2.getType(),
                PolynomialEvaluationType.INTEGRAL_EVALUATION);
    }
}
