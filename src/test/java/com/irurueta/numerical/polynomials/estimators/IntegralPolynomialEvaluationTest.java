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
import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

class IntegralPolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final int MIN_ORDER = 1;
    private static final int MAX_ORDER = 5;

    @Test
    void testConstructor() {
        var eval = new IntegralPolynomialEvaluation();

        // check default values
        assertEquals(0.0, eval.getEvaluation(), 0.0);
        assertEquals(0.0, eval.getX(), 0.0);
        assertNull(eval.getConstants());
        assertEquals(1, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval.getType());

        // test constructor with values
        final var randomizer = new UniformRandomizer();

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval = new IntegralPolynomialEvaluation(x, evaluation, constants, order);

        // check default values
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
        assertEquals(x, eval.getX(), 0.0);
        assertSame(constants, eval.getConstants());
        assertEquals(order, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval.getType());

        eval = new IntegralPolynomialEvaluation(x, evaluation, constants);

        // check default values
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
        assertEquals(x, eval.getX(), 0.0);
        assertSame(constants, eval.getConstants());
        assertEquals(1, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval.getType());

        eval = new IntegralPolynomialEvaluation(x, evaluation, order);

        // check default values
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
        assertEquals(x, eval.getX(), 0.0);
        assertNull(eval.getConstants());
        assertEquals(order, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval.getType());

        eval = new IntegralPolynomialEvaluation(x, evaluation);

        // check default values
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
        assertEquals(x, eval.getX(), 0.0);
        assertNull(eval.getConstants());
        assertEquals(1, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval.getType());
    }

    @Test
    void testGetSetX() {
        final var eval = new IntegralPolynomialEvaluation();

        // check default value
        assertEquals(0.0, eval.getX(), 0.0);

        // set new value
        final var randomizer = new UniformRandomizer();

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setX(x);

        // check correctness
        assertEquals(x, eval.getX(), 0.0);
    }

    @Test
    void testGetSetConstants() {
        final var eval = new IntegralPolynomialEvaluation();

        // check default value
        assertNull(eval.getConstants());

        // set new value
        final var randomizer = new UniformRandomizer();

        final var constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setConstant(constants);

        // check correctness
        assertSame(constants, eval.getConstants());
    }

    @Test
    void testGetSetIntegralOrder() {
        final var eval = new IntegralPolynomialEvaluation();

        // check default value
        assertEquals(1, eval.getIntegralOrder());

        // set new value
        final var randomizer = new UniformRandomizer();

        final var order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setIntegralOrder(order);

        // check correctness
        assertEquals(order, eval.getIntegralOrder());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> eval.setIntegralOrder(0));
    }

    @Test
    void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final var randomizer = new UniformRandomizer();

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        final var eval1 = new IntegralPolynomialEvaluation(x, evaluation, constants, order);

        // check default values
        assertEquals(evaluation, eval1.getEvaluation(), 0.0);
        assertEquals(x, eval1.getX(), 0.0);
        assertSame(constants, eval1.getConstants());
        assertEquals(order, eval1.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval1.getType());

        // serialize and deserialize
        final var bytes = SerializationHelper.serialize(eval1);
        final var eval2 = SerializationHelper.<IntegralPolynomialEvaluation>deserialize(bytes);

        // check correctness
        assertEquals(evaluation, eval2.getEvaluation(), 0.0);
        assertEquals(x, eval2.getX(), 0.0);
        assertArrayEquals(constants, eval2.getConstants(), 0.0);
        assertEquals(order, eval2.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_EVALUATION, eval2.getType());
    }
}
