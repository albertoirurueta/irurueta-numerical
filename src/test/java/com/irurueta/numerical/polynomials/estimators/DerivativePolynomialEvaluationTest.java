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

class DerivativePolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;
    private static final int MIN_ORDER = 1;
    private static final int MAX_ORDER = 4;

    @Test
    void testConstructor() {
        // test empty constructor
        var eval = new DerivativePolynomialEvaluation();

        // check initial values
        assertEquals(1, eval.getDerivativeOrder());
        assertEquals(0.0, eval.getX(), 0.0);
        assertEquals(0.0, eval.getEvaluation(), 0.0);
        assertEquals(PolynomialEvaluationType.DERIVATIVE_EVALUATION, eval.getType());

        // test constructor with values
        final var randomizer = new UniformRandomizer();

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval = new DerivativePolynomialEvaluation(x, evaluation, order);

        // check correctness
        assertEquals(order, eval.getDerivativeOrder());
        assertEquals(x, eval.getX(), 0.0);
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
        assertEquals(PolynomialEvaluationType.DERIVATIVE_EVALUATION, eval.getType());

        // test constructor with values (without order)
        eval = new DerivativePolynomialEvaluation(x, evaluation);

        // check correctness
        assertEquals(1, eval.getDerivativeOrder());
        assertEquals(x, eval.getX(), 0.0);
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
        assertEquals(PolynomialEvaluationType.DERIVATIVE_EVALUATION, eval.getType());
    }

    @Test
    void testGetSetEvaluation() {
        final var eval = new DerivativePolynomialEvaluation();

        // check default value
        assertEquals(0.0, eval.getEvaluation(), 0.0);

        // set new value
        final var randomizer = new UniformRandomizer();
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setEvaluation(evaluation);

        // check correctness
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
    }

    @Test
    void testGetSetX() {
        final var eval = new DerivativePolynomialEvaluation();

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
    void testGetSetDerivativeOrder() {
        final var eval = new DerivativePolynomialEvaluation();

        // check default value
        assertEquals(1, eval.getDerivativeOrder());

        // set new value
        final var randomizer = new UniformRandomizer();
        final var order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setDerivativeOrder(order);

        // check correctness
        assertEquals(order, eval.getDerivativeOrder());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> eval.setDerivativeOrder(0));
    }

    @Test
    void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final var eval1 = new DerivativePolynomialEvaluation();

        // set new values
        final var randomizer = new UniformRandomizer();
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval1.setEvaluation(evaluation);
        eval1.setX(x);
        eval1.setDerivativeOrder(order);

        // check correctness
        assertEquals(evaluation, eval1.getEvaluation(), 0.0);
        assertEquals(x, eval1.getX(), 0.0);
        assertEquals(order, eval1.getDerivativeOrder());

        // serialize and deserialize
        final var bytes = SerializationHelper.serialize(eval1);
        final var eval2 = SerializationHelper.<DerivativePolynomialEvaluation>deserialize(bytes);

        // check correctness
        assertEquals(evaluation, eval2.getEvaluation(), 0.0);
        assertEquals(x, eval2.getX(), 0.0);
        assertEquals(order, eval2.getDerivativeOrder());
    }
}
