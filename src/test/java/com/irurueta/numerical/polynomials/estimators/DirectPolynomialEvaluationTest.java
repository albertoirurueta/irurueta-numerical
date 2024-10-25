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

import static org.junit.jupiter.api.Assertions.assertEquals;

class DirectPolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    @Test
    void testConstructor() {
        // test empty constructor
        var eval = new DirectPolynomialEvaluation();

        // check initial values
        assertEquals(0.0, eval.getX(), 0.0);
        assertEquals(0.0, eval.getEvaluation(), 0.0);
        assertEquals(PolynomialEvaluationType.DIRECT_EVALUATION, eval.getType());

        // test constructor with values
        final var randomizer = new UniformRandomizer();

        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        eval = new DirectPolynomialEvaluation(x, evaluation);

        // check correctness
        assertEquals(x, eval.getX(), 0.0);
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
    }

    @Test
    void testGetSetX() {
        final var eval = new DirectPolynomialEvaluation();

        // check initial value
        assertEquals(0.0, eval.getX(), 0.0);

        // set new value
        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setX(x);

        // check correctness
        assertEquals(x, eval.getX(), 0.0);
    }

    @Test
    void testGetSetEvaluation() {
        final var eval = new DirectPolynomialEvaluation();

        // check initial value
        assertEquals(0.0, eval.getEvaluation(), 0.0);

        // set new value
        final var randomizer = new UniformRandomizer();
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setEvaluation(evaluation);

        // check correctness
        assertEquals(evaluation, eval.getEvaluation(), 0.0);
    }

    @Test
    void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final var eval1 = new DirectPolynomialEvaluation();

        // set new values
        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval1.setX(x);
        final var evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval1.setEvaluation(evaluation);

        // check correctness
        assertEquals(x, eval1.getX(), 0.0);
        assertEquals(evaluation, eval1.getEvaluation(), 0.0);

        // serialize and deserialize
        final var bytes = SerializationHelper.serialize(eval1);
        final var eval2 = SerializationHelper.<DirectPolynomialEvaluation>deserialize(bytes);

        // check correctness
        assertEquals(eval2.getX(), x, 0.0);
        assertEquals(eval2.getEvaluation(), evaluation, 0.0);
    }
}
