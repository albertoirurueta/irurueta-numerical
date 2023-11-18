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

import static org.junit.Assert.assertEquals;

public class DirectPolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    @Test
    public void testConstructor() {
        // test empty constructor
        DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation();

        // check initial values
        assertEquals(0.0, eval.getX(), 0.0);
        assertEquals(0.0, eval.getEvaluation(), 0.0);
        assertEquals(PolynomialEvaluationType.DIRECT_EVALUATION,
                eval.getType());

        // test constructor with values
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);

        eval = new DirectPolynomialEvaluation(x, evaluation);

        // check correctness
        assertEquals(eval.getX(), x, 0.0);
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
    }

    @Test
    public void testGetSetX() {
        final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation();

        // check initial value
        assertEquals(0.0, eval.getX(), 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setX(x);

        // check correctness
        assertEquals(eval.getX(), x, 0.0);
    }

    @Test
    public void testGetSetEvaluation() {
        final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation();

        // check initial value
        assertEquals(0.0, eval.getEvaluation(), 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        eval.setEvaluation(evaluation);

        // check correctness
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
    }

    @Test
    public void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final DirectPolynomialEvaluation eval1 = new DirectPolynomialEvaluation();

        // set new values
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval1.setX(x);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        eval1.setEvaluation(evaluation);

        // check correctness
        assertEquals(eval1.getX(), x, 0.0);
        assertEquals(eval1.getEvaluation(), evaluation, 0.0);

        // serialize and deserialize
        final byte[] bytes = SerializationHelper.serialize(eval1);
        final DirectPolynomialEvaluation eval2 = SerializationHelper.deserialize(bytes);

        // check correctness
        assertEquals(eval2.getX(), x, 0.0);
        assertEquals(eval2.getEvaluation(), evaluation, 0.0);
    }
}
