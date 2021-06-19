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
import static org.junit.Assert.fail;

public class DerivativePolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;
    private static final int MIN_ORDER = 1;
    private static final int MAX_ORDER = 4;

    @Test
    public void testConstructor() {
        // test empty constructor
        DerivativePolynomialEvaluation eval =
                new DerivativePolynomialEvaluation();

        // check initial values
        assertEquals(eval.getDerivativeOrder(), 1);
        assertEquals(eval.getX(), 0.0, 0.0);
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.DERIVATIVE_EVALUATION);

        // test constructor with values
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval = new DerivativePolynomialEvaluation(x, evaluation, order);

        // check correctness
        assertEquals(eval.getDerivativeOrder(), order);
        assertEquals(eval.getX(), x, 0.0);
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.DERIVATIVE_EVALUATION);

        // test constructor with values (without order)
        eval = new DerivativePolynomialEvaluation(x, evaluation);

        // check correctness
        assertEquals(eval.getDerivativeOrder(), 1);
        assertEquals(eval.getX(), x, 0.0);
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getType(),
                PolynomialEvaluationType.DERIVATIVE_EVALUATION);
    }

    @Test
    public void testGetSetEvaluation() {
        final DerivativePolynomialEvaluation eval =
                new DerivativePolynomialEvaluation();

        // check default value
        assertEquals(eval.getEvaluation(), 0.0, 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        eval.setEvaluation(evaluation);

        // check correctness
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
    }

    @Test
    public void testGetSetX() {
        final DerivativePolynomialEvaluation eval =
                new DerivativePolynomialEvaluation();

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
    public void testGetSetDerivativeOrder() {
        final DerivativePolynomialEvaluation eval =
                new DerivativePolynomialEvaluation();

        // check default value
        assertEquals(eval.getDerivativeOrder(), 1);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setDerivativeOrder(order);

        // check correctness
        assertEquals(eval.getDerivativeOrder(), order);

        // Force IllegalArgumentException
        try {
            eval.setDerivativeOrder(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final DerivativePolynomialEvaluation eval1 =
                new DerivativePolynomialEvaluation();

        // set new values
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval1.setEvaluation(evaluation);
        eval1.setX(x);
        eval1.setDerivativeOrder(order);

        // check correctness
        assertEquals(eval1.getEvaluation(), evaluation, 0.0);
        assertEquals(eval1.getX(), x, 0.0);
        assertEquals(eval1.getDerivativeOrder(), order);

        // serialize and deserialize
        final byte[] bytes = SerializationHelper.serialize(eval1);
        final DerivativePolynomialEvaluation eval2 = SerializationHelper.deserialize(bytes);

        // check correctness
        assertEquals(eval2.getEvaluation(), evaluation, 0.0);
        assertEquals(eval2.getX(), x, 0.0);
        assertEquals(eval2.getDerivativeOrder(), order);
    }
}
