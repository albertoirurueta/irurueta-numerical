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

public class IntegralIntervalPolynomialEvaluationTest {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final int MIN_ORDER = 1;
    private static final int MAX_ORDER = 5;

    @Test
    public void testConstructor() {
        IntegralIntervalPolynomialEvaluation eval =
                new IntegralIntervalPolynomialEvaluation();

        // check default values
        assertEquals(0.0, eval.getEvaluation(), 0.0);
        assertEquals(0.0, eval.getStartX(), 0.0);
        assertEquals(0.0, eval.getEndX(), 0.0);
        assertEquals(1, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_INTERVAL,
                eval.getType());

        // test constructor with values
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double startX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        eval = new IntegralIntervalPolynomialEvaluation(startX, endX,
                evaluation, order);

        // check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getStartX(), startX, 0.0);
        assertEquals(eval.getEndX(), endX, 0.0);
        assertEquals(eval.getIntegralOrder(), order);
        assertEquals(PolynomialEvaluationType.INTEGRAL_INTERVAL,
                eval.getType());


        eval = new IntegralIntervalPolynomialEvaluation(startX, endX,
                evaluation);

        // check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getStartX(), startX, 0.0);
        assertEquals(eval.getEndX(), endX, 0.0);
        assertEquals(1, eval.getIntegralOrder());
        assertEquals(PolynomialEvaluationType.INTEGRAL_INTERVAL,
                eval.getType());
    }

    @Test
    public void testGetSetStartX() {
        final IntegralIntervalPolynomialEvaluation eval =
                new IntegralIntervalPolynomialEvaluation();

        // check default value
        assertEquals(0.0, eval.getStartX(), 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double startX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        eval.setStartX(startX);

        // check correctness
        assertEquals(eval.getStartX(), startX, 0.0);
    }

    @Test
    public void testGetSetEndX() {
        final IntegralIntervalPolynomialEvaluation eval =
                new IntegralIntervalPolynomialEvaluation();

        // check default value
        assertEquals(0.0, eval.getEndX(), 0.0);

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setEndX(endX);

        // check correctness
        assertEquals(eval.getEndX(), endX, 0.0);
    }

    @Test
    public void testGetSetIntegralOrder() {
        final IntegralIntervalPolynomialEvaluation eval =
                new IntegralIntervalPolynomialEvaluation();

        // check default value
        assertEquals(1, eval.getIntegralOrder());

        // set new value
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setIntegralOrder(order);

        // check correctness
        assertEquals(eval.getIntegralOrder(), order);
    }

    @Test
    public void testSerializeDeserialize() throws IOException, ClassNotFoundException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double startX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE,
                MAX_RANDOM_VALUE);
        final int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);

        final IntegralIntervalPolynomialEvaluation eval1 = new IntegralIntervalPolynomialEvaluation(startX, endX,
                evaluation, order);

        // check default values
        assertEquals(eval1.getEvaluation(), evaluation, 0.0);
        assertEquals(eval1.getStartX(), startX, 0.0);
        assertEquals(eval1.getEndX(), endX, 0.0);
        assertEquals(eval1.getIntegralOrder(), order);
        assertEquals(PolynomialEvaluationType.INTEGRAL_INTERVAL,
                eval1.getType());

        // serialize and deserialize
        final byte[] bytes = SerializationHelper.serialize(eval1);
        final IntegralIntervalPolynomialEvaluation eval2 = SerializationHelper.deserialize(bytes);

        // check correctness
        assertEquals(eval2.getEvaluation(), evaluation, 0.0);
        assertEquals(eval2.getStartX(), startX, 0.0);
        assertEquals(eval2.getEndX(), endX, 0.0);
        assertEquals(eval2.getIntegralOrder(), order);
        assertEquals(PolynomialEvaluationType.INTEGRAL_INTERVAL,
                eval2.getType());
    }
}
