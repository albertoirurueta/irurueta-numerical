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

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.*;

class LMSEPolynomialEstimatorTest implements PolynomialEstimatorListener {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final int MIN_DEGREE = 1;
    private static final int MAX_DEGREE = 5;

    private static final double ABSOLUTE_ERROR = 1e-8;

    private int estimateStart;
    private int estimateEnd;

    @Test
    void testConstructor() {
        // empty constructor
        var estimator = new LMSEPolynomialEstimator();

        // check correctness
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // constructor with degree
        estimator = new LMSEPolynomialEstimator(2);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LMSEPolynomialEstimator(0));

        // constructor with evaluations
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        estimator = new LMSEPolynomialEstimator(evaluations);

        // check correctness
        assertEquals(1, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // constructor with listener
        estimator = new LMSEPolynomialEstimator(this);

        // check correctness
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // constructor with degree and evaluations
        estimator = new LMSEPolynomialEstimator(2, evaluations);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LMSEPolynomialEstimator(0, evaluations));

        // constructor with degree and listener
        estimator = new LMSEPolynomialEstimator(2, this);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LMSEPolynomialEstimator(0, this));

        // constructor with evaluations and listener
        estimator = new LMSEPolynomialEstimator(evaluations, this);

        // check correctness
        assertEquals(1, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // constructor with degree, evaluations and listener
        estimator = new LMSEPolynomialEstimator(2, evaluations, this);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LMSEPolynomialEstimator(0, evaluations,
                this));
    }

    @Test
    void testIsSetLMSESolutionAllowed() throws LockedException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default value
        assertFalse(estimator.isLMSESolutionAllowed());

        // set new value
        estimator.setLMSESolutionAllowed(true);

        // check correctness
        assertTrue(estimator.isLMSESolutionAllowed());
    }

    @Test
    void testGetSetDegree() throws LockedException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default value
        assertEquals(1, estimator.getDegree());

        // set new value
        estimator.setDegree(2);

        // check correctness
        assertEquals(2, estimator.getDegree());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegree(0));
    }

    @Test
    void testGetSetEvaluations() throws LockedException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default value
        assertNull(estimator.getEvaluations());

        // set new value
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        estimator.setEvaluations(evaluations);

        // check correctness
        assertSame(evaluations, estimator.getEvaluations());
    }

    @Test
    void testSetDegreeAndEvaluations() throws LockedException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());

        // set new values
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        estimator.setDegreeAndEvaluations(2, evaluations);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeAndEvaluations(0, evaluations));
    }

    @Test
    void testIsReady() throws LockedException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default value
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
    }

    @Test
    void testGetMinNumberOfEvaluations() {
        final var randomizer = new UniformRandomizer();

        final var degree = randomizer.nextInt(MIN_DEGREE, MAX_DEGREE);
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(degree), degree + 1);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEstimator.getMinNumberOfEvaluations(0));

        final var estimator = new LMSEPolynomialEstimator(degree);
        assertEquals(degree + 1, estimator.getMinNumberOfEvaluations());
    }

    @Test
    void testGetSetListener() throws LockedException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default value
        assertNull(estimator.getListener());

        // set new value
        estimator.setListener(this);

        // check correctness
        assertSame(this, estimator.getListener());
    }

    @Test
    void testEstimateWithDirectEvaluationsNoLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);


        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectEvaluationsLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectAndDerivativeEvaluationsNoLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations() - 1; i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }
        final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        final var value = polynomial.evaluateDerivative(x);
        final var eval = new DerivativePolynomialEvaluation(x, value, 1);
        evaluations.add(eval);

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectAndDerivativeEvaluationLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluateDerivative(x);
            final var eval = new DerivativePolynomialEvaluation(x, value, 1);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithIntegralEvaluationsNoLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var integral = polynomial.integrationAndReturnNew(constant);
            final var value = integral.evaluate(x);

            final var eval = new IntegralPolynomialEvaluation(x, value, new double[]{constant}, 1);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithIntegralEvaluationsLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var integral = polynomial.integrationAndReturnNew(constant);
            final var value = integral.evaluate(x);

            final var eval = new IntegralPolynomialEvaluation(x, value, new double[]{constant}, 1);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithIntegralIntervalEvaluationsNoLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.integrateInterval(startX, endX);

            final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, value, 1);
            eval.setConstants(new double[]{0.0});
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polyParams, polynomial2.getPolyParams(), ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithIntegralIntervalEvaluationsLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.integrateInterval(startX, endX);

            final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, value, 1);
            eval.setConstants(new double[]{0.0});
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);


        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectEvaluationsNoLMSEAllowedSecondDegree() throws LockedException, NotReadyException,
            PolynomialEstimationException {

        final var estimator = new LMSEPolynomialEstimator(2);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 2nd degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectEvaluationsLMSEAllowedSecondDegree() throws LockedException, NotReadyException,
            PolynomialEstimationException {

        final var estimator = new LMSEPolynomialEstimator(2);
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 2nd degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polyParams, polynomial2.getPolyParams(), ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectAndSecondOrderDerivativeEvaluationsNoLMSEAllowedSecondDegree() throws LockedException,
            NotReadyException, PolynomialEstimationException {

        final var estimator = new LMSEPolynomialEstimator(2);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 2nd degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations() - 2; i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        var value = polynomial.evaluateDerivative(x);
        var eval = new DerivativePolynomialEvaluation(x, value, 1);
        evaluations.add(eval);

        x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        value = polynomial.evaluateSecondDerivative(x);
        eval = new DerivativePolynomialEvaluation(x, value, 2);
        evaluations.add(eval);

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);


        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithDirectAndSecondOrderDerivativeEvaluationLMSEAllowedSecondDegree() throws LockedException,
            NotReadyException, PolynomialEstimationException {

        final var estimator = new LMSEPolynomialEstimator(2);
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 2nd degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
        }

        var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        var value = polynomial.evaluateDerivative(x);
        var eval = new DerivativePolynomialEvaluation(x, value, 1);
        evaluations.add(eval);

        x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        value = polynomial.evaluateSecondDerivative(x);
        eval = new DerivativePolynomialEvaluation(x, value, 2);
        evaluations.add(eval);

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithSecondOrderIntegralEvaluationsNoLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var integral = polynomial.nthIntegrationAndReturnNew(2, constants);
            final var value = integral.evaluate(x);

            final var eval = new IntegralPolynomialEvaluation(x, value, constants, 2);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithSecondOrderIntegralEvaluationLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var integral = polynomial.nthIntegrationAndReturnNew(2, constants);
            final var value = integral.evaluate(x);

            final var eval = new IntegralPolynomialEvaluation(x, value, constants, 2);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);


        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithSecondOrderIntegralIntervalEvaluationsNoLMSEAllowed() throws LockedException,
            NotReadyException, PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var value = polynomial.nthOrderIntegrateInterval(startX, endX, 2, constants);

            final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, value, 2);
            eval.setConstants(constants);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    void testEstimateWithSecondOrderIntegralIntervalEvaluationLMSEAllowed() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var value = polynomial.nthOrderIntegrateInterval(startX, endX, 2, constants);

            final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, value, 2);
            eval.setConstants(constants);
            evaluations.add(eval);
        }

        estimator.setEvaluations(evaluations);

        assertTrue(estimator.isReady());
        assertSame(polyParams, polynomial.getPolyParams());

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final var polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polyParams, polynomial2.getPolyParams(), ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }


    private void reset() {
        estimateStart = estimateEnd = 0;
    }

    @Override
    public void onEstimateStart(final PolynomialEstimator estimator) {
        estimateStart++;
        checkIsLocked(estimator);
    }

    @Override
    public void onEstimateEnd(final PolynomialEstimator estimator) {
        estimateEnd++;
        checkIsLocked(estimator);
    }

    private void checkIsLocked(final PolynomialEstimator estimator) {
        assertTrue(estimator.isLocked());

        // Force LockedException
        assertThrows(LockedException.class, () -> estimator.setDegree(2));
        assertThrows(LockedException.class, () -> estimator.setEvaluations(null));
        assertThrows(LockedException.class, () -> estimator.setListener(null));
        assertThrows(LockedException.class, estimator::estimate);
        assertThrows(LockedException.class, () -> ((LMSEPolynomialEstimator) estimator).setLMSESolutionAllowed(true));
    }
}
