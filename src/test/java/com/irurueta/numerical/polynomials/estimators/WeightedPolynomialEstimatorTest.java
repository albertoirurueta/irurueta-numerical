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

class WeightedPolynomialEstimatorTest implements PolynomialEstimatorListener {

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
        var estimator = new WeightedPolynomialEstimator();

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // test constructor with degree
        estimator = new WeightedPolynomialEstimator(2);

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(0));

        // test constructor with evaluations and weights
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        final var weights = new double[1];
        estimator = new WeightedPolynomialEstimator(evaluations, weights);

        // check correctness
        assertSame(weights, estimator.getWeights());
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(null, weights));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(evaluations, null));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(wrongEvaluations, weights));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(evaluations, new double[2]));

        // test constructor with listener
        estimator = new WeightedPolynomialEstimator(this);

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // test constructor with degree, evaluations and weights
        estimator = new WeightedPolynomialEstimator(2, evaluations, weights);

        // check correctness
        assertSame(weights, estimator.getWeights());
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, null,
                weights));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, evaluations,
                null));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, wrongEvaluations,
                weights));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, evaluations,
                new double[2]));

        // test degree and listener
        estimator = new WeightedPolynomialEstimator(2, this);

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(0, this));

        // test evaluations, weights and listener
        estimator = new WeightedPolynomialEstimator(evaluations, weights, this);

        // check correctness
        assertSame(weights, estimator.getWeights());
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(null, weights,
                this));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(evaluations, null,
                this));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(wrongEvaluations, weights,
                this));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(evaluations, new double[2],
                this));

        // test constructor with degree, evaluations, weights and listener
        estimator = new WeightedPolynomialEstimator(2, evaluations, weights, this);

        // check correctness
        assertSame(weights, estimator.getWeights());
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(this, estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR, estimator.getType());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, null,
                weights));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, evaluations,
                null));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, wrongEvaluations,
                weights));
        assertThrows(IllegalArgumentException.class, () -> new WeightedPolynomialEstimator(2, evaluations,
                new double[2]));
    }

    @Test
    void testGetSetMaxEvaluations() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default value
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS, estimator.getMaxEvaluations());

        // set new value
        estimator.setMaxEvaluations(100);

        // check correctness
        assertEquals(100, estimator.getMaxEvaluations());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setMaxEvaluations(1));
    }

    @Test
    void testIsSetSortWeightsEnabled() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default value
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS, estimator.isSortWeightsEnabled());

        // set new value
        estimator.setSortWeightsEnabled(!WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);

        // check correctness
        assertEquals(estimator.isSortWeightsEnabled(), !WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
    }

    @Test
    void testGetSetEvaluationsAndWeights() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

        assertNull(estimator.getWeights());
        assertNull(estimator.getEvaluations());

        // set new values
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        final var weights = new double[1];
        estimator.setEvaluationsAndWeights(evaluations, weights);

        // check correctness
        assertSame(evaluations, estimator.getEvaluations());
        assertSame(weights, estimator.getWeights());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluations(evaluations));
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeAndEvaluations(2, evaluations));

        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluationsAndWeights(null,
                weights));
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluationsAndWeights(evaluations,
                null));
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluationsAndWeights(wrongEvaluations,
                weights));
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluationsAndWeights(evaluations,
                new double[2]));
    }

    @Test
    void testGetSetDegree() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

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
    void testGetSetDegreeEvaluationsAndWeights() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

        assertEquals(WeightedPolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertNull(estimator.getWeights());
        assertNull(estimator.getEvaluations());

        // set new values
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        final var weights = new double[1];
        estimator.setDegreeEvaluationsAndWeights(2, evaluations, weights);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertSame(evaluations, estimator.getEvaluations());
        assertSame(weights, estimator.getWeights());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluations(evaluations));
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeAndEvaluations(2, evaluations));
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeEvaluationsAndWeights(0,
                evaluations, weights));

        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeEvaluationsAndWeights(2,
                null, weights));
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeEvaluationsAndWeights(2,
                evaluations, null));
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeEvaluationsAndWeights(2,
                wrongEvaluations, weights));
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegreeEvaluationsAndWeights(2,
                evaluations, new double[2]));
    }

    @Test
    void testIsReady() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

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
        final var weights = new double[2];
        for (var i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);
            weights[i] = 1.0;
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

        assertTrue(estimator.isReady());
    }

    @Test
    void testGetMinNumberOfEvaluations() {
        final var randomizer = new UniformRandomizer();

        final var degree = randomizer.nextInt(MIN_DEGREE, MAX_DEGREE);
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(degree), degree + 1);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialEstimator.getMinNumberOfEvaluations(0));

        final var estimator = new WeightedPolynomialEstimator(degree);
        assertEquals(degree + 1, estimator.getMinNumberOfEvaluations());
    }

    @Test
    void testGetSetListener() throws LockedException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default value
        assertNull(estimator.getListener());

        // set new value
        estimator.setListener(this);

        // check correctness
        assertSame(this, estimator.getListener());
    }

    @Test
    void testEstimateWithDirectEvaluations() throws LockedException, NotReadyException, PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

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
        final var weights = new double[2 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
    void testEstimateWithDirectAndDerivativeEvaluations() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

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
        final var weights = new double[4 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluateDerivative(x);
            final var eval = new DerivativePolynomialEvaluation(x, value, 1);
            evaluations.add(eval);
            weights[2 * estimator.getMinNumberOfEvaluations() + i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
    void testEstimateWithIntegralEvaluations() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

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
        final var weights = new double[2 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var integral = polynomial.integrationAndReturnNew(constant);
            final var value = integral.evaluate(x);

            final var eval = new IntegralPolynomialEvaluation(x, value, new double[]{constant}, 1);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
    void testEstimateWithIntegralIntervalEvaluations() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

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
        final var weights = new double[2 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.integrateInterval(startX, endX);

            final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, value, 1);
            eval.setConstants(new double[]{0.0});
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
    void testEstimateWithDirectEvaluationsSecondDegree() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator(2);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        final var weights = new double[2 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var value = polynomial.evaluate(x);

            final var eval = new DirectPolynomialEvaluation(x, value);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
    void testEstimateWithSecondOrderIntegralEvaluationsSecondDegree() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator(2);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        assertThrows(NotReadyException.class, estimator::estimate);

        // create random 1st degree polynomial
        final var randomizer = new UniformRandomizer();
        final var polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        final var weights = new double[2 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var integral = polynomial.nthIntegrationAndReturnNew(2, constants);
            final var value = integral.evaluate(x);

            final var eval = new IntegralPolynomialEvaluation(x, value, constants, 2);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
    void testEstimateWithSecondOrderIntegralIntervalEvaluationsSecondDegree() throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final var estimator = new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

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
        final var weights = new double[2 * estimator.getMinNumberOfEvaluations()];
        for (var i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var value = polynomial.nthOrderIntegrateInterval(startX, endX, 2, constants);

            final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, value, 2);
            eval.setConstants(constants);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

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
        assertThrows(LockedException.class, () -> ((WeightedPolynomialEstimator) estimator).setEvaluationsAndWeights(
                null, null));
        assertThrows(LockedException.class, () -> ((WeightedPolynomialEstimator) estimator).
                setDegreeEvaluationsAndWeights(2, null, null));
        assertThrows(LockedException.class, () -> ((WeightedPolynomialEstimator) estimator).setMaxEvaluations(1));
        assertThrows(LockedException.class, () -> ((WeightedPolynomialEstimator) estimator).setSortWeightsEnabled(
                false));
        assertThrows(LockedException.class, () -> estimator.setListener(null));
        assertThrows(LockedException.class, estimator::estimate);
    }

    private void reset() {
        estimateStart = estimateEnd = 0;
    }
}
