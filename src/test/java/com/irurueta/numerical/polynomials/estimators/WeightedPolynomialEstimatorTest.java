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
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

public class WeightedPolynomialEstimatorTest implements
        PolynomialEstimatorListener {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final int MIN_DEGREE = 1;
    private static final int MAX_DEGREE = 5;

    private static final double ABSOLUTE_ERROR = 1e-8;

    private int estimateStart;
    private int estimateEnd;

    @Test
    public void testConstructor() {
        // empty constructor
        WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());


        // test constructor with degree
        estimator = new WeightedPolynomialEstimator(2);

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with evaluations and weights
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        final double[] weights = new double[1];
        estimator = new WeightedPolynomialEstimator(evaluations, weights);

        // check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());


        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations, null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(wrongEvaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations,
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with listener
        estimator = new WeightedPolynomialEstimator(this);

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());


        // test constructor with degree, evaluations and weights
        estimator = new WeightedPolynomialEstimator(2, evaluations, weights);

        // check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());


        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(2,
                    null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations,
                    null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(2, wrongEvaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations,
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test degree and listener
        estimator = new WeightedPolynomialEstimator(2, this);

        // check correctness
        assertNull(estimator.getWeights());
        assertFalse(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test evaluations, weights and listener
        estimator = new WeightedPolynomialEstimator(evaluations, weights, this);

        // check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());


        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(
                    null, weights, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations,
                    null, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(wrongEvaluations, weights,
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(evaluations,
                    new double[2], this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with degree, evaluations, weights and listener
        estimator = new WeightedPolynomialEstimator(2, evaluations, weights,
                this);

        // check correctness
        assertSame(estimator.getWeights(), weights);
        assertTrue(estimator.areWeightsAvailable());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());


        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new WeightedPolynomialEstimator(2,
                    null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations,
                    null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(2, wrongEvaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new WeightedPolynomialEstimator(2, evaluations,
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testGetSetMaxEvaluations() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default value
        assertEquals(WeightedPolynomialEstimator.DEFAULT_MAX_EVALUATIONS,
                estimator.getMaxEvaluations());

        // set new value
        estimator.setMaxEvaluations(100);

        // check correctness
        assertEquals(100, estimator.getMaxEvaluations());

        // Force IllegalArgumentException
        try {
            estimator.setMaxEvaluations(1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsSetSortWeightsEnabled() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default value
        assertEquals(WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS,
                estimator.isSortWeightsEnabled());

        // set new value
        estimator.setSortWeightsEnabled(
                !WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);

        // check correctness
        assertEquals(estimator.isSortWeightsEnabled(),
                !WeightedPolynomialEstimator.DEFAULT_SORT_WEIGHTS);
    }

    @Test
    public void testGetSetEvaluationsAndWeights() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        assertNull(estimator.getWeights());
        assertNull(estimator.getEvaluations());

        // set new values
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        final double[] weights = new double[1];
        estimator.setEvaluationsAndWeights(evaluations, weights);

        // check correctness
        assertSame(estimator.getEvaluations(), evaluations);
        assertSame(estimator.getWeights(), weights);

        // Force IllegalArgumentException
        try {
            estimator.setEvaluations(evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setDegreeAndEvaluations(2, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        try {
            estimator.setEvaluationsAndWeights(null,
                    weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setEvaluationsAndWeights(evaluations, null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setEvaluationsAndWeights(wrongEvaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setEvaluationsAndWeights(evaluations, new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetDegree() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default value
        assertEquals(1, estimator.getDegree());

        // set new value
        estimator.setDegree(2);

        // check correctness
        assertEquals(2, estimator.getDegree());

        // Force IllegalArgumentException
        try {
            estimator.setDegree(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetDegreeEvaluationsAndWeights() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        assertEquals(WeightedPolynomialEstimator.MIN_DEGREE,
                estimator.getDegree());
        assertNull(estimator.getWeights());
        assertNull(estimator.getEvaluations());

        // set new values
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        final double[] weights = new double[1];
        estimator.setDegreeEvaluationsAndWeights(2, evaluations, weights);

        // check correctness
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertSame(estimator.getWeights(), weights);

        // Force IllegalArgumentException
        try {
            estimator.setEvaluations(evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setDegreeAndEvaluations(2, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setDegreeEvaluationsAndWeights(0, evaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        try {
            estimator.setDegreeEvaluationsAndWeights(2,
                    null, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setDegreeEvaluationsAndWeights(2, evaluations,
                    null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setDegreeEvaluationsAndWeights(2, wrongEvaluations, weights);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setDegreeEvaluationsAndWeights(2, evaluations,
                    new double[2]);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsReady() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default value
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());

        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights = new double[2];
        for (int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double value = polynomial.evaluate(x);

            final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                    value);
            evaluations.add(eval);
            weights[i] = 1.0;
        }

        estimator.setEvaluationsAndWeights(evaluations, weights);

        assertTrue(estimator.isReady());
    }

    @Test
    public void testGetMinNumberOfEvaluations() {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int degree = randomizer.nextInt(MIN_DEGREE, MAX_DEGREE);
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(degree),
                degree + 1);

        // Force IllegalArgumentException
        try {
            assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(0),
                    degree + 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator(degree);
        assertEquals(estimator.getMinNumberOfEvaluations(), degree + 1);
    }

    @Test
    public void testGetSetListener() throws LockedException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default value
        assertNull(estimator.getListener());

        // set new value
        estimator.setListener(this);

        // check correctness
        assertSame(estimator.getListener(), this);
    }

    @Test
    public void testEstimateWithDirectEvaluations()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double value = polynomial.evaluate(x);

            final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                    value);
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
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    public void testEstimateWithDirectAndDerivativeEvaluations()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[4 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double value = polynomial.evaluate(x);

            final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                    value);
            evaluations.add(eval);

            weights[i] = randomizer.nextDouble(0.5, 1.0);
        }
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double value = polynomial.evaluateDerivative(x);
            final DerivativePolynomialEvaluation eval =
                    new DerivativePolynomialEvaluation(x, value, 1);
            evaluations.add(eval);
            weights[2 * estimator.getMinNumberOfEvaluations() + i] =
                    randomizer.nextDouble(0.5, 1.0);
        }


        estimator.setEvaluationsAndWeights(evaluations, weights);

        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);

        estimator.setListener(this);
        reset();

        assertEquals(0, estimateStart);
        assertEquals(0, estimateEnd);

        // estimate
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    public void testEstimateWithIntegralEvaluations()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double constant = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final Polynomial integral = polynomial.integrationAndReturnNew(constant);
            final double value = integral.evaluate(x);

            final IntegralPolynomialEvaluation eval =
                    new IntegralPolynomialEvaluation(x, value,
                            new double[]{constant}, 1);
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
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    public void testEstimateWithIntegralIntervalEvaluations()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double startX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double value = polynomial.integrateInterval(startX, endX);

            final IntegralIntervalPolynomialEvaluation eval =
                    new IntegralIntervalPolynomialEvaluation(startX, endX,
                            value, 1);
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
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    public void testEstimateWithDirectEvaluationsSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator(2);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double value = polynomial.evaluate(x);

            final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                    value);
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
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    public void testEstimateWithSecondOrderIntegralEvaluationsSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator(2);

        // check default values
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(2, polynomial.getDegree());

        assertEquals(3, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final Polynomial integral = polynomial.nthIntegrationAndReturnNew(2,
                    constants);
            final double value = integral.evaluate(x);

            final IntegralPolynomialEvaluation eval =
                    new IntegralPolynomialEvaluation(x, value, constants, 2);
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
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
        assertEquals(1, estimateStart);
        assertEquals(1, estimateEnd);
    }

    @Test
    public void testEstimateWithSecondOrderIntegralIntervalEvaluationsSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        final WeightedPolynomialEstimator estimator =
                new WeightedPolynomialEstimator();

        // check default values
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertEquals(50, estimator.getMaxEvaluations());

        // Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }


        // create random 1st degree polynomial
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final Polynomial polynomial = new Polynomial(polyParams);

        assertEquals(1, polynomial.getDegree());

        assertEquals(2, estimator.getMinNumberOfEvaluations());
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        final double[] weights =
                new double[2 * estimator.getMinNumberOfEvaluations()];
        for (int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            final double startX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            final double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final double value = polynomial.nthOrderIntegrateInterval(startX, endX, 2,
                    constants);

            final IntegralIntervalPolynomialEvaluation eval =
                    new IntegralIntervalPolynomialEvaluation(startX, endX,
                            value, 2);
            eval.setConstants(constants);
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
        final Polynomial polynomial2 = estimator.estimate();

        // check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                ABSOLUTE_ERROR);
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
        try {
            estimator.setDegree(2);
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        }
        try {
            ((WeightedPolynomialEstimator) estimator).setEvaluationsAndWeights(
                    null, null);
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        }
        try {
            ((WeightedPolynomialEstimator) estimator).
                    setDegreeEvaluationsAndWeights(2, null, null);
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        }
        try {
            ((WeightedPolynomialEstimator) estimator).setMaxEvaluations(1);
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        }
        try {
            ((WeightedPolynomialEstimator) estimator).setSortWeightsEnabled(
                    false);
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        }
        try {
            estimator.setListener(null);
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        }
        try {
            estimator.estimate();
            fail("LockedException expected but not thrown");
        } catch (final LockedException ignore) {
        } catch (final Exception ignore) {
            fail("LockedException expected but not thrown");
        }
    }

    private void reset() {
        estimateStart = estimateEnd = 0;
    }
}
