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
import com.irurueta.numerical.robust.RobustEstimatorException;
import com.irurueta.numerical.robust.RobustEstimatorMethod;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.*;

class PROSACPolynomialRobustEstimatorTest implements PolynomialRobustEstimatorListener {

    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;

    private static final double ABSOLUTE_ERROR = 1e-8;

    private static final int PERCENTAGE_OUTLIER = 20;

    private static final int MIN_EVALUATIONS = 500;
    private static final int MAX_EVALUATIONS = 1000;

    private static final double STD_ERROR = 100.0;

    private static final double MIN_SCORE_ERROR = -0.3;
    private static final double MAX_SCORE_ERROR = 0.3;

    private static final int TIMES = 10;

    private int estimateStart;
    private int estimateEnd;
    private int estimateNextIteration;
    private int estimateProgressChange;

    @Test
    void testConstructor() {
        // test empty constructor
        var estimator = new PROSACPolynomialRobustEstimator();

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(PolynomialEstimator.MIN_DEGREE),
                estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // test constructor with degree
        estimator = new PROSACPolynomialRobustEstimator(2);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(2), estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(0));

        // test constructor with evaluations
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new PROSACPolynomialRobustEstimator(evaluations);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(PolynomialEstimator.MIN_DEGREE),
                estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(wrongEvaluations));

        // test constructor with listener
        estimator = new PROSACPolynomialRobustEstimator(this);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(PolynomialEstimator.MIN_DEGREE),
                estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // test constructor with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new PROSACPolynomialRobustEstimator(2, evaluations);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(2), estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(0, evaluations));
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(2,
                wrongEvaluations));

        // test constructor with degree and listener
        estimator = new PROSACPolynomialRobustEstimator(2, this);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(2), estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(0, this));

        // test constructor with evaluations and listener
        estimator = new PROSACPolynomialRobustEstimator(evaluations, this);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(PolynomialEstimator.MIN_DEGREE),
                estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(wrongEvaluations,
                this));

        // test constructor with degree, evaluations and listener
        estimator = new PROSACPolynomialRobustEstimator(2, evaluations, this);

        // check correctness
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(2), estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(0, evaluations,
                this));
        assertThrows(IllegalArgumentException.class, () -> new PROSACPolynomialRobustEstimator(2,
                wrongEvaluations, this));
    }

    @Test
    void testGetSetThreshold() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, estimator.getThreshold(), 0.0);

        // set new value
        estimator.setThreshold(1.0);

        // check correctness
        assertEquals(1.0, estimator.getThreshold(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setThreshold(0.0));
    }

    @Test
    void testGetSetEvaluations() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getEvaluations());

        // set new value
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator.setEvaluations(evaluations);

        // check correctness
        assertSame(estimator.getEvaluations(), evaluations);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluations(null));
        final var wrong = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> estimator.setEvaluations(wrong));
    }

    @Test
    void testGetSetListener() {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getListener());

        // set new value
        estimator.setListener(this);

        // check correctness
        assertSame(this, estimator.getListener());
    }

    @Test
    void testGetSetProgressDelta() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);

        // set new value
        estimator.setProgressDelta(0.5f);

        // check correctness
        assertEquals(0.5, estimator.getProgressDelta(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setProgressDelta(-1.0f));
        assertThrows(IllegalArgumentException.class, () -> estimator.setProgressDelta(2.0f));
    }

    @Test
    void testGetSetConfidence() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);

        // set new value
        estimator.setConfidence(0.5);

        // check correctness
        assertEquals(0.5, estimator.getConfidence(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setConfidence(-1.0));
        assertThrows(IllegalArgumentException.class, () -> estimator.setConfidence(2.0));
    }

    @Test
    void testGetSetMaxIterations() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());

        // set new value
        estimator.setMaxIterations(10);

        // check correctness
        assertEquals(10, estimator.getMaxIterations());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setMaxIterations(0));
    }

    @Test
    void testIsSetGeometricDistanceUsed() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());

        // set new value
        estimator.setGeometricDistanceUsed(!PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);

        //check correctness
        assertEquals(estimator.isGeometricDistanceUsed(), !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
    }

    @Test
    void testGetSetDegree() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());

        // set new value
        estimator.setDegree(2);

        // check correctness
        assertEquals(2, estimator.getDegree());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setDegree(0));
    }

    @Test
    void testGetSetQualityScores() throws LockedException {
        final var estimator = new PROSACPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getQualityScores());

        // set new value
        final var scores = new double[2];
        estimator.setQualityScores(scores);

        // check correctness
        assertSame(scores, estimator.getQualityScores());

        // Force IllegalArgumentException
        final var wrong = new double[1];
        assertThrows(IllegalArgumentException.class, () -> estimator.setQualityScores(wrong));
    }

    @Test
    void testEstimateDirectEvaluationsAlgebraicDistance() throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (var t = 0; t < TIMES; t++) {
            final var estimator = new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // create random 1st degree polynomial
            final var randomizer = new UniformRandomizer();
            final var polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var polynomial = new Polynomial(polyParams);

            final var numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, MAX_EVALUATIONS);
            final var errorRandomizer = new GaussianRandomizer(0.0, STD_ERROR);
            final var evaluations = new ArrayList<PolynomialEvaluation>();
            final var qualityScores = new double[numEvaluations];
            for (var i = 0; i < numEvaluations; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new DirectPolynomialEvaluation(x, valueWithError);
                evaluations.add(eval);
            }

            estimator.setEvaluations(evaluations);
            estimator.setQualityScores(qualityScores);

            estimator.setListener(this);
            reset();

            assertEquals(0, estimateStart);
            assertEquals(0, estimateEnd);
            assertEquals(0, estimateNextIteration);
            assertEquals(0, estimateProgressChange);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final var polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polyParams, polynomial2.getPolyParams(), ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    void testEstimateDirectAndDerivativeEvaluationsAlgebraicDistance() throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (var t = 0; t < TIMES; t++) {
            final var estimator = new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // create random 1st degree polynomial
            final var randomizer = new UniformRandomizer();
            final var polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var polynomial = new Polynomial(polyParams);

            var numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, MAX_EVALUATIONS);
            if (numEvaluations % 2 != 0) {
                numEvaluations++;
            }

            final var errorRandomizer = new GaussianRandomizer(0.0, STD_ERROR);
            final var evaluations = new ArrayList<PolynomialEvaluation>();
            final var qualityScores = new double[numEvaluations];
            var j = 0;
            for (var i = 0; i < numEvaluations / 2; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[j] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new DirectPolynomialEvaluation(x, valueWithError);
                evaluations.add(eval);
                j++;
            }
            for (var i = 0; i < numEvaluations / 2; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.evaluateDerivative(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[j] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new DerivativePolynomialEvaluation(x, valueWithError, 1);
                evaluations.add(eval);
                j++;
            }

            estimator.setEvaluations(evaluations);
            estimator.setQualityScores(qualityScores);

            estimator.setListener(this);
            reset();

            assertEquals(0, estimateStart);
            assertEquals(0, estimateEnd);
            assertEquals(0, estimateNextIteration);
            assertEquals(0, estimateProgressChange);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final var polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    void testEstimateIntegralEvaluationsAlgebraicDistance() throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (var t = 0; t < TIMES; t++) {
            final var estimator = new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // create random 1st degree polynomial
            final var randomizer = new UniformRandomizer();
            final var polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var polynomial = new Polynomial(polyParams);

            final var numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, MAX_EVALUATIONS);
            final var errorRandomizer = new GaussianRandomizer(0.0, STD_ERROR);
            final var evaluations = new ArrayList<PolynomialEvaluation>();
            final var qualityScores = new double[numEvaluations];
            for (var i = 0; i < numEvaluations; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var constant = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var integral = polynomial.integrationAndReturnNew(constant);
                final var value = integral.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new IntegralPolynomialEvaluation(x, valueWithError, new double[]{constant},
                        1);
                evaluations.add(eval);
            }

            estimator.setEvaluations(evaluations);
            estimator.setQualityScores(qualityScores);

            estimator.setListener(this);
            reset();

            assertEquals(0, estimateStart);
            assertEquals(0, estimateEnd);
            assertEquals(0, estimateNextIteration);
            assertEquals(0, estimateProgressChange);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final var polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    void testEstimateIntegralIntervalEvaluationsAlgebraicDistance() throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (var t = 0; t < TIMES; t++) {
            final var estimator = new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // create random 1st degree polynomial
            final var randomizer = new UniformRandomizer();
            final var polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var polynomial = new Polynomial(polyParams);

            final var numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, MAX_EVALUATIONS);
            final var errorRandomizer = new GaussianRandomizer(0.0, STD_ERROR);
            final var evaluations = new ArrayList<PolynomialEvaluation>();
            final var qualityScores = new double[numEvaluations];
            for (var i = 0; i < numEvaluations; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final var startX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.integrateInterval(startX, endX);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new IntegralIntervalPolynomialEvaluation(startX, endX, valueWithError, 1);
                evaluations.add(eval);
            }

            estimator.setEvaluations(evaluations);
            estimator.setQualityScores(qualityScores);

            estimator.setListener(this);
            reset();

            assertEquals(0, estimateStart);
            assertEquals(0, estimateEnd);
            assertEquals(0, estimateNextIteration);
            assertEquals(0, estimateProgressChange);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final var polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams, ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    void testEstimateDirectEvaluationsGeometricDistance() throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (var t = 0; t < TIMES; t++) {
            final var estimator = new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertTrue(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // create random 1st degree polynomial
            final var randomizer = new UniformRandomizer();
            final var polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var polynomial = new Polynomial(polyParams);

            final var numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, MAX_EVALUATIONS);
            final var errorRandomizer = new GaussianRandomizer(0.0, STD_ERROR);
            final var evaluations = new ArrayList<PolynomialEvaluation>();
            final var qualityScores = new double[numEvaluations];
            for (var i = 0; i < numEvaluations; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new DirectPolynomialEvaluation(x, valueWithError);
                evaluations.add(eval);
            }

            estimator.setEvaluations(evaluations);
            estimator.setQualityScores(qualityScores);

            estimator.setListener(this);
            reset();

            assertEquals(0, estimateStart);
            assertEquals(0, estimateEnd);
            assertEquals(0, estimateNextIteration);
            assertEquals(0, estimateProgressChange);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final var polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polyParams, polynomial2.getPolyParams(), ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    void testEstimateDirectAndDerivativeEvaluationsGeometricDistance() throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (var t = 0; t < TIMES; t++) {
            final var estimator = new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertTrue(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            assertThrows(NotReadyException.class, estimator::estimate);

            // create random 1st degree polynomial
            final var randomizer = new UniformRandomizer();
            final var polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final var polynomial = new Polynomial(polyParams);

            var numEvaluations = randomizer.nextInt(MIN_EVALUATIONS, MAX_EVALUATIONS);
            if (numEvaluations % 2 != 0) {
                numEvaluations++;
            }

            final var errorRandomizer = new GaussianRandomizer(0.0, STD_ERROR);
            final var evaluations = new ArrayList<PolynomialEvaluation>();
            final var qualityScores = new double[numEvaluations];
            var j = 0;
            for (var i = 0; i < numEvaluations / 2; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new DirectPolynomialEvaluation(x, valueWithError);
                evaluations.add(eval);
                j++;
            }
            for (var i = 0; i < numEvaluations / 2; i++) {
                final var scoreError = randomizer.nextDouble(MIN_SCORE_ERROR, MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final var x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var value = polynomial.evaluateDerivative(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final var error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[j] = 1.0 / (1.0 + Math.abs(error)) + scoreError;
                } else {
                    valueWithError = value;
                }

                final var eval = new DerivativePolynomialEvaluation(x, valueWithError, 1);
                evaluations.add(eval);
                j++;
            }

            estimator.setEvaluations(evaluations);
            estimator.setQualityScores(qualityScores);

            estimator.setListener(this);
            reset();

            assertEquals(0, estimateStart);
            assertEquals(0, estimateEnd);
            assertEquals(0, estimateNextIteration);
            assertEquals(0, estimateProgressChange);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final var polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polyParams, polynomial2.getPolyParams(), ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    private void reset() {
        estimateStart = estimateEnd = estimateNextIteration = estimateProgressChange = 0;
    }

    @Override
    public void onEstimateStart(final PolynomialRobustEstimator estimator) {
        estimateStart++;
    }

    @Override
    public void onEstimateEnd(final PolynomialRobustEstimator estimator) {
        estimateEnd++;
    }

    @Override
    public void onEstimateNextIteration(final PolynomialRobustEstimator estimator, final int iteration) {
        estimateNextIteration++;
    }

    @Override
    public void onEstimateProgressChange(final PolynomialRobustEstimator estimator, final float progress) {
        estimateProgressChange++;
    }
}
