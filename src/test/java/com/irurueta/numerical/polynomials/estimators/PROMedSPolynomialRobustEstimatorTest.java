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
import org.junit.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

public class PROMedSPolynomialRobustEstimatorTest implements
        PolynomialRobustEstimatorListener {

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
    public void testConstructor() {
        // test empty constructor
        PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());


        // test constructor with degree
        estimator = new PROMedSPolynomialRobustEstimator(2);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROMedSPolynomialRobustEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test constructor with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new PROMedSPolynomialRobustEstimator(evaluations);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        estimator = null;
        try {
            estimator = new PROMedSPolynomialRobustEstimator(wrongEvaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test constructor with listener
        estimator = new PROMedSPolynomialRobustEstimator(this);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // test constructor with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new PROMedSPolynomialRobustEstimator(2, evaluations);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROMedSPolynomialRobustEstimator(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new PROMedSPolynomialRobustEstimator(2, wrongEvaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test constructor with degree and listener
        estimator = new PROMedSPolynomialRobustEstimator(2, this);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROMedSPolynomialRobustEstimator(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test constructor with evaluations and listener
        estimator = new PROMedSPolynomialRobustEstimator(evaluations, this);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROMedSPolynomialRobustEstimator(wrongEvaluations, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test constructor with degree, evaluations and listener
        estimator = new PROMedSPolynomialRobustEstimator(2, evaluations, this);

        // check correctness
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROMedSPolynomialRobustEstimator(0, evaluations,
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new PROMedSPolynomialRobustEstimator(2, wrongEvaluations,
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        assertNull(estimator);
    }

    @Test
    public void testGetSetThreshold() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertEquals(PROMedSPolynomialRobustEstimator.DEFAULT_STOP_THRESHOLD,
                estimator.getStopThreshold(), 0.0);

        // set new value
        estimator.setStopThreshold(1.0);

        // check correctness
        assertEquals(1.0, estimator.getStopThreshold(), 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setStopThreshold(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetEvaluations() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getEvaluations());

        // set new value
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator.setEvaluations(evaluations);

        // check correctness
        assertSame(estimator.getEvaluations(), evaluations);

        // Force IllegalArgumentException
        try {
            estimator.setEvaluations(null);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setEvaluations(new ArrayList<PolynomialEvaluation>());
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetListener() {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getListener());

        // set new value
        estimator.setListener(this);

        // check correctness
        assertSame(estimator.getListener(), this);
    }

    @Test
    public void testGetSetProgressDelta() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);

        // set new value
        estimator.setProgressDelta(0.5f);

        // check correctness
        assertEquals(0.5, estimator.getProgressDelta(), 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setProgressDelta(-1.0f);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setProgressDelta(2.0f);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetConfidence() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);

        // set new value
        estimator.setConfidence(0.5);

        // check correctness
        assertEquals(0.5, estimator.getConfidence(), 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setConfidence(-1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator.setConfidence(2.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetMaxIterations() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());

        // set new value
        estimator.setMaxIterations(10);

        // check correctness
        assertEquals(10, estimator.getMaxIterations());

        // Force IllegalArgumentException
        try {
            estimator.setMaxIterations(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsSetGeometricDistanceUsed() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());

        // set new value
        estimator.setGeometricDistanceUsed(
                !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);

        // check correctness
        assertEquals(estimator.isGeometricDistanceUsed(),
                !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
    }

    @Test
    public void testGetSetDegree() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());

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
    public void testGetSetQualityScores() throws LockedException {
        final PROMedSPolynomialRobustEstimator estimator =
                new PROMedSPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getQualityScores());

        // set new value
        final double[] scores = new double[2];
        estimator.setQualityScores(scores);

        // check correctness
        assertSame(estimator.getQualityScores(), scores);

        // Force IllegalArgumentException
        final double[] wrong = new double[1];
        try {
            estimator.setQualityScores(wrong);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testEstimateDirectEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROMedSPolynomialRobustEstimator estimator =
                    new PROMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

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

            final int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS,
                    MAX_EVALUATIONS);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            final List<PolynomialEvaluation> evaluations = new ArrayList<>();
            final double[] qualityScores = new double[numEvaluations];
            for (int i = 0; i < numEvaluations; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                        valueWithError);
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
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateDirectAndDerivativeEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROMedSPolynomialRobustEstimator estimator =
                    new PROMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

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

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS,
                    MAX_EVALUATIONS);
            if (numEvaluations % 2 != 0) {
                numEvaluations++;
            }

            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            final List<PolynomialEvaluation> evaluations = new ArrayList<>();
            final double[] qualityScores = new double[numEvaluations];
            int j = 0;
            for (int i = 0; i < numEvaluations / 2; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[j] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                        valueWithError);
                evaluations.add(eval);
                j++;
            }
            for (int i = 0; i < numEvaluations / 2; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.evaluateDerivative(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[j] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final DerivativePolynomialEvaluation eval =
                        new DerivativePolynomialEvaluation(x, valueWithError,
                                1);
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
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateIntegralEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROMedSPolynomialRobustEstimator estimator =
                    new PROMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

            // Force NotReadyException
            try {
                estimator.estimate();
                fail("NotReadyException expected but not thrown");
            } catch (NotReadyException ignore) {
            }

            // create random 1st degree polynomial
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());
            final double[] polyParams = new double[2];
            randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

            final Polynomial polynomial = new Polynomial(polyParams);

            final int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS,
                    MAX_EVALUATIONS);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            final List<PolynomialEvaluation> evaluations = new ArrayList<>();
            final double[] qualityScores = new double[numEvaluations];
            for (int i = 0; i < numEvaluations; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double constant = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final Polynomial integral = polynomial.integrationAndReturnNew(
                        constant);
                final double value = integral.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final IntegralPolynomialEvaluation eval =
                        new IntegralPolynomialEvaluation(x, valueWithError,
                                new double[]{constant}, 1);
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
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateIntegralIntervalEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROMedSPolynomialRobustEstimator estimator =
                    new PROMedSPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertFalse(estimator.isGeometricDistanceUsed());

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

            final int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS,
                    MAX_EVALUATIONS);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            final List<PolynomialEvaluation> evaluations = new ArrayList<>();
            final double[] qualityScores = new double[numEvaluations];
            for (int i = 0; i < numEvaluations; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final double startX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double endX = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.integrateInterval(startX, endX);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final IntegralIntervalPolynomialEvaluation eval =
                        new IntegralIntervalPolynomialEvaluation(startX, endX,
                                valueWithError, 1);
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
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateDirectEvaluationsGeometricDistance()
            throws LockedException, NotReadyException, RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROMedSPolynomialRobustEstimator estimator =
                    new PROMedSPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertTrue(estimator.isGeometricDistanceUsed());

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

            final int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS,
                    MAX_EVALUATIONS);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            final List<PolynomialEvaluation> evaluations = new ArrayList<>();
            final double[] qualityScores = new double[numEvaluations];
            for (int i = 0; i < numEvaluations; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[i] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                        valueWithError);
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
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateDirectAndDerivativeEvaluationsGeometricDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROMedSPolynomialRobustEstimator estimator =
                    new PROMedSPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            // check default values
            assertEquals(1, estimator.getDegree());
            assertFalse(estimator.isReady());
            assertTrue(estimator.isGeometricDistanceUsed());

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

            int numEvaluations = randomizer.nextInt(MIN_EVALUATIONS,
                    MAX_EVALUATIONS);
            if (numEvaluations % 2 != 0) {
                numEvaluations++;
            }

            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, STD_ERROR);
            final List<PolynomialEvaluation> evaluations = new ArrayList<>();
            final double[] qualityScores = new double[numEvaluations];
            int j = 0;
            for (int i = 0; i < numEvaluations / 2; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.evaluate(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[i] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x,
                        valueWithError);
                evaluations.add(eval);
                j++;
            }
            for (int i = 0; i < numEvaluations / 2; i++) {
                final double scoreError = randomizer.nextDouble(MIN_SCORE_ERROR,
                        MAX_SCORE_ERROR);
                qualityScores[j] = 1.0 + scoreError;

                final double x = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double value = polynomial.evaluateDerivative(x);

                final double valueWithError;
                if (randomizer.nextInt(0, 100) < PERCENTAGE_OUTLIER) {
                    // evaluation is outlier
                    final double error = errorRandomizer.nextDouble();
                    valueWithError = value + error;
                    qualityScores[j] = 1.0 / (1.0 + Math.abs(error)) +
                            scoreError;
                } else {
                    valueWithError = value;
                }

                final DerivativePolynomialEvaluation eval =
                        new DerivativePolynomialEvaluation(x, valueWithError,
                                1);
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
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(1, estimateStart);
            assertEquals(1, estimateEnd);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    private void reset() {
        estimateStart = estimateEnd = estimateNextIteration =
                estimateProgressChange = 0;
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
    public void onEstimateNextIteration(final PolynomialRobustEstimator estimator,
                                        final int iteration) {
        estimateNextIteration++;
    }

    @Override
    public void onEstimateProgressChange(final PolynomialRobustEstimator estimator,
                                         final float progress) {
        estimateProgressChange++;
    }
}
