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

@SuppressWarnings("Duplicates")
public class PROSACPolynomialRobustEstimatorTest implements
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
        PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());


        // test constructor with degree
        estimator = new PROSACPolynomialRobustEstimator(2);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROSACPolynomialRobustEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new PROSACPolynomialRobustEstimator(evaluations);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvals = new ArrayList<>();
        estimator = null;
        try {
            estimator = new PROSACPolynomialRobustEstimator(wrongEvals);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with listener
        estimator = new PROSACPolynomialRobustEstimator(this);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());


        // test constructor with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = new PROSACPolynomialRobustEstimator(2, evaluations);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROSACPolynomialRobustEstimator(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new PROSACPolynomialRobustEstimator(2, wrongEvals);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with degree and listener
        estimator = new PROSACPolynomialRobustEstimator(2, this);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertNull(estimator.getEvaluations());
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROSACPolynomialRobustEstimator(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with evaluations and listener
        estimator = new PROSACPolynomialRobustEstimator(evaluations, this);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(
                        PolynomialEstimator.MIN_DEGREE));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROSACPolynomialRobustEstimator(wrongEvals, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test constructor with degree, evaluations and listener
        estimator = new PROSACPolynomialRobustEstimator(2, evaluations, this);

        // check correctness
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);
        assertEquals(estimator.getMethod(), RobustEstimatorMethod.PROSAC);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(estimator.getMinNumberOfEvaluations(),
                PolynomialEstimator.getMinNumberOfEvaluations(2));
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new PROSACPolynomialRobustEstimator(0, evaluations,
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = new PROSACPolynomialRobustEstimator(2, wrongEvals,
                    this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testGetSetThreshold() throws LockedException {
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(estimator.getThreshold(),
                PROSACPolynomialRobustEstimator.DEFAULT_THRESHOLD, 0.0);

        // set new value
        estimator.setThreshold(1.0);

        // check correctness
        assertEquals(estimator.getThreshold(), 1.0, 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setThreshold(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetEvaluations() throws LockedException {
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

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
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertNull(estimator.getListener());

        // set new value
        estimator.setListener(this);

        // check correctness
        assertSame(estimator.getListener(), this);
    }

    @Test
    public void testGetSetProgressDelta() throws LockedException {
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(estimator.getProgressDelta(),
                PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA, 0.0);

        // set new value
        estimator.setProgressDelta(0.5f);

        // check correctness
        assertEquals(estimator.getProgressDelta(), 0.5, 0.0);

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
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(estimator.getConfidence(),
                PolynomialRobustEstimator.DEFAULT_CONFIDENCE, 0.0);

        // set new value
        estimator.setConfidence(0.5);

        // check correctness
        assertEquals(estimator.getConfidence(), 0.5, 0.0);

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
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(estimator.getMaxIterations(),
                PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS);

        // set new value
        estimator.setMaxIterations(10);

        // check correctness
        assertEquals(estimator.getMaxIterations(), 10);

        // Force IllegalArgumentException
        try {
            estimator.setMaxIterations(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsSetGeometricDistanceUsed() throws LockedException {
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(estimator.isGeometricDistanceUsed(),
                PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);

        // set new value
        estimator.setGeometricDistanceUsed(
                !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);

        //check correctness
        assertEquals(estimator.isGeometricDistanceUsed(),
                !PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE);
    }

    @Test
    public void testGetSetDegree() throws LockedException {
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

        // check default value
        assertEquals(estimator.getDegree(), PolynomialEstimator.MIN_DEGREE);

        // set new value
        estimator.setDegree(2);

        // check correctness
        assertEquals(estimator.getDegree(), 2);

        // Force IllegalArgumentException
        try {
            estimator.setDegree(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetQualityScores() throws LockedException {
        final PROSACPolynomialRobustEstimator estimator =
                new PROSACPolynomialRobustEstimator();

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
            final PROSACPolynomialRobustEstimator estimator =
                    new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(estimator.getDegree(), 1);
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

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateDirectAndDerivativeEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROSACPolynomialRobustEstimator estimator =
                    new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(estimator.getDegree(), 1);
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

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateIntegralEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROSACPolynomialRobustEstimator estimator =
                    new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(estimator.getDegree(), 1);
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

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateIntegralIntervalEvaluationsAlgebraicDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROSACPolynomialRobustEstimator estimator =
                    new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);

            // check default values
            assertEquals(estimator.getDegree(), 1);
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

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateDirectEvaluationsGeometricDistance()
            throws LockedException, NotReadyException, RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROSACPolynomialRobustEstimator estimator =
                    new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            // check default values
            assertEquals(estimator.getDegree(), 1);
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

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
            assertTrue(estimateNextIteration > 0);
            assertTrue(estimateProgressChange >= 0);
        }
    }

    @Test
    public void testEstimateDirectAndDerivativeEvaluationsGeometricDistance()
            throws LockedException, NotReadyException,
            RobustEstimatorException {

        for (int t = 0; t < TIMES; t++) {
            final PROSACPolynomialRobustEstimator estimator =
                    new PROSACPolynomialRobustEstimator();
            estimator.setListener(this);
            estimator.setGeometricDistanceUsed(true);

            // check default values
            assertEquals(estimator.getDegree(), 1);
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

            assertEquals(estimateStart, 0);
            assertEquals(estimateEnd, 0);
            assertEquals(estimateNextIteration, 0);
            assertEquals(estimateProgressChange, 0);
            assertTrue(estimator.isReady());
            assertFalse(estimator.isLocked());

            // estimate
            final Polynomial polynomial2 = estimator.estimate();

            // check correctness
            assertArrayEquals(polynomial2.getPolyParams(), polyParams,
                    ABSOLUTE_ERROR);
            assertEquals(estimateStart, 1);
            assertEquals(estimateEnd, 1);
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
