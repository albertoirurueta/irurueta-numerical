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

import com.irurueta.numerical.robust.RobustEstimatorMethod;
import org.junit.*;

import java.util.ArrayList;
import java.util.List;

import static com.irurueta.numerical.polynomials.estimators.PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA;
import static org.junit.Assert.*;

public class PolynomialRobustEstimatorTest implements
        PolynomialRobustEstimatorListener {

    @Test
    public void testCreate() {
        // test empty creator
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create();

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // test creator with degree
        estimator = PolynomialRobustEstimator.create(2);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with listener
        estimator = PolynomialRobustEstimator.create(this);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testCreteRANSAC() {
        // test creator with method
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with evaluations and method
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvals = new ArrayList<>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.RANSAC);

        // check
        assertTrue(estimator instanceof RANSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.RANSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testCreateLMedS() {
        // test creator with method
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with evaluations and method
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvals = new ArrayList<>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvals, this,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

        // test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.LMEDS);

        // check
        assertTrue(estimator instanceof LMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvals, this,
                    RobustEstimatorMethod.LMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testCreateMSAC() {
        // test creator with method
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());


        // test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations and method
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());


        // test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.MSAC);

        // check
        assertTrue(estimator instanceof MSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations, this,
                    RobustEstimatorMethod.MSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testCreatePROSAC() {
        // test empty creator
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());


        // test creator with degree
        estimator = PolynomialRobustEstimator.create(2,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with listener
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());


        // test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.PROSAC);

        // check
        assertTrue(estimator instanceof PROSACPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations, this,
                    RobustEstimatorMethod.PROSAC);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testCreatePROMedS() {
        // test empty creator
        PolynomialRobustEstimator estimator =
                PolynomialRobustEstimator.create(RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());


        // test creator with degree
        estimator = PolynomialRobustEstimator.create(2,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        final List<PolynomialEvaluation> wrongEvaluations = new ArrayList<>();
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with listener
        estimator = PolynomialRobustEstimator.create(this,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());


        // test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, this,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
                estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE,
                estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS,
                estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE,
                estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(wrongEvaluations, this,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertTrue(estimator instanceof PROMedSPolynomialRobustEstimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA,
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
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        estimator = null;
        try {
            estimator = PolynomialRobustEstimator.create(0, evaluations, this,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            estimator = PolynomialRobustEstimator.create(2, wrongEvaluations, this,
                    RobustEstimatorMethod.PROMEDS);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Override
    public void onEstimateStart(final PolynomialRobustEstimator estimator) {
    }

    @Override
    public void onEstimateEnd(final PolynomialRobustEstimator estimator) {
    }

    @Override
    public void onEstimateNextIteration(final PolynomialRobustEstimator estimator,
                                        final int iteration) {
    }

    @Override
    public void onEstimateProgressChange(final PolynomialRobustEstimator estimator,
                                         final float progress) {
    }
}
