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
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static com.irurueta.numerical.polynomials.estimators.PolynomialRobustEstimator.DEFAULT_PROGRESS_DELTA;
import static org.junit.jupiter.api.Assertions.*;

class PolynomialRobustEstimatorTest implements PolynomialRobustEstimatorListener {

    @Test
    void testCreate() {
        // test empty creator
        var estimator = PolynomialRobustEstimator.create();

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // test creator with degree
        estimator = PolynomialRobustEstimator.create(2);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0));

        // test creator with evaluations
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations));

        // test creator with listener
        estimator = PolynomialRobustEstimator.create(this);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2,
                wrongEvaluations));

        // test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, this));

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                this));

        // test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                this));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                this));
    }

    @Test
    void testCreteRANSAC() {
        // test creator with method
        var estimator = PolynomialRobustEstimator.create(RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0,
                RobustEstimatorMethod.RANSAC));

        // test creator with evaluations and method
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final var wrongEvals = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvals,
                RobustEstimatorMethod.RANSAC));

        // test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                RobustEstimatorMethod.RANSAC));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvals,
                RobustEstimatorMethod.RANSAC));

        // test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, this,
                RobustEstimatorMethod.RANSAC));

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvals, this,
                RobustEstimatorMethod.RANSAC));

        // test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this, RobustEstimatorMethod.RANSAC);

        // check
        assertInstanceOf(RANSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.RANSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                this, RobustEstimatorMethod.RANSAC));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvals,
                this, RobustEstimatorMethod.RANSAC));
    }

    @Test
    void testCreateLMedS() {
        // test creator with method
        var estimator = PolynomialRobustEstimator.create(RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0,
                RobustEstimatorMethod.LMEDS));

        // test creator with evaluations and method
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        final var wrongEvals = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvals,
                RobustEstimatorMethod.LMEDS));

        // test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                RobustEstimatorMethod.LMEDS));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvals,
                RobustEstimatorMethod.LMEDS));

        // test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, this,
                RobustEstimatorMethod.LMEDS));

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvals, this,
                RobustEstimatorMethod.LMEDS));

        // test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this, RobustEstimatorMethod.LMEDS);

        // check
        assertInstanceOf(LMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.LMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                this, RobustEstimatorMethod.LMEDS));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvals,
                this, RobustEstimatorMethod.LMEDS));
    }

    @Test
    void testCreateMSAC() {
        // test creator with method
        var estimator = PolynomialRobustEstimator.create(RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // test creator with degree and method
        estimator = PolynomialRobustEstimator.create(2, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0,
                RobustEstimatorMethod.MSAC));

        // test creator with evaluations and method
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                RobustEstimatorMethod.MSAC));

        // test creator with listener and method
        estimator = PolynomialRobustEstimator.create(this, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // test creator with degree, evaluations and method
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                RobustEstimatorMethod.MSAC));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                RobustEstimatorMethod.MSAC));

        // test creator with degree, listener and method
        estimator = PolynomialRobustEstimator.create(2, this, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, this,
                RobustEstimatorMethod.MSAC));

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                this, RobustEstimatorMethod.MSAC));

        // test creator with degree, evaluations, listener and method
        estimator = PolynomialRobustEstimator.create(2, evaluations, this, RobustEstimatorMethod.MSAC);

        // check
        assertInstanceOf(MSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertTrue(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.MSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                this, RobustEstimatorMethod.MSAC));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                this, RobustEstimatorMethod.MSAC));
    }

    @Test
    void testCreatePROSAC() {
        // test empty creator
        var estimator = PolynomialRobustEstimator.create(RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // test creator with degree
        estimator = PolynomialRobustEstimator.create(2, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0,
                RobustEstimatorMethod.PROSAC));

        // test creator with evaluations
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                RobustEstimatorMethod.PROSAC));

        // test creator with listener
        estimator = PolynomialRobustEstimator.create(this, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                RobustEstimatorMethod.PROSAC));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                RobustEstimatorMethod.PROSAC));

        // test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, this,
                RobustEstimatorMethod.PROSAC));

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                this, RobustEstimatorMethod.PROSAC));

        // test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this, RobustEstimatorMethod.PROSAC);

        // check
        assertInstanceOf(PROSACPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROSAC, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                this, RobustEstimatorMethod.PROSAC));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                this, RobustEstimatorMethod.PROSAC));
    }

    @Test
    void testCreatePROMedS() {
        // test empty creator
        var estimator = PolynomialRobustEstimator.create(RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // test creator with degree
        estimator = PolynomialRobustEstimator.create(2, RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0,
                RobustEstimatorMethod.PROMEDS));

        // test creator with evaluations
        final var evaluations = new ArrayList<PolynomialEvaluation>();
        evaluations.add(new DirectPolynomialEvaluation());
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(evaluations, RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertSame(estimator.getEvaluations(), evaluations);
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        final var wrongEvaluations = new ArrayList<PolynomialEvaluation>();
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                RobustEstimatorMethod.PROMEDS));

        // test creator with listener
        estimator = PolynomialRobustEstimator.create(this, RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(PolynomialEstimator.MIN_DEGREE, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // test creator with degree and evaluations
        evaluations.add(new DirectPolynomialEvaluation());
        estimator = PolynomialRobustEstimator.create(2, evaluations, RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                RobustEstimatorMethod.PROMEDS));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                RobustEstimatorMethod.PROMEDS));

        // test creator with degree and listener
        estimator = PolynomialRobustEstimator.create(2, this, RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertNull(estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, this,
                RobustEstimatorMethod.PROMEDS));

        // test creator with evaluations and listener
        estimator = PolynomialRobustEstimator.create(evaluations, this, RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(1, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(wrongEvaluations,
                this, RobustEstimatorMethod.PROMEDS));

        // test creator with degree, evaluations and listener
        estimator = PolynomialRobustEstimator.create(2, evaluations, this,
                RobustEstimatorMethod.PROMEDS);

        // check
        assertInstanceOf(PROMedSPolynomialRobustEstimator.class, estimator);
        assertSame(evaluations, estimator.getEvaluations());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertSame(this, estimator.getListener());
        assertFalse(estimator.isLocked());
        assertEquals(DEFAULT_PROGRESS_DELTA, estimator.getProgressDelta(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_CONFIDENCE, estimator.getConfidence(), 0.0);
        assertEquals(PolynomialRobustEstimator.DEFAULT_MAX_ITERATIONS, estimator.getMaxIterations());
        assertEquals(PolynomialRobustEstimator.DEFAULT_USE_GEOMETRIC_DISTANCE, estimator.isGeometricDistanceUsed());
        assertEquals(2, estimator.getDegree());
        assertFalse(estimator.isReady());
        assertNull(estimator.getQualityScores());
        assertEquals(RobustEstimatorMethod.PROMEDS, estimator.getMethod());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(0, evaluations,
                this, RobustEstimatorMethod.PROMEDS));
        assertThrows(IllegalArgumentException.class, () -> PolynomialRobustEstimator.create(2, wrongEvaluations,
                this, RobustEstimatorMethod.PROMEDS));
    }

    @Override
    public void onEstimateStart(final PolynomialRobustEstimator estimator) {
        // no action needed
    }

    @Override
    public void onEstimateEnd(final PolynomialRobustEstimator estimator) {
        // no action needed
    }

    @Override
    public void onEstimateNextIteration(final PolynomialRobustEstimator estimator, final int iteration) {
        // no action needed
    }

    @Override
    public void onEstimateProgressChange(final PolynomialRobustEstimator estimator, final float progress) {
        // no action needed
    }
}
