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

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class PolynomialEstimatorTest implements PolynomialEstimatorListener {

    @Test
    public void testCreate() {
        // default type
        PolynomialEstimator estimator = PolynomialEstimator.create();

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree
        estimator = PolynomialEstimator.create(2);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        estimator = PolynomialEstimator.create(evaluations);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with listener
        estimator = PolynomialEstimator.create(this);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree and evaluations
        estimator = PolynomialEstimator.create(2, evaluations);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree and listener
        estimator = PolynomialEstimator.create(2, this);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with evaluations and listener
        estimator = PolynomialEstimator.create(evaluations, this);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree, evaluations and listener
        estimator = PolynomialEstimator.create(2, evaluations, this);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());
    }

    @Test
    public void testCreateLMSE() {
        // default type
        PolynomialEstimator estimator = PolynomialEstimator.create(
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree
        estimator = PolynomialEstimator.create(2,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        estimator = PolynomialEstimator.create(evaluations,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with listener
        estimator = PolynomialEstimator.create(this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree and evaluations
        estimator = PolynomialEstimator.create(2, evaluations,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree and listener
        estimator = PolynomialEstimator.create(2, this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with evaluations and listener
        estimator = PolynomialEstimator.create(evaluations, this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree, evaluations and listener
        estimator = PolynomialEstimator.create(2, evaluations, this,
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof LMSEPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR,
                estimator.getType());
    }

    @Test
    public void testCreateWeighted() {
        // default type
        PolynomialEstimator estimator = PolynomialEstimator.create(
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree
        estimator = PolynomialEstimator.create(2,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with evaluations
        final List<PolynomialEvaluation> evaluations = new ArrayList<>();
        estimator = PolynomialEstimator.create(evaluations,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with listener
        estimator = PolynomialEstimator.create(this,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree and evaluations
        estimator = PolynomialEstimator.create(2, evaluations,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree and listener
        estimator = PolynomialEstimator.create(2, this,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with evaluations and listener
        estimator = PolynomialEstimator.create(evaluations, this,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(1, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(2, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());

        // default type with degree, evaluations and listener
        estimator = PolynomialEstimator.create(2, evaluations, this,
                PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR);

        // check correctness
        assertTrue(estimator instanceof WeightedPolynomialEstimator);
        assertEquals(2, estimator.getDegree());
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(3, estimator.getMinNumberOfEvaluations());
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertEquals(PolynomialEstimatorType.WEIGHTED_POLYNOMIAL_ESTIMATOR,
                estimator.getType());
    }

    @Override
    public void onEstimateStart(final PolynomialEstimator estimator) {
    }

    @Override
    public void onEstimateEnd(final PolynomialEstimator estimator) {
    }
}
