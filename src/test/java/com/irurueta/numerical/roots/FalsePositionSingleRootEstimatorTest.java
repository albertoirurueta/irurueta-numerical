/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.roots;

import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class FalsePositionSingleRootEstimatorTest {

    public static final double MIN_EVAL_POINT = 0.0;
    public static final double MAX_EVAL_POINT = 1.0;

    public static final double MIN_TOLERANCE = 3e-8;
    public static final double MAX_TOLERANCE = 1e-5;

    private double constant;
    private double root1;
    private double root2;
    private double root3;

    private final SingleDimensionFunctionEvaluatorListener constantPolynomial = point -> constant;
    private final SingleDimensionFunctionEvaluatorListener firstDegreePolynomial = point -> (point - root1);
    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomial =
            point -> (point - root1) * (point - root2);
    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomialWithTwoComplexConjugateRoots =
            point -> (point * point + Math.abs(root1));
    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomial =
            point -> (point - root1) * (point - root2) * (point - root3);
    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithTripleRoot =
            point -> (point - root1) * (point - root1) * (point - root1);

    @Test
    void testConstructor() throws NotAvailableException, InvalidBracketRangeException {

        final var randomizer = new UniformRandomizer();
        final var minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var maxEvalPoint = randomizer.nextDouble(minEvalPoint, MAX_EVAL_POINT);
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        // Testing 1st constructor
        var estimator = new FalsePositionSingleRootEstimator();
        assertNotNull(estimator);

        assertThrows(NotAvailableException.class, estimator::getListener);
        assertEquals(FalsePositionSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, estimator.getMaxEvaluationPoint(),
                0.0);
        assertEquals(FalsePositionSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, estimator.getMinEvaluationPoint(),
                0.0);
        assertThrows(NotAvailableException.class, estimator::getRoot);
        assertEquals(FalsePositionSingleRootEstimator.DEFAULT_TOLERANCE, estimator.getTolerance(), 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());

        // Testing 2nd constructor
        estimator = new FalsePositionSingleRootEstimator(constantPolynomial, minEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(estimator);

        assertEquals(constantPolynomial, estimator.getListener());
        assertEquals(maxEvalPoint, estimator.getMaxEvaluationPoint(), 0.0);
        assertEquals(minEvalPoint, estimator.getMinEvaluationPoint(), 0.0);
        assertThrows(NotAvailableException.class, estimator::getRoot);
        assertEquals(tolerance, estimator.getTolerance(), 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertFalse(estimator.isRootAvailable());

        // Force InvalidBracketRangeException
        assertThrows(InvalidBracketRangeException.class, () -> new FalsePositionSingleRootEstimator(constantPolynomial,
                maxEvalPoint, minEvalPoint, tolerance));

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new FalsePositionSingleRootEstimator(constantPolynomial,
                minEvalPoint, maxEvalPoint, -tolerance));
    }

    @Test
    void testGetSetListenerAvailabilityAndIsReady() throws LockedException, NotAvailableException {

        final var estimator = new FalsePositionSingleRootEstimator();

        // check default values
        assertThrows(NotAvailableException.class, estimator::getListener);
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());

        // set listener
        estimator.setListener(constantPolynomial);
        // check correctness
        assertEquals(constantPolynomial, estimator.getListener());
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isReady());
    }

    @Test
    void testSetBracketGetEvaluationPointsAndAvailability() throws NotAvailableException, LockedException,
            InvalidBracketRangeException {

        final var randomizer = new UniformRandomizer();
        final var minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        final var maxEvalPoint = randomizer.nextDouble(minEvalPoint, MAX_EVAL_POINT);

        final var estimator = new FalsePositionSingleRootEstimator();

        // check default values
        assertTrue(estimator.isBracketAvailable());
        assertEquals(FalsePositionSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, estimator.getMinEvaluationPoint(),
                0.0);
        assertEquals(FalsePositionSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, estimator.getMaxEvaluationPoint(),
                0.0);

        // set new values
        estimator.setBracket(minEvalPoint, maxEvalPoint);
        // check correctness
        assertTrue(estimator.isBracketAvailable());
        assertEquals(minEvalPoint, estimator.getMinEvaluationPoint(), 0.0);
        assertEquals(maxEvalPoint, estimator.getMaxEvaluationPoint(), 0.0);

        // Force InvalidBracketRangeException
        assertThrows(InvalidBracketRangeException.class, () -> estimator.setBracket(maxEvalPoint, minEvalPoint));
    }

    @Test
    void testGetSetTolerance() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final var estimator = new FalsePositionSingleRootEstimator();

        // check default values
        assertEquals(FalsePositionSingleRootEstimator.DEFAULT_TOLERANCE, estimator.getTolerance(), 0.0);

        // set new value
        estimator.setTolerance(tolerance);
        // check correctness
        assertEquals(estimator.getTolerance(), tolerance, 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setTolerance(-tolerance));
    }

    @Test
    void testEstimate() throws LockedException, NotReadyException, InvalidBracketRangeException,
            RootEstimationException, NotAvailableException {
        final var randomizer = new UniformRandomizer();
        constant = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        root1 = randomizer.nextDouble(MIN_EVAL_POINT, 0.2 * MAX_EVAL_POINT);
        root2 = randomizer.nextDouble(0.4 * MAX_EVAL_POINT, 0.6 * MAX_EVAL_POINT);
        root3 = randomizer.nextDouble(0.8 * MAX_EVAL_POINT, MAX_EVAL_POINT);

        // instantiate estimator with brackets for accuracy (otherwise estimation
        // might fail
        final var estimator = new FalsePositionSingleRootEstimator();

        // test constant polynomial (has no root)
        estimator.setListener(constantPolynomial);
        assertFalse(estimator.isLocked());
        assertThrows(RootEstimationException.class, () -> estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT));
        assertFalse(estimator.isLocked());
        assertThrows(RootEstimationException.class, estimator::estimate);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isRootAvailable());
        assertThrows(NotAvailableException.class, estimator::getRoot);

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // test 1st degree polynomial
        estimator.setListener(firstDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, 2.0 * estimator.getTolerance());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // test 2nd degree polynomial
        estimator.setListener(secondDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, 0.5 * (root1 + root2));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, 2.0 * estimator.getTolerance());

        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root1 + root2), 0.5 * (root2 + root3));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root2, 2.0 * estimator.getTolerance());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // test 2nd degree polynomial with two complex conjugate roots
        estimator.setListener(secondDegreePolynomialWithTwoComplexConjugateRoots);
        assertFalse(estimator.isLocked());
        assertThrows(RootEstimationException.class, () -> estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT));
        assertFalse(estimator.isLocked());
        assertThrows(RootEstimationException.class, estimator::estimate);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isRootAvailable());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // test 3rd degree polynomial
        // we need to properly set bracketing for each root and then refine the
        // result using estimate method
        estimator.setListener(thirdDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, 0.5 * (root1 + root2));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, 2.0 * estimator.getTolerance());

        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root1 + root2), 0.5 * (root2 + root3));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root2, estimator.getTolerance());

        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root2 + root3), MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root3, 2.0 * estimator.getTolerance());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // test 3rd degree polynomial with triple root
        estimator.setListener(thirdDegreePolynomialWithTripleRoot);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.9 * root1, 1.1 * root1);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, 2.0 * estimator.getTolerance());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);

        // test 3rd degree polynomial with 1 real root and 2 conjugate complex
        // roots
        estimator.setListener(thirdDegreePolynomialWithTripleRoot);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.9 * root1, 1.1 * root1);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, 2.0 * estimator.getTolerance());
    }
}
