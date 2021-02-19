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
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class BisectionSingleRootEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;

    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 1e-5;

    private double constant;
    private double root1;
    private double root2;
    private double root3;

    private final SingleDimensionFunctionEvaluatorListener constantPolynomial;
    private final SingleDimensionFunctionEvaluatorListener firstDegreePolynomial;
    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomial;
    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomialWithTwoComplexConjugateRoots;
    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomial;
    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithDoubleRoot;
    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithTripleRoot;
    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots;

    public BisectionSingleRootEstimatorTest() {

        constantPolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return constant;
            }
        };

        firstDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return point - root1;
            }
        };

        secondDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return (point - root1) * (point - root2);
            }
        };

        secondDegreePolynomialWithTwoComplexConjugateRoots = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return point * point + Math.abs(root1);
            }
        };

        thirdDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return (point - root1) * (point - root2) * (point - root3);
            }
        };

        thirdDegreePolynomialWithDoubleRoot = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return (point - root1) * (point - root1) * (point - root2);
            }
        };

        thirdDegreePolynomialWithTripleRoot = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return (point - root1) * (point - root1) * (point - root1);
            }
        };

        thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point - root1) * (point * point + Math.abs(root2));
                    }
                };
    }

    @Test
    public void testConstructor() throws NotAvailableException,
            InvalidBracketRangeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);

        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        BisectionSingleRootEstimator estimator;


        // Test 1st constructor
        estimator = new BisectionSingleRootEstimator();
        assertNotNull(estimator);

        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertEquals(estimator.getMaxEvaluationPoint(),
                BisectionSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 0.0);
        assertEquals(estimator.getMinEvaluationPoint(),
                BisectionSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertEquals(estimator.getTolerance(),
                BisectionSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());


        // Test 2nd constructor
        estimator = new BisectionSingleRootEstimator(constantPolynomial,
                minEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(estimator);

        assertEquals(estimator.getListener(), constantPolynomial);
        assertEquals(estimator.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        assertEquals(estimator.getMinEvaluationPoint(), minEvalPoint, 0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertEquals(estimator.getTolerance(), tolerance, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertFalse(estimator.isRootAvailable());

        // Force InvalidBracketRangeException
        estimator = null;
        try {
            estimator = new BisectionSingleRootEstimator(constantPolynomial,
                    maxEvalPoint, minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }

        // Force IllegalArgumentException
        try {
            estimator = new BisectionSingleRootEstimator(constantPolynomial,
                    minEvalPoint, maxEvalPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);
    }

    @Test
    public void testGetSetListenerAndAvailability() throws LockedException,
            NotAvailableException {

        final BisectionSingleRootEstimator estimator =
                new BisectionSingleRootEstimator();

        // check default values
        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());

        // set listener
        estimator.setListener(constantPolynomial);
        // check correctness
        assertEquals(estimator.getListener(), constantPolynomial);
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isReady());
    }

    @Test
    public void testGetSetBracketGetEvaluationPointsAndAvailability()
            throws NotAvailableException, LockedException,
            InvalidBracketRangeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);

        final BisectionSingleRootEstimator estimator =
                new BisectionSingleRootEstimator();

        // check default values
        assertTrue(estimator.isBracketAvailable());
        assertEquals(estimator.getMinEvaluationPoint(),
                BisectionSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(estimator.getMaxEvaluationPoint(),
                BisectionSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 0.0);

        // set new values
        estimator.setBracket(minEvalPoint, maxEvalPoint);
        // check correctness
        assertTrue(estimator.isBracketAvailable());
        assertEquals(estimator.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(estimator.getMaxEvaluationPoint(), maxEvalPoint, 0.0);

        // Force InvalidBracketRangeException
        try {
            estimator.setBracket(maxEvalPoint, minEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
    }

    @Test
    public void testGetSetTolerance() throws LockedException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        final BisectionSingleRootEstimator estimator =
                new BisectionSingleRootEstimator();

        // check default values
        assertEquals(estimator.getTolerance(),
                BisectionSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);

        // set new value
        estimator.setTolerance(tolerance);

        // check correctness
        assertEquals(estimator.getTolerance(), tolerance, 0.0);

        // Force IllegalArgumentException
        try {
            estimator.setTolerance(-tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testEstimate() throws LockedException, NotReadyException,
            InvalidBracketRangeException, RootEstimationException,
            NotAvailableException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        constant = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        root1 = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        root2 = randomizer.nextDouble(root1, MAX_EVAL_POINT);
        root3 = randomizer.nextDouble(root2, MAX_EVAL_POINT);

        // instantiate estimator with brackets for accuracy (otherwise estimation
        // might fail
        final BisectionSingleRootEstimator estimator =
                new BisectionSingleRootEstimator();

        // test constant polynomial (has no root)
        estimator.setListener(constantPolynomial);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (final RootEstimationException ignore) {
        }
        assertFalse(estimator.isLocked());
        try {
            estimator.estimate();
            fail("RootEstimationException expected but not thrown");
        } catch (final RootEstimationException ignore) {
        }
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isRootAvailable());
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }


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
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());


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
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());

        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root1 + root2), MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root2, estimator.getTolerance());


        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);


        // test 2nd degree polynomial with two complex conjugate roots
        estimator.setListener(
                secondDegreePolynomialWithTwoComplexConjugateRoots);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (final RootEstimationException ignore) {
        }
        assertFalse(estimator.isLocked());
        try {
            estimator.estimate();
            fail("RootEstimationException expected but not thrown");
        } catch (final RootEstimationException ignore) {
        }
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
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());

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
        assertEquals(estimator.getRoot(), root3, estimator.getTolerance());


        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);


        // test 3rd degree polynomial with double root
        // we need to properly set bracketing for each root and then refine the
        // result using estimate method
        estimator.setListener(thirdDegreePolynomialWithDoubleRoot);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root1 + root2), MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root2, estimator.getTolerance());


        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);


        // test 3rd degree polynomial with triple root
        estimator.setListener(thirdDegreePolynomialWithTripleRoot);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());


        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);


        // test 3rd degree polynomial with 1 real root and 2 conjugate complex
        // roots
        estimator.setListener(thirdDegreePolynomialWithTripleRoot);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT,
                0.5 * (MIN_EVAL_POINT + root2));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
    }
}
