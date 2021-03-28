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

public class NewtonRaphsonSingleRootEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;

    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 1e-5;

    private double constant;
    private double root1;
    private double root2;
    private double root3;

    private final SingleDimensionFunctionEvaluatorListener constantPolynomial;
    private final SingleDimensionFunctionEvaluatorListener derivativeContantPolynomial;

    private final SingleDimensionFunctionEvaluatorListener firstDegreePolynomial;
    private final SingleDimensionFunctionEvaluatorListener derivativeFirstDegreePolynomial;

    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomial;
    private final SingleDimensionFunctionEvaluatorListener derivativeSecondDegreePolynomial;

    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomialWithDoubleRoot;
    private final SingleDimensionFunctionEvaluatorListener derivativeSecondDegreePolynomialWithDoubleRoot;

    private final SingleDimensionFunctionEvaluatorListener secondDegreePolynomialWithTwoComplexConjugateRoots;
    private final SingleDimensionFunctionEvaluatorListener derivativeSecondDegreePolynomialWithTwoComplexConjugateRoots;

    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomial;
    private final SingleDimensionFunctionEvaluatorListener derivativeThirdDegreePolynomial;

    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithDoubleRoot;
    private final SingleDimensionFunctionEvaluatorListener derivativeThirdDegreePolynomialWithDoubleRoot;

    private final SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots;
    private final SingleDimensionFunctionEvaluatorListener derivativeThirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots;


    public NewtonRaphsonSingleRootEstimatorTest() {

        constantPolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return constant;
            }
        };

        derivativeContantPolynomial =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return 0.0;
                    }
                };

        firstDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return (point - root1);
            }
        };

        derivativeFirstDegreePolynomial =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return 1.0;
                    }
                };

        secondDegreePolynomial =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point - root1) * (point - root2);
                    }
                };

        derivativeSecondDegreePolynomial =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return 2.0 * point - root1 - root2;
                    }
                };

        secondDegreePolynomialWithDoubleRoot =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point - root1) * (point - root1);
                    }
                };

        derivativeSecondDegreePolynomialWithDoubleRoot =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return 2.0 * (point - root1);
                    }
                };

        secondDegreePolynomialWithTwoComplexConjugateRoots =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point * point + Math.abs(root1));
                    }
                };

        derivativeSecondDegreePolynomialWithTwoComplexConjugateRoots =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return 2.0 * point;
                    }
                };

        thirdDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(final double point) {
                return (point - root1) * (point - root2) * (point - root3);
            }
        };

        derivativeThirdDegreePolynomial =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point - root2) * (point - root3) +
                                (point - root1) * (2.0 * point - root2 - root3);
                    }
                };

        thirdDegreePolynomialWithDoubleRoot =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point - root1) * (point - root1) * (point - root2);
                    }
                };

        derivativeThirdDegreePolynomialWithDoubleRoot =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return 2.0 * (point - root1) * (point - root2) +
                                (point - root1) * (point - root1);
                    }
                };

        thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point - root1) * (point * point + Math.abs(root2));
                    }
                };

        derivativeThirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots =
                new SingleDimensionFunctionEvaluatorListener() {

                    @Override
                    public double evaluate(final double point) {
                        return (point * point + Math.abs(root2)) +
                                2.0 * point * (point - root1);
                    }
                };
    }

    @Test
    public void testConstructor() throws NotAvailableException, InvalidBracketRangeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);
        final double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);

        NewtonRaphsonSingleRootEstimator estimator;

        // testing 1st constructor
        estimator = new NewtonRaphsonSingleRootEstimator();
        assertNotNull(estimator);

        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        try {
            estimator.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertEquals(estimator.getMaxEvaluationPoint(),
                NewtonRaphsonSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 0.0);
        assertEquals(estimator.getMinEvaluationPoint(),
                NewtonRaphsonSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertEquals(estimator.getTolerance(),
                NewtonRaphsonSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());


        // Test 2nd constructor
        estimator = new NewtonRaphsonSingleRootEstimator(constantPolynomial,
                minEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(estimator);

        assertEquals(estimator.getListener(), constantPolynomial);
        try {
            estimator.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
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
        assertFalse(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());

        // Force InvalidBracketRangeException
        estimator = null;
        try {
            estimator = new NewtonRaphsonSingleRootEstimator(constantPolynomial,
                    maxEvalPoint, minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
        // Force IllegalArgumentException
        try {
            estimator = new NewtonRaphsonSingleRootEstimator(constantPolynomial,
                    minEvalPoint, maxEvalPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);


        // test 3rd constructor
        estimator = new NewtonRaphsonSingleRootEstimator(constantPolynomial,
                derivativeContantPolynomial, minEvalPoint, maxEvalPoint,
                tolerance);
        assertNotNull(estimator);

        assertEquals(estimator.getListener(), constantPolynomial);
        assertEquals(estimator.getDerivativeListener(),
                derivativeContantPolynomial);
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
        assertTrue(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertFalse(estimator.isRootAvailable());

        // Force InvalidBracketRangeException
        estimator = null;
        try {
            estimator = new NewtonRaphsonSingleRootEstimator(constantPolynomial,
                    derivativeContantPolynomial, maxEvalPoint, minEvalPoint,
                    tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (final InvalidBracketRangeException ignore) {
        }
        // Force IllegalArgumentException
        try {
            estimator = new NewtonRaphsonSingleRootEstimator(constantPolynomial,
                    derivativeContantPolynomial, minEvalPoint, maxEvalPoint,
                    -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(estimator);

    }

    @Test
    public void testGetSetListenerAndDerivativeListenerAvailabilityAndIsReady()
            throws LockedException, NotAvailableException {

        final NewtonRaphsonSingleRootEstimator estimator =
                new NewtonRaphsonSingleRootEstimator();

        // check default values
        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.isListenerAvailable());
        try {
            estimator.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (final NotAvailableException ignore) {
        }
        assertFalse(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isReady());


        // set listener
        estimator.setListener(constantPolynomial);
        // check correctness
        assertEquals(estimator.getListener(), constantPolynomial);
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());

        // set derivative listener
        estimator.setDerivativeListener(derivativeContantPolynomial);
        // check correctness
        assertEquals(estimator.getDerivativeListener(),
                derivativeContantPolynomial);
        assertTrue(estimator.isDerivativeListenerAvailable());
        // because both delegate are available...
        assertTrue(estimator.isReady());
    }

    @Test
    public void testSetBracketGetEvaluationPointsAndAvailability()
            throws NotAvailableException, LockedException,
            InvalidBracketRangeException {

        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT,
                MAX_EVAL_POINT);
        final double maxEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);

        final NewtonRaphsonSingleRootEstimator estimator =
                new NewtonRaphsonSingleRootEstimator();

        // check default values
        assertTrue(estimator.isBracketAvailable());
        assertEquals(estimator.getMinEvaluationPoint(),
                NewtonRaphsonSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(estimator.getMaxEvaluationPoint(),
                NewtonRaphsonSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 0.0);

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

        final NewtonRaphsonSingleRootEstimator estimator =
                new NewtonRaphsonSingleRootEstimator();

        // check default values
        assertEquals(estimator.getTolerance(),
                NewtonRaphsonSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);

        // set new value
        estimator.setTolerance(tolerance);
        // Check correctness
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
        root1 = randomizer.nextDouble(MIN_EVAL_POINT, 0.2 * MAX_EVAL_POINT);
        root2 = randomizer.nextDouble(0.4 * MAX_EVAL_POINT, 0.6 * MAX_EVAL_POINT);
        root3 = randomizer.nextDouble(0.8 * MAX_EVAL_POINT, MAX_EVAL_POINT);

        // instantiate estimator with brackets for accuracy (otherwise estimation
        // might fail)
        final NewtonRaphsonSingleRootEstimator estimator =
                new NewtonRaphsonSingleRootEstimator();

        // test constant polynomial
        estimator.setListener(constantPolynomial);
        estimator.setDerivativeListener(derivativeContantPolynomial);
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
        estimator.setDerivativeListener(derivativeFirstDegreePolynomial);
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
        // we need to properly set bracketing for each root and then refine the
        // result using estimate method
        estimator.setListener(secondDegreePolynomial);
        estimator.setDerivativeListener(derivativeSecondDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.9 * root1, 1.1 * root1);
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


        // test 2n degree polynomial with double root
        estimator.setListener(secondDegreePolynomialWithDoubleRoot);
        estimator.setDerivativeListener(derivativeSecondDegreePolynomialWithDoubleRoot);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (final RootEstimationException ignore) {
        }
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);


        // test 2nd degree with two complex conjugate roots
        estimator.setListener(secondDegreePolynomialWithTwoComplexConjugateRoots);
        estimator.setDerivativeListener(derivativeSecondDegreePolynomialWithTwoComplexConjugateRoots);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (final RootEstimationException ignore) {
        }
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
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
        estimator.setDerivativeListener(derivativeThirdDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.9 * root1, 1.1 * root1);
        // reset bracket
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

        estimator.computeBracket(0.9 * root3, 1.1 * root3);
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
        estimator.setDerivativeListener(derivativeThirdDegreePolynomialWithDoubleRoot);
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
        estimator.setListener(thirdDegreePolynomial);
        estimator.setDerivativeListener(derivativeThirdDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.9 * root1, 1.1 * root1);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());

        // reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);


        // test third degree polynomial with 1 real root and 2 conjugate complex
        // roots
        estimator.setListener(
                thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots);
        estimator.setDerivativeListener(
                derivativeThirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots);
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
