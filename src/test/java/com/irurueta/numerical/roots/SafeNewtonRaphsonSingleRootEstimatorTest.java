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

import com.irurueta.numerical.*;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.*;

public class SafeNewtonRaphsonSingleRootEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;
    
    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 1e-5;
    
    private double constant;
    private double root1;
    private double root2;
    private double root3;
    
    private SingleDimensionFunctionEvaluatorListener constantPolynomial;
    private SingleDimensionFunctionEvaluatorListener derivativeContantPolynomial;
    
    private SingleDimensionFunctionEvaluatorListener firstDegreePolynomial;
    private SingleDimensionFunctionEvaluatorListener derivativeFirstDegreePolynomial;
    
    private SingleDimensionFunctionEvaluatorListener secondDegreePolynomial;
    private SingleDimensionFunctionEvaluatorListener derivativeSecondDegreePolynomial;

    private SingleDimensionFunctionEvaluatorListener secondDegreePolynomialWithTwoComplexConjugateRoots;
    private SingleDimensionFunctionEvaluatorListener derivativeSecondDegreePolynomialWithTwoComplexConjugateRoots;
    
    private SingleDimensionFunctionEvaluatorListener thirdDegreePolynomial;
    private SingleDimensionFunctionEvaluatorListener derivativeThirdDegreePolynomial;
    
    private SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithDoubleRoot;
    private SingleDimensionFunctionEvaluatorListener derivativeThirdDegreePolynomialWithDoubleRoot;

    private SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots;
    private SingleDimensionFunctionEvaluatorListener derivativeThirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots;
    
    
    public SafeNewtonRaphsonSingleRootEstimatorTest() {
        constantPolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return constant;
            }
        };
        
        derivativeContantPolynomial = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return 0.0;
            }
        };
        
        firstDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1);
            }
        };
        
        derivativeFirstDegreePolynomial = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return 1.0;
            }
        };
        
        secondDegreePolynomial = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point - root2);
            }
        };
        
        derivativeSecondDegreePolynomial = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return 2.0 * point - root1 - root2;
            }
        };
        
        secondDegreePolynomialWithTwoComplexConjugateRoots = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point * point + Math.abs(root1));
            }
        };
        
        derivativeSecondDegreePolynomialWithTwoComplexConjugateRoots = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return 2.0 * point;
            }
        };
        
        thirdDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point - root2) * (point - root3);
            }
        };
        
        derivativeThirdDegreePolynomial = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root2) * (point - root3) +
                        (point - root1) * (2.0 * point - root2 - root3);
            }
        };
        
        thirdDegreePolynomialWithDoubleRoot = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point - root1) * (point - root2);
            }
        };
        
        derivativeThirdDegreePolynomialWithDoubleRoot = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return 2.0 * (point - root1) * (point - root2) + 
                        (point - root1) * (point - root1);
            }
        };

        thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point * point + Math.abs(root2));
            }
        };
        
        derivativeThirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots = 
                new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point * point + Math.abs(root2)) + 
                        2.0 * point * (point - root1);
            }
        };        
    }

    @BeforeClass
    public static void setUpClass() { }

    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() throws NotAvailableException,
            InvalidBracketRangeException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(minEvalPoint, 
                MAX_EVAL_POINT);
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        SafeNewtonRaphsonSingleRootEstimator estimator;
        
        //testing 1st constructor
        estimator = new SafeNewtonRaphsonSingleRootEstimator();
        assertNotNull(estimator);
        
        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            estimator.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getMaxEvaluationPoint(),
                SafeNewtonRaphsonSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 
                0.0);
        assertEquals(estimator.getMinEvaluationPoint(),
                SafeNewtonRaphsonSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 
                0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getTolerance(),
                SafeNewtonRaphsonSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());
        
        
        //Test 2nd constructor
        estimator = new SafeNewtonRaphsonSingleRootEstimator(constantPolynomial,
                minEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(estimator);
        
        assertEquals(estimator.getListener(), constantPolynomial);
        try {
            estimator.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        assertEquals(estimator.getMinEvaluationPoint(), minEvalPoint, 0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getTolerance(), tolerance, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());
        
        //Force InvalidBracketRangeException
        estimator = null;
        try {
            estimator = new SafeNewtonRaphsonSingleRootEstimator(constantPolynomial,
                    maxEvalPoint, minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException ignore) { }
        //Force IllegalArgumentException
        try {
            estimator = new SafeNewtonRaphsonSingleRootEstimator(constantPolynomial,
                    minEvalPoint, maxEvalPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        
        
        //test 3rd constructor
        estimator = new SafeNewtonRaphsonSingleRootEstimator(constantPolynomial,
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
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getTolerance(), tolerance, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertTrue(estimator.isListenerAvailable());
        assertTrue(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertFalse(estimator.isRootAvailable());
        
        //Force InvalidBracketRangeException
        estimator = null;
        try {
            estimator = new SafeNewtonRaphsonSingleRootEstimator(constantPolynomial,
                    derivativeContantPolynomial, maxEvalPoint, minEvalPoint, 
                    tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException ignore) { }
        //Force IllegalArgumentException
        try {
            estimator = new SafeNewtonRaphsonSingleRootEstimator(constantPolynomial,
                    derivativeContantPolynomial, minEvalPoint, maxEvalPoint, 
                    -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);

    }
    
    @Test
    public void testGetSetListenerAndDerivativeListenerAvailabilityAndIsReady() 
            throws LockedException, NotAvailableException {
        
        SafeNewtonRaphsonSingleRootEstimator estimator =
                new SafeNewtonRaphsonSingleRootEstimator();
        
        //check default values
        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isListenerAvailable());
        try {
            estimator.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isDerivativeListenerAvailable());
        assertFalse(estimator.isReady());
        
        
        //set listener
        estimator.setListener(constantPolynomial);
        //check correctness
        assertEquals(estimator.getListener(), constantPolynomial);
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());
        
        //set derivative listener
        estimator.setDerivativeListener(derivativeContantPolynomial);
        //check correctness
        assertEquals(estimator.getDerivativeListener(), 
                derivativeContantPolynomial);
        assertTrue(estimator.isDerivativeListenerAvailable());
        //because both delegate are available...
        assertTrue(estimator.isReady());
    }
    
    @Test
    public void testSetBracketGetEvaluationPointsAndAvailability() 
            throws NotAvailableException, LockedException, 
            InvalidBracketRangeException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(minEvalPoint, 
                MAX_EVAL_POINT);
        
        SafeNewtonRaphsonSingleRootEstimator estimator = 
                new SafeNewtonRaphsonSingleRootEstimator();
        
        //check default values
        assertTrue(estimator.isBracketAvailable());
        assertEquals(estimator.getMinEvaluationPoint(),
                SafeNewtonRaphsonSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 
                0.0);
        assertEquals(estimator.getMaxEvaluationPoint(),
                SafeNewtonRaphsonSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 
                0.0);
        
        //set new values
        estimator.setBracket(minEvalPoint, maxEvalPoint);
        //check correctness
        assertTrue(estimator.isBracketAvailable());
        assertEquals(estimator.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(estimator.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        
        //Force InvalidBracketRangeException
        try {
            estimator.setBracket(maxEvalPoint, minEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException ignore) { }
    }
    
    @Test
    public void testGetSetTolerance() throws LockedException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        SafeNewtonRaphsonSingleRootEstimator estimator = 
                new SafeNewtonRaphsonSingleRootEstimator();
        
        //check default values
        assertEquals(estimator.getTolerance(),
                SafeNewtonRaphsonSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);
        
        //set new value
        estimator.setTolerance(tolerance);
        //Check correctness
        assertEquals(estimator.getTolerance(), tolerance, 0.0);
        
        //Force IllegalArgumentException
        try {
            estimator.setTolerance(-tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException, 
        InvalidBracketRangeException, RootEstimationException, NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        constant = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        root1 = randomizer.nextDouble(MIN_EVAL_POINT, 0.2 * MAX_EVAL_POINT);
        root2 = randomizer.nextDouble(0.4 * MAX_EVAL_POINT, 0.6 * MAX_EVAL_POINT);
        root3 = randomizer.nextDouble(0.8 * MAX_EVAL_POINT, MAX_EVAL_POINT);
        
        //instantiate estimator with brackets for accuracy (otherwise estimation
        //might fail)
        SafeNewtonRaphsonSingleRootEstimator estimator = 
                new SafeNewtonRaphsonSingleRootEstimator();
        
        //test constant polynomial
        estimator.setListener(constantPolynomial);
        estimator.setDerivativeListener(derivativeContantPolynomial);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (RootEstimationException ignore) { }
        assertFalse(estimator.isLocked());
        try {
            estimator.estimate();
            fail("RootEstimationException expected but not thrown");
        } catch (RootEstimationException ignore) { }
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isRootAvailable());
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //test 1st degree polynomial
        estimator.setListener(firstDegreePolynomial);
        estimator.setDerivativeListener(derivativeFirstDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //test 2nd degree polynomial
        //we need to properly set bracketing for each root and then refine the
        //result using estimate method
        estimator.setListener(secondDegreePolynomial);
        estimator.setDerivativeListener(derivativeSecondDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, 0.5 * (root1 + root2));
        estimator.setBracket(MIN_EVAL_POINT, 0.5 * (root1 + root2));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
        
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root1 + root2), MAX_EVAL_POINT);
        estimator.setBracket(0.5 * (root1 + root2), MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root2, estimator.getTolerance());
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        /*
        //test 2n degree polynomial with double root
        estimator.setListener(secondDegreePolynomialWithDoubleRoot);
        estimator.setDerivativeListener(derivativeSecondDegreePolynomialWithDoubleRoot);
        assertFalse(estimator.isLocked());
        try{
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        }catch(RootEstimationException e){}
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        */
        
        
        
        //test 2nd degree with two complex conjugate roots
        estimator.setListener(secondDegreePolynomialWithTwoComplexConjugateRoots);
        estimator.setDerivativeListener(derivativeSecondDegreePolynomialWithTwoComplexConjugateRoots);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (RootEstimationException ignore) { }
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        try {
            estimator.estimate();
            fail("RootEstimationException expected but not thrown");
        } catch (RootEstimationException ignore) { }
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isRootAvailable());
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        
        
        //test 3rd degree polynomial
        //we need to properly set bracketing for each root and then refine the
        //result using estimate method
        estimator.setListener(thirdDegreePolynomial);
        estimator.setDerivativeListener(derivativeThirdDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, 0.5 * (root1 + root2));
        //reset bracket
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
        assertFalse(estimator.isLocked());
        
        estimator.computeBracket(0.5 * (root1 + root2), 0.5 * (root2 + root3));
        estimator.setBracket(0.5 * (root1 + root2), 0.5 * (root2 + root3));
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
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        
        
        //test 3rd degree polynomial with double root
        //we need to properly set bracketing for each root and then refine the
        //result using estimate method
        estimator.setListener(thirdDegreePolynomialWithDoubleRoot);
        estimator.setDerivativeListener(derivativeThirdDegreePolynomialWithDoubleRoot);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(0.5 * (root1 + root2), MAX_EVAL_POINT);
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root2, estimator.getTolerance());
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        
        
        //test 3rd degree polynomial with triple root
        estimator.setListener(thirdDegreePolynomial);
        estimator.setDerivativeListener(derivativeThirdDegreePolynomial);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, 0.5 * (root1 + root2));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        
        
        //test third degree polynomial with 1 real root and 2 conjugate complex
        //roots
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
