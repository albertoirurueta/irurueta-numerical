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

public class SecantSingleRootEstimatorTest {
    
    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;
    
    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 1e-5;

    private static final int TIMES = 10;
    
    private double constant;
    private double root1;
    private double root2;
    private double root3;
    
    private SingleDimensionFunctionEvaluatorListener constantPolynomial;
    private SingleDimensionFunctionEvaluatorListener firstDegreePolynomial;
    private SingleDimensionFunctionEvaluatorListener secondDegreePolynomial;
    private SingleDimensionFunctionEvaluatorListener secondDegreePolynomialWithTwoComplexConjugateRoots;
    private SingleDimensionFunctionEvaluatorListener thirdDegreePolynomial;
    private SingleDimensionFunctionEvaluatorListener thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots;    
    
    public SecantSingleRootEstimatorTest() {
        
        constantPolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return constant;
            }
        };
        
        firstDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return point - root1;
            }
        };
        
        secondDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point - root2);
            }
        };
        
        secondDegreePolynomialWithTwoComplexConjugateRoots = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return point * point + Math.abs(root1);
            }
        };
        
        thirdDegreePolynomial = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point - root2) * (point - root3);
            }
        };
        
        thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) {
                return (point - root1) * (point * point + Math.abs(root2));
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
    public void testConstructor() throws NotAvailableException, InvalidBracketRangeException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(minEvalPoint, 
                MAX_EVAL_POINT);
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        SecantSingleRootEstimator estimator;
        
        
        //testing 1st constructor
        estimator = new SecantSingleRootEstimator();
        assertNotNull(estimator);
        
        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getMaxEvaluationPoint(),
                SecantSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 0.0);
        assertEquals(estimator.getMinEvaluationPoint(),
                SecantSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getTolerance(), 
                SecantSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertFalse(estimator.isRootAvailable());
        
        
        
        //testing 2nd constructor
        estimator = new SecantSingleRootEstimator(constantPolynomial, 
                minEvalPoint, maxEvalPoint, tolerance);
        assertNotNull(estimator);
        
        assertEquals(estimator.getListener(), constantPolynomial);
        assertEquals(estimator.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        assertEquals(estimator.getMinEvaluationPoint(), minEvalPoint, 0.0);
        try {
            estimator.getRoot();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(estimator.getTolerance(), tolerance, 0.0);
        assertTrue(estimator.isBracketAvailable());
        assertTrue(estimator.isListenerAvailable());
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertFalse(estimator.isRootAvailable());
        
        //Force InvalidBracketRangeException
        estimator = null;
        try {
            estimator = new SecantSingleRootEstimator(constantPolynomial,
                    maxEvalPoint, minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException ignore) { }
        
        //Force IllegalArgumentException
        try {
            estimator = new SecantSingleRootEstimator(constantPolynomial,
                    minEvalPoint, maxEvalPoint, -tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
    }
    
    @Test
    public void testGetSetListenerAvailabilityAndIsReady() 
            throws LockedException, NotAvailableException {
        
        SecantSingleRootEstimator estimator = new SecantSingleRootEstimator();
        
        //check default values
        try {
            estimator.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(estimator.isListenerAvailable());
        assertFalse(estimator.isReady());
        
        //set listener
        estimator.setListener(constantPolynomial);
        //check correctness
        assertEquals(estimator.getListener(), constantPolynomial);
        assertTrue(estimator.isListenerAvailable());
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
        
        SecantSingleRootEstimator estimator = new SecantSingleRootEstimator();
        
        //check default values
        assertTrue(estimator.isBracketAvailable());
        assertEquals(estimator.getMinEvaluationPoint(),
                SecantSingleRootEstimator.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(estimator.getMaxEvaluationPoint(),
                SecantSingleRootEstimator.DEFAULT_MAX_EVAL_POINT, 0.0);
        
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
        
        SecantSingleRootEstimator estimator = new SecantSingleRootEstimator();
        
        //check default values
        assertEquals(estimator.getTolerance(),
                SecantSingleRootEstimator.DEFAULT_TOLERANCE, 0.0);
        
        //set new value
        estimator.setTolerance(tolerance);
        //check correctness
        assertEquals(estimator.getTolerance(), tolerance, 0.0);
        
        //Force IllegalArgumentException
        try {
            estimator.setTolerance(-tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testEstimate() throws LockedException, NotReadyException, 
            InvalidBracketRangeException, RootEstimationException,
            NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        constant = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        root1 = randomizer.nextDouble(MIN_EVAL_POINT, 0.2 * MAX_EVAL_POINT);
        root2 = randomizer.nextDouble(0.4 * MAX_EVAL_POINT, 0.6 * MAX_EVAL_POINT);
        root3 = randomizer.nextDouble(0.8 * MAX_EVAL_POINT, MAX_EVAL_POINT);

        //instantiate estimator with brackets for accuracy (otherwise estimation
        //might fail
        SecantSingleRootEstimator estimator = 
                new SecantSingleRootEstimator();
        
        //test constant polynomial (has no root)
        estimator.setListener(constantPolynomial);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (RootEstimationException ignore){}
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

        
        
        
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        
        //test 2nd degree polynomial with two complex conjugate roots
        estimator.setListener(
                secondDegreePolynomialWithTwoComplexConjugateRoots);
        assertFalse(estimator.isLocked());
        try {
            estimator.computeBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
            fail("RootEstimationException expected but not thrown");
        } catch (RootEstimationException ignore){}
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
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
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
            if (Math.abs(estimator.getRoot() - root3) > estimator.getTolerance()) {
                continue;
            }
            assertEquals(estimator.getRoot(), root3, estimator.getTolerance());

            numValid++;
            break;
        }

        assertTrue(numValid > 0);

        
                
        //reset bracket
        estimator.setBracket(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        
        //test 3rd degree polynomial with 1 real root and 2 conjugate complex 
        //roots
        estimator.setListener(thirdDegreePolynomialWithOneRealRootAndTwoComplexConjugateRoots);
        assertFalse(estimator.isLocked());
        estimator.computeBracket(MIN_EVAL_POINT, 
                0.5 * (root1 + root2));
        assertFalse(estimator.isLocked());
        estimator.estimate();
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isRootAvailable());
        assertEquals(estimator.getRoot(), root1, estimator.getTolerance());
      
    }    
}
