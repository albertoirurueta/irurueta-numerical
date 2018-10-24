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
import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

@SuppressWarnings("Duplicates")
public class LMSEPolynomialEstimatorTest implements PolynomialEstimatorListener {
    
    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;
    
    private static final int MIN_DEGREE = 1;
    private static final int MAX_DEGREE = 5;
    
    private static final double ABSOLUTE_ERROR = 1e-8;
    
    private int estimateStart;
    private int estimateEnd;
    
    public LMSEPolynomialEstimatorTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }

    @Test
    public void testConstructor() {
        //empty constructor
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check correctness
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(), 
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //constructor with degree
        estimator = new LMSEPolynomialEstimator(2);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMSEPolynomialEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        
        
        //constructor with evaluations
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        estimator = new LMSEPolynomialEstimator(evaluations);
        
        //check correctness
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //constructor with listener
        estimator = new LMSEPolynomialEstimator(this);
        
        //check correctness
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //constructor with degree and evaluations
        estimator = new LMSEPolynomialEstimator(2, evaluations);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertNull(estimator.getListener());
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMSEPolynomialEstimator(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        
        
        //constructor with degree and listener
        estimator = new LMSEPolynomialEstimator(2, this);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        assertNull(estimator.getEvaluations());
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMSEPolynomialEstimator(0, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
        
        
        //constructor with evaluations and listener
        estimator = new LMSEPolynomialEstimator(evaluations, this);
        
        //check correctness
        assertEquals(estimator.getDegree(), 1);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        
        //constructor with degree, evaluations and listener
        estimator = new LMSEPolynomialEstimator(2, evaluations, this);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        assertFalse(estimator.isReady());
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        assertFalse(estimator.isLocked());
        assertSame(estimator.getListener(), this);
        assertFalse(estimator.isLMSESolutionAllowed());
        assertEquals(estimator.getType(),
                PolynomialEstimatorType.LMSE_POLYNOMIAL_ESTIMATOR);
        
        //Force IllegalArgumentException
        estimator = null;
        try {
            estimator = new LMSEPolynomialEstimator(0, evaluations, this);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(estimator);
    }
    
    @Test
    public void testIsSetLMSESolutionAllowed() throws LockedException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default value
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //set new value
        estimator.setLMSESolutionAllowed(true);
        
        //check correctness
        assertTrue(estimator.isLMSESolutionAllowed());
    }
    
    @Test
    public void testGetSetDegree() throws LockedException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default value
        assertEquals(estimator.getDegree(), 1);
        
        //set new value
        estimator.setDegree(2);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        
        //Force IllegalArgumentException
        try {
            estimator.setDegree(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testGetSetEvaluations() throws LockedException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default value
        assertNull(estimator.getEvaluations());
        
        //set new value
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        estimator.setEvaluations(evaluations);
        
        //check correctness
        assertSame(estimator.getEvaluations(), evaluations);
    }
    
    @Test
    public void testSetDegreeAndEvaluations() throws LockedException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertNull(estimator.getEvaluations());
        
        //set new values
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        estimator.setDegreeAndEvaluations(2, evaluations);
        
        //check correctness
        assertEquals(estimator.getDegree(), 2);
        assertSame(estimator.getEvaluations(), evaluations);
        
        //Force IllegalArgumentException
        try {
            estimator.setDegreeAndEvaluations(0, evaluations);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testIsReady() throws LockedException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default value
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = 
                new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
    }
    
    @Test
    public void testGetMinNumberOfEvaluations() {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int degree = randomizer.nextInt(MIN_DEGREE, MAX_DEGREE);
        assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(degree), 
                degree + 1);
        
        //Force IllegalArgumentException
        try {
            assertEquals(PolynomialEstimator.getMinNumberOfEvaluations(0), degree + 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator(degree);
        assertEquals(estimator.getMinNumberOfEvaluations(), degree + 1);
    }
    
    @Test
    public void testGetSetListener() throws LockedException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default value
        assertNull(estimator.getListener());
        
        //set new value
        estimator.setListener(this);
        
        //check correctness
        assertSame(estimator.getListener(), this);
    }
    
    @Test
    public void testEstimateWithDirectEvaluationsNoLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithDirectEvaluationsLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }
    
    @Test
    public void testEstimateWithDirectAndDerivativeEvaluationsNoLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations() - 1; i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        double value = polynomial.evaluateDerivative(x);
        DerivativePolynomialEvaluation eval = 
                new DerivativePolynomialEvaluation(x, value, 1);
        evaluations.add(eval);

        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithDirectAndDerivativeEvaluationLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2*estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        for(int i = 0; i < 2*estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluateDerivative(x);
            DerivativePolynomialEvaluation eval = 
                    new DerivativePolynomialEvaluation(x, value, 1);
            evaluations.add(eval);
        }

        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }
    
    @Test
    public void testEstimateWithIntegralEvaluationsNoLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double constant = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            Polynomial integral = polynomial.integrationAndReturnNew(constant);
            double value = integral.evaluate(x);
            
            IntegralPolynomialEvaluation eval = 
                    new IntegralPolynomialEvaluation(x, value, 
                            new double[]{constant}, 1);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithIntegralEvaluationsLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2*estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double constant = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            Polynomial integral = polynomial.integrationAndReturnNew(constant);
            double value = integral.evaluate(x);
            
            IntegralPolynomialEvaluation eval = 
                    new IntegralPolynomialEvaluation(x, value, 
                            new double[]{constant}, 1);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }
    
    @Test
    public void testEstimateWithIntegralIntervalEvaluationsNoLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);            
            double value = polynomial.integrateInterval(startX, endX);
            
            IntegralIntervalPolynomialEvaluation eval = 
                    new IntegralIntervalPolynomialEvaluation(startX, endX, 
                            value, 1);
            eval.setConstants(new double[]{0.0});
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithIntegralIntervalEvaluationsLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2*estimator.getMinNumberOfEvaluations(); i++) {
            double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);            
            double value = polynomial.integrateInterval(startX, endX);
            
            IntegralIntervalPolynomialEvaluation eval = 
                    new IntegralIntervalPolynomialEvaluation(startX, endX, 
                            value, 1);
            eval.setConstants(new double[]{0.0});
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }
    
    @Test
    public void testEstimateWithDirectEvaluationsNoLMSEAllowedSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator(2);
        
        //check default values
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        //create random 2nd degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);
        
        assertEquals(polynomial.getDegree(), 2);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithDirectEvaluationsLMSEAllowedSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator(2);
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        //create random 2nd degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);
        
        assertEquals(polynomial.getDegree(), 2);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);     
    }

    @Test
    public void testEstimateWithDirectAndSecondOrderDerivativeEvaluationsNoLMSEAllowedSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator(2);
        
        //check default values
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        //create random 2nd degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);
        
        assertEquals(polynomial.getDegree(), 2);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations() - 2; i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        double value = polynomial.evaluateDerivative(x);
        DerivativePolynomialEvaluation eval = 
                new DerivativePolynomialEvaluation(x, value, 1);
        evaluations.add(eval);
        
        x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        value = polynomial.evaluateSecondDerivative(x);
        eval = new DerivativePolynomialEvaluation(x, value, 2);
        evaluations.add(eval);
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);       
    }

    @Test
    public void testEstimateWithDirectAndSecondOrderDerivativeEvaluationLMSEAllowedSecondDegree()
            throws LockedException, NotReadyException,
            PolynomialEstimationException {
        
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator(2);
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 2);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        //create random 2nd degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[3];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);
        
        assertEquals(polynomial.getDegree(), 2);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 3);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double value = polynomial.evaluate(x);
            
            DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation(x, 
                    value);
            evaluations.add(eval);
        }
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        double value = polynomial.evaluateDerivative(x);
        DerivativePolynomialEvaluation eval = 
                new DerivativePolynomialEvaluation(x, value, 1);
        evaluations.add(eval);
        
        x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        value = polynomial.evaluateSecondDerivative(x);
        eval = new DerivativePolynomialEvaluation(x, value, 2);
        evaluations.add(eval);
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);       
    }
    
    @Test
    public void testEstimateWithSecondOrderIntegralEvaluationsNoLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            Polynomial integral = polynomial.nthIntegrationAndReturnNew(2, 
                    constants);
            double value = integral.evaluate(x);
            
            IntegralPolynomialEvaluation eval = 
                    new IntegralPolynomialEvaluation(x, value, constants, 2);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithSecondOrderIntegralEvaluationLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2 * estimator.getMinNumberOfEvaluations(); i++) {
            double x = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            Polynomial integral = polynomial.nthIntegrationAndReturnNew(2, 
                    constants);
            double value = integral.evaluate(x);
            
            IntegralPolynomialEvaluation eval = 
                    new IntegralPolynomialEvaluation(x, value, constants, 2);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }
    
    @Test
    public void testEstimateWithSecondOrderIntegralIntervalEvaluationsNoLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertFalse(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < estimator.getMinNumberOfEvaluations(); i++) {
            double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);  
            double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            
            double value = polynomial.nthOrderIntegrateInterval(startX, endX, 2,
                    constants);
            
            IntegralIntervalPolynomialEvaluation eval = 
                    new IntegralIntervalPolynomialEvaluation(startX, endX, 
                            value, 2);
            eval.setConstants(constants);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }

    @Test
    public void testEstimateWithSecondOrderIntegralIntervalEvaluationLMSEAllowed() 
            throws LockedException, NotReadyException, 
            PolynomialEstimationException {
        LMSEPolynomialEstimator estimator = new LMSEPolynomialEstimator();
        estimator.setLMSESolutionAllowed(true);
        
        //check default values
        assertEquals(estimator.getDegree(), 1);
        assertFalse(estimator.isReady());
        assertTrue(estimator.isLMSESolutionAllowed());
        
        //Force NotReadyException
        try {
            estimator.estimate();
            fail("NotReadyException expected but not thrown");
        } catch (NotReadyException ignore) { }
        
        
        //create random 1st degree polynomial
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double[] polyParams = new double[2];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        Polynomial polynomial = new Polynomial(polyParams);        
        
        assertEquals(polynomial.getDegree(), 1);
        
        assertEquals(estimator.getMinNumberOfEvaluations(), 2);
        List<PolynomialEvaluation> evaluations = new ArrayList<>();
        for(int i = 0; i < 2*estimator.getMinNumberOfEvaluations(); i++) {
            double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
            double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);  
            double[] constants = new double[2];
            randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            
            double value = polynomial.nthOrderIntegrateInterval(startX, endX, 2,
                    constants);
            
            IntegralIntervalPolynomialEvaluation eval = 
                    new IntegralIntervalPolynomialEvaluation(startX, endX, 
                            value, 2);
            eval.setConstants(constants);
            evaluations.add(eval);
        }
        
        estimator.setEvaluations(evaluations);
        
        assertTrue(estimator.isReady());
        assertSame(polynomial.getPolyParams(), polyParams);
        
        estimator.setListener(this);
        reset();
        
        assertEquals(estimateStart, 0);
        assertEquals(estimateEnd, 0);
        
        
        //estimate
        Polynomial polynomial2 = estimator.estimate();
        
        //check correctness
        assertArrayEquals(polynomial2.getPolyParams(), polyParams, 
                ABSOLUTE_ERROR);
        assertEquals(estimateStart, 1);
        assertEquals(estimateEnd, 1);
    }
    
    
    private void reset() {
        estimateStart = estimateEnd = 0;
    }

    @Override
    public void onEstimateStart(PolynomialEstimator estimator) {
        estimateStart++;
        checkIsLocked(estimator);
    }

    @Override
    public void onEstimateEnd(PolynomialEstimator estimator) {
        estimateEnd++;
        checkIsLocked(estimator);
    }
    
    private void checkIsLocked(PolynomialEstimator estimator) {
        assertTrue(estimator.isLocked());
        
        //Force LockedException
        try {
            estimator.setDegree(2);
            fail("LockedException expected but not thrown");
        } catch (LockedException ignore) { }
        try {
            estimator.setEvaluations(null);
            fail("LockedException expected but not thrown");
        } catch (LockedException ignore) { }
        try {
            estimator.setListener(null);
            fail("LockedException expected but not thrown");
        } catch (LockedException ignore) { }
        try {
            estimator.estimate();
            fail("LockedException expected but not thrown");
        } catch (LockedException ignore) {
        } catch (Exception ignore) {
            fail("LockedException expected but not thrown");
        }
        try {
            ((LMSEPolynomialEstimator)estimator).setLMSESolutionAllowed(true);
            fail("LockedException expected but not thrown");
        } catch (LockedException ignore) { }
    }
}
