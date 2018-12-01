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
package com.irurueta.numerical.optimization;

import com.irurueta.numerical.*;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.*;

public class QuasiNewtonMultiOptimizerTest {
    
    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 4;
    
    private static final double MIN_EVAL_POINT = -1e3;
    private static final double MAX_EVAL_POINT = 1e3;
    
    private static final double MIN_OFFSET = -1e3;
    private static final double MAX_OFFSET = 1e3;
    
    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;
    
    private static final double MIN_TOLERANCE = 3e-8;
    private static final double MAX_TOLERANCE = 3e-6;
    
    private static final double ABSOLUTE_ERROR = 1e-6;
    
    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;
    
    private MultiDimensionFunctionEvaluatorListener listener;
    private GradientFunctionEvaluatorListener gradientListener;
    
    public QuasiNewtonMultiOptimizerTest() {
        listener = new MultiDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double[] point) {
                int dims = Math.min(Math.min(point.length, minimum.length), 
                        width.length);
                
                double value = 1.0;
                for(int i = 0; i < dims; i++){
                    value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
                }
                
                value += offset;
                return value;
            }
        };
        
        gradientListener = new GradientFunctionEvaluatorListener() {

            @Override
            public void evaluateGradient(double[] params, double[] result) {
                int dims = Math.min(Math.min(params.length, minimum.length), 
                        width.length);
                
                double value;
                for(int j = 0; j < dims; j++){
                    value = 1.0;
                    for(int i = 0; i < dims; i++){
                        if(i != j){
                            value *= Math.pow(params[i] - minimum[i], 2.0) / 
                                    width[i];
                        }else{
                            value *= 2.0 * (params[i] - minimum[i]) / width[i];
                        }
                    }
                    
                    result[j] = value;
                }                
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
    public void testConstructor() throws NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        double[] startPoint = new double[ndims];
        minimum = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        width = new double[ndims];
        randomizer.fill(width, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        QuasiNewtonMultiOptimizer optimizer;
        
        //Test 1st constructor
        optimizer = new QuasiNewtonMultiOptimizer();
        assertNotNull(optimizer);
        
        assertEquals(optimizer.getTolerance(), 
                QuasiNewtonMultiOptimizer.DEFAULT_TOLERANCE, 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        assertEquals(optimizer.getIterations(), 0);
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            optimizer.getGradientListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(optimizer.isGradientListenerAvailable());
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isReady());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(optimizer.isLocked());
        
        
        //Test 2nd constructor
        optimizer = new QuasiNewtonMultiOptimizer(listener, gradientListener,
                tolerance);
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        assertEquals(optimizer.getIterations(), 0);
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isReady());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(optimizer.isLocked());
        
        //Force IllegalArgumentException
        optimizer = null;
        try {
            optimizer = new QuasiNewtonMultiOptimizer(listener, 
                    gradientListener, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(optimizer);
        
        
        //Test 3rd constructor
        optimizer = new QuasiNewtonMultiOptimizer(listener, gradientListener,
                startPoint, tolerance);
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 
                ABSOLUTE_ERROR);
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
        assertTrue(optimizer.isListenerAvailable());
        assertTrue(optimizer.isReady());
        assertFalse(optimizer.isResultAvailable());
        assertEquals(optimizer.getIterations(), 0);
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        assertFalse(optimizer.isLocked());
        
        //Force IllegalArgumentException
        optimizer = null;
        try {
            optimizer = new QuasiNewtonMultiOptimizer(listener, 
                    gradientListener, startPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(optimizer);        
    }
    
    @Test
    public void testIsReady() throws LockedException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] startPoint = new double[nDims];
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer();
        
        //optimizer will be ready when listener, gradient listener and start
        //point are provided
        assertFalse(optimizer.isReady());
        
        //set start point
        optimizer.setStartPoint(startPoint);
        
        //optimizer is not ready yet
        assertFalse(optimizer.isReady());
        
        //set listener
        optimizer.setListener(listener);
        
        //optimizer is not ready yet
        assertFalse(optimizer.isReady());
        
        //set gradient listener
        optimizer.setGradientListener(gradientListener);
        
        //optimizer is now ready
        assertTrue(optimizer.isReady());
    }
    
    @Test
    public void testGetSetTolerance() throws LockedException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer();
        
        //get default tolerance
        assertEquals(optimizer.getTolerance(), 
                QuasiNewtonMultiOptimizer.DEFAULT_TOLERANCE, 0.0);
        
        //set new tolerance
        optimizer.setTolerance(tolerance);
        
        //check correctness
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
    }
    
    @Test
    public void testGetSetStartPointAndAvailability() throws LockedException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] startPoint = new double[nDims];
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer();
        
        //initially start point is not available
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        //set start point
        optimizer.setStartPoint(startPoint);
        
        //check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 
                ABSOLUTE_ERROR);
        
        //set start point
        optimizer.setStartPoint(startPoint);
        
        //check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 
                ABSOLUTE_ERROR);
    }
    
    @Test
    public void testGetSetGradientListenerAndAvailability() 
            throws LockedException, NotAvailableException {
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer();
        
        //initially gradient listener is not avaialble
        assertFalse(optimizer.isGradientListenerAvailable());
        try {
            optimizer.getGradientListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        //set new listener
        optimizer.setGradientListener(gradientListener);
        
        //check correctness
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.getGradientListener(), gradientListener);
    }
    
    @Test
    public void testGetSetListenerAndAvailability() throws LockedException, 
        NotAvailableException {
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer();
        
        //initially listener is not available
        assertFalse(optimizer.isListenerAvailable());
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        //set new listener
        optimizer.setListener(listener);
        
        //check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
    }
    
    @Test
    public void testIsLocked() {
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer();
        
        assertFalse(optimizer.isLocked());
    }
    
    @Test
    public void testMinimize() throws Throwable {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] startPoint = new double[ndims];
        minimum = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        QuasiNewtonMultiOptimizer optimizer = new QuasiNewtonMultiOptimizer(
                listener, gradientListener, startPoint, 
                QuasiNewtonMultiOptimizer.DEFAULT_TOLERANCE);
        
        //check status
        assertTrue(optimizer.isReady());
        assertFalse(optimizer.isLocked());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException ignore) { }
        
        //minimize
        optimizer.minimize();
        
        //check status
        assertTrue(optimizer.isReady());
        assertFalse(optimizer.isLocked());
        assertTrue(optimizer.isResultAvailable());
        assertTrue(optimizer.getIterations() >= 0);
        
        //pick result
        double[] result = optimizer.getResult();
        double valueAtResult = optimizer.getEvaluationAtResult();
        
        //check correctness of result
        
        //compute function value at true minimum
        double valueAtMinimum = listener.evaluate(minimum);
        
        //compare values
        assertNotNull(result);
        assertEquals(valueAtResult, valueAtMinimum, ABSOLUTE_ERROR);
    }
}
