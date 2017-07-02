/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.optimization.PowellMultiOptimizer
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 25, 2012
 */
package com.irurueta.numerical.optimization;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Arrays;
import java.util.Random;
import org.junit.*;
import static org.junit.Assert.*;

public class PowellMultiOptimizerTest {
    
    public static final int MIN_DIMS = 2;
    public static final int MAX_DIMS = 4;
    
    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;
    
    public static final double MIN_DIRECTION = -1.0;
    public static final double MAX_DIRECTION = 1.0;
    
    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    public static final double MIN_TOLERANCE = 3e-8;
    public static final double MAX_TOLERANCE = 3e-6;
    
    public static final double ABSOLUTE_ERROR = 1e-6;
    
    public static final int TIMES = 100;
    
    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;
    
    private MultiDimensionFunctionEvaluatorListener listener;
    
    public PowellMultiOptimizerTest() {        
        listener = new MultiDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double[] point) throws Throwable {
                int dims = Math.min(Math.min(point.length, minimum.length), 
                        width.length);
                
                double value = 1.0;
                for (int i = 0; i < dims; i++) {
                    value *= Math.pow(point[i] - minimum[i], 2.0) / width[i];
                }
                
                value += offset;
                
                return value;
            }
        };
    }

    @BeforeClass
    public static void setUpClass() throws Exception { }

    @AfterClass
    public static void tearDownClass() throws Exception { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() throws WrongSizeException, 
            NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        minimum = new double[ndims];
        double[] point = new double[ndims];
        Matrix setsOfDirections = new Matrix(ndims, ndims);
        
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        PowellMultiOptimizer optimizer;
        
        //Test 1st constructor
        optimizer = new PowellMultiOptimizer();
        assertNotNull(optimizer);
        
        try {
            optimizer.getDirections();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.areDirectionsAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isDirectionAvailable());
        assertFalse(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), 
                PowellMultiOptimizer.DEFAULT_TOLERANCE, 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isLocked());
        
        
        //Test 2nd constructor
        optimizer = new PowellMultiOptimizer(listener, tolerance);
        assertNotNull(optimizer);
        
        try {
            optimizer.getDirections();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.areDirectionsAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isDirectionAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertEquals(optimizer.getListener(), listener);
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isLocked());
        
        //Force IllegalArgumentException
        optimizer = null;
        try {
            optimizer = new PowellMultiOptimizer(listener, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(optimizer);
        
        
        
        //Test 3rd constructor
        optimizer = new PowellMultiOptimizer(listener, point, setsOfDirections, 
                tolerance);
        assertNotNull(optimizer);
        
        assertEquals(optimizer.getDirections(), setsOfDirections);
        assertTrue(optimizer.areDirectionsAvailable());
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isDirectionAvailable());
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), point, 0.0);
        assertEquals(optimizer.getListener(), listener);
        assertTrue(optimizer.isListenerAvailable());
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isLocked());
    }
    
    @Test
    public void testGetSetDirectionsAndAvailability() 
            throws WrongSizeException, LockedException, NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        Matrix directions = new Matrix(nDims, nDims);
        double[] point = new double[nDims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        
        try {
            optimizer.getDirections();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.areDirectionsAvailable());
        
        //set directions
        optimizer.setPointAndDirections(point, directions);
        //check correctness
        assertEquals(optimizer.getDirections(), directions);
        assertTrue(optimizer.areDirectionsAvailable());        
    }
    
    @Test
    public void testGetSetDirection() throws LockedException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] point = new double[nDims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        double[] direction = new double[nDims];
        
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isDirectionAvailable());
        
        //set direction
        optimizer.setStartPointAndDirection(point, direction);
        
        //Check correctness
        assertArrayEquals(optimizer.getDirection(), direction, 0.0);
        assertTrue(optimizer.isDirectionAvailable());
    }
    
    @Test
    public void testIsReady() throws LockedException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int dims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] point = new double[dims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        double[] direction = new double[dims];
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        assertFalse(optimizer.isReady());
        
        //set listener
        optimizer.setListener(listener);
        assertFalse(optimizer.isReady());
        
        //set start point
        optimizer.setStartPointAndDirection(point, direction);
        assertTrue(optimizer.isReady());
    }
    
    @Test
    public void testGetSetTolerance() throws LockedException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        
        //get tolerance
        optimizer.setTolerance(tolerance);
        //check correctness
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        
        //Force IllegalArgumentException
        try {
            optimizer.setTolerance(-tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetStartPointAndAvailability() throws LockedException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int dims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] startPoint = new double[dims];
        double[] direction = new double[dims];
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        
        //get start point
        try {
            optimizer.getStartPoint();
            fail("NotAvailableExceptio expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isStartPointAvailable());
        
        //set start point
        optimizer.setStartPointAndDirection(startPoint, direction);
        //check correctness
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
        assertTrue(optimizer.isStartPointAvailable());
    }
    
    @Test
    public void testGetSetListenerAndAvailability() throws LockedException {
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        
        //get listener
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isListenerAvailable());
        
        //set listener
        optimizer.setListener(listener);
        assertTrue(optimizer.isListenerAvailable());
    }
    
    @Test
    public void testIsLocked() {
        
        PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
        
        assertFalse(optimizer.isLocked());
    }
    
    @Test
    public void testMinimize() throws LockedException, NotReadyException, 
        OptimizationException, NotAvailableException, Throwable {
        
        for (int i = 0; i < TIMES; i++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
            minimum = new double[ndims];
            double[] point = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
            double[] startPoint = Arrays.copyOf(point, ndims);
        
            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        
            PowellMultiOptimizer optimizer = new PowellMultiOptimizer();
            optimizer.setListener(listener);
            optimizer.setStartPoint(startPoint);
        
            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());
        
            //minimize
            optimizer.minimize();
        
            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());
        
            //get result
            double[] result = optimizer.getResult();
        
            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
            assertFalse(optimizer.isLocked());
        
            //check correctness of estimated result by comparing the difference of
            //function value at true minimum and the estimated minimum. If both are
            //very similar, then the algorithm will have converged
        
            //value at estimated result
            double value1 = listener.evaluate(result);
            //returned value by optimizer (must be equal to value1)
            double value2 = optimizer.getEvaluationAtResult();
            //value at true minimum
            double value3 = listener.evaluate(minimum);
            //value at start point
            double value4 = listener.evaluate(startPoint);
        
            assertEquals(value1, value2, 0.0);
            assertEquals(value1, value3, ABSOLUTE_ERROR);
        
            //check that indeed function has been minimized
            assertTrue(value1 <= value4);
        }
    }
}
