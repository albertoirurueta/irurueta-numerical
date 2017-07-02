/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.optimization.DerivativeConjugateGradientMultiOptimizer
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 20, 2012
 */
package com.irurueta.numerical.optimization;

import com.irurueta.numerical.GradientFunctionEvaluatorListener;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.*;
import org.junit.*;

public class DerivativeConjugateGradientMultiOptimizerTest {
    
    public static final int MIN_DIMS = 2;
    public static final int MAX_DIMS = 4;
    
    public static final double MIN_EVAL_POINT = -10.0;
    public static final double MAX_EVAL_POINT = 10.0;
    
    public static final double MIN_OFFSET = 0.0;
    public static final double MAX_OFFSET = 1.0;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    public static final double MIN_TOLERANCE = 3e-8;
    public static final double MAX_TOLERANCE = 3e-6;
    
    public static final double ABSOLUTE_ERROR = 1e-3;
    
    public static final int PROBABILITY = 50;
    
    public static final double MIN_DIRECTION = -1.0;
    public static final double MAX_DIRECTION = 1.0;
    
    public static final int TIMES = 100;
    
    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;
    
    private MultiDimensionFunctionEvaluatorListener listener;
    private GradientFunctionEvaluatorListener gradientListener;
    
    public DerivativeConjugateGradientMultiOptimizerTest() {
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
        
        gradientListener = new GradientFunctionEvaluatorListener() {

            @Override
            public void evaluateGradient(double[] params, double[] result) 
                    throws Throwable {
                int dims = Math.min(Math.min(params.length, minimum.length), 
                        width.length);
                
                double value;
                for (int j = 0; j < dims; j++) {
                    value = 1.0;
                    for (int i = 0; i < dims; i++) {
                        if (i != j) {
                            value *= Math.pow(params[i] - minimum[i], 2.0) / 
                                    width[i];
                        } else {
                            value *= 2.0 * (params[i] - minimum[i]) / width[i];
                        }
                    }
                    
                    result[j] = value;
                }
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
    public void testConstructor() throws NotAvailableException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        boolean usePolakRibiere = randomizer.nextBoolean();
        
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        double[] startPoint = new double[ndims];
        minimum = new double[ndims];
        width = new double[ndims];
        double[] direction = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);
        
        DerivativeConjugateGradientMultiOptimizer optimizer;
        
        //test 1st constructor
        optimizer = new DerivativeConjugateGradientMultiOptimizer();
        assertNotNull(optimizer);
        assertFalse(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), 
                DerivativeConjugateGradientMultiOptimizer.DEFAULT_TOLERANCE, 
                0.0);
        try {
            optimizer.getGradientListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.isPolakRibiereEnabled(), 
                DerivativeConjugateGradientMultiOptimizer.
                DEFAULT_USE_POLAK_RIBIERE);
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isDirectionAvailable());
        try {
            optimizer.getDirection();
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
        
        
        
        //test 2nd constructor
        optimizer = new DerivativeConjugateGradientMultiOptimizer(listener,
                gradientListener, startPoint, direction, tolerance, 
                usePolakRibiere);
        assertNotNull(optimizer);
        assertTrue(optimizer.isReady());
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
        assertEquals(optimizer.isPolakRibiereEnabled(), usePolakRibiere);
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
        assertTrue(optimizer.isDirectionAvailable());
        assertArrayEquals(optimizer.getDirection(), direction, 0.0);
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
        
        
        //Force IllegalArgumentException
        optimizer = null;
        try {
            optimizer = new DerivativeConjugateGradientMultiOptimizer(listener,
                    gradientListener, startPoint, direction, -tolerance, 
                    usePolakRibiere);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(optimizer);
    }
    
    @Test
    public void testIsReady() throws LockedException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] startPoint = new double[ndims];
        minimum = new double[ndims];
        width = new double[ndims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        assertFalse(optimizer.isReady());
        
        //set listener
        optimizer.setListener(listener);
        
        assertFalse(optimizer.isReady());
        
        //set gradient listener
        optimizer.setGradientListener(gradientListener);
        
        assertFalse(optimizer.isReady());
        
        //set start point
        optimizer.setStartPoint(startPoint);
        
        //once both listeners and start point are available, optimizer becomes
        //ready
        assertTrue(optimizer.isReady());
    }
    
    @Test
    public void testGetSetTolerance() throws LockedException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        //get tolerance
        assertEquals(optimizer.getTolerance(),
                DerivativeConjugateGradientMultiOptimizer.DEFAULT_TOLERANCE, 
                0.0);
        
        //set new tolerance
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
    public void testGetSetGradientListenerAndAvailability() 
        throws LockedException, NotAvailableException {
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        //Get gradient listener
        try {
            optimizer.getGradientListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isGradientListenerAvailable());
        
        //set gradient listener
        optimizer.setGradientListener(gradientListener);
        
        //check correctness
        assertEquals(optimizer.getGradientListener(), gradientListener);
        assertTrue(optimizer.isGradientListenerAvailable());
    }
    
    @Test
    public void testGetSetUsePolakRibiere() throws LockedException {
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        //get initial status
        assertEquals(optimizer.isPolakRibiereEnabled(),
                DerivativeConjugateGradientMultiOptimizer.
                DEFAULT_USE_POLAK_RIBIERE);
        
        //disable
        optimizer.setUsePolakRibiere(false);
        
        //check correctness
        assertFalse(optimizer.isPolakRibiereEnabled());
        
        //enable
        optimizer.setUsePolakRibiere(true);
        
        //check correctness
        assertTrue(optimizer.isPolakRibiereEnabled());
    }
    
    @Test
    public void testGetSetStartPointAndAvailability() throws LockedException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        double[] startPoint = new double[nDims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        //get start point
        assertFalse(optimizer.isStartPointAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        
        //set start point
        optimizer.setStartPoint(startPoint);
        
        //check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
    }
    
    @Test
    public void testGetSetStartPointAndDirectionAndAvailability() 
            throws LockedException, NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int nDims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);

        double[] startPoint = new double[nDims];
        randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);        
        double[] direction = new double[nDims];
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        assertFalse(optimizer.isStartPointAvailable());
        assertFalse(optimizer.isDirectionAvailable());
        try {
            optimizer.getStartPoint();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getDirection();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        
        //set start point and direction
        optimizer.setStartPointAndDirection(startPoint, direction);
        
        //check correctness
        assertTrue(optimizer.isStartPointAvailable());
        assertTrue(optimizer.isDirectionAvailable());
        assertArrayEquals(optimizer.getStartPoint(), startPoint, 0.0);
        assertArrayEquals(optimizer.getDirection(), direction, 0.0);
        
        //Force IllegalArgumentException
        double[] wrongStartPoint = new double[nDims + 1];
        try {
            optimizer.setStartPointAndDirection(wrongStartPoint, direction);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetListenerAndAvailability() throws LockedException, 
        NotAvailableException {
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        //check initial status
        assertFalse(optimizer.isListenerAvailable());
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        
        //set listener
        optimizer.setListener(listener);
        
        //check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), listener);
    }
    
    @Test
    public void testIsLocked() {
        
        DerivativeConjugateGradientMultiOptimizer optimizer =
                new DerivativeConjugateGradientMultiOptimizer();
        
        assertFalse(optimizer.isLocked());
    }
    
    @Test
    public void testMinimize() throws LockedException, NotReadyException, 
        OptimizationException, NotAvailableException, Throwable {
        
        //we repeat the process because depending on the start point this
        //algorithm is not very accurate
        for (int i = 0; i < TIMES; i++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
            double[] startPoint = new double[ndims];
            minimum = new double[ndims];
            width = new double[ndims];
            randomizer.fill(startPoint, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);            
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
            DerivativeConjugateGradientMultiOptimizer optimizer =
                    new DerivativeConjugateGradientMultiOptimizer();
            optimizer.setListener(listener);
            optimizer.setGradientListener(gradientListener);
            optimizer.setStartPoint(startPoint);
        
            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertFalse(optimizer.isResultAvailable());        
            try {
                optimizer.getResult();
                fail("NotAvailableException expected but not thrown");
            } catch (NotAvailableException e) { }
        
            //minimize
            optimizer.minimize();
        
            //check correctness of result
            assertFalse(optimizer.isLocked());
            assertTrue(optimizer.isReady());
            assertTrue(optimizer.isResultAvailable());
        
            double[] result = optimizer.getResult();
        
            //check that function at estimated result and at tue minimum have 
            //almost the same value
            double valueMin = listener.evaluate(minimum);
            double valueResult = listener.evaluate(result);
        
            assertEquals(valueMin, valueResult, ABSOLUTE_ERROR);
        }
    }
}
