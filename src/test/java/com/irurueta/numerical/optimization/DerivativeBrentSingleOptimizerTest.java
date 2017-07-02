/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.optimization.DerivativeBrentSingleOptimizer
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 19, 2012
 */
package com.irurueta.numerical.optimization;

import com.irurueta.numerical.InvalidBracketRangeException;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.*;
import org.junit.*;

public class DerivativeBrentSingleOptimizerTest {
    
    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;
    
    public static final double MIN_TOLERANCE = 3e-8;
    public static final double MAX_TOLERANCE = 1e-5;
    
    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    private double minimum;
    private double offset;
    private double width;
    
    private SingleDimensionFunctionEvaluatorListener listener;
    private SingleDimensionFunctionEvaluatorListener derivativeListener;
    
    public DerivativeBrentSingleOptimizerTest() {
        listener = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) throws Throwable {
                return (point - minimum) * (point - minimum) / width + offset;
            }
        };
          
        derivativeListener = new SingleDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double point) throws Throwable {
                return 2.0 / width * (point - minimum);
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
    public void testConstructor() throws NotAvailableException, 
        InvalidBracketRangeException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double middleEvalPoint = randomizer.nextDouble(minEvalPoint, 
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(middleEvalPoint, 
                MAX_EVAL_POINT);
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        DerivativeBrentSingleOptimizer optimizer;
        
        //testing 1st constructor
        optimizer = new DerivativeBrentSingleOptimizer();
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), 
                DerivativeBrentSingleOptimizer.DEFAULT_TOLERANCE, 0.0);
        assertFalse(optimizer.isReady());
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(optimizer.getMinEvaluationPoint(), 
                BracketedSingleOptimizer.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(),
                BracketedSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(),
                BracketedSingleOptimizer.DEFAULT_MAX_EVAL_POINT, 0.0);
        try {
            optimizer.getEvaluationAtMin();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtMiddle();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtMax();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.areBracketEvaluationsAvailable());
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
        try {
            optimizer.getListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.isListenerAvailable());
        
        
        //testing 2nd constructor
        optimizer = new DerivativeBrentSingleOptimizer(listener, 
                derivativeListener, minEvalPoint, middleEvalPoint, maxEvalPoint,
                tolerance);
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertTrue(optimizer.isReady());
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(optimizer.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(), middleEvalPoint, 
                0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        try {
            optimizer.getEvaluationAtMin();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtMiddle();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtMax();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.areBracketEvaluationsAvailable());
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
        assertEquals(optimizer.getDerivativeListener(), derivativeListener);
        assertTrue(optimizer.isDerivativeListenerAvailable());
        
        //Force InvalidBracketRangeException
        optimizer = null;
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, maxEvalPoint, middleEvalPoint,
                    minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, minEvalPoint, maxEvalPoint, 
                    middleEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, maxEvalPoint, minEvalPoint,
                    middleEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, middleEvalPoint, minEvalPoint, 
                    maxEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, middleEvalPoint, maxEvalPoint,
                    minEvalPoint, tolerance);
            fail("IsnvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }

        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, maxEvalPoint, middleEvalPoint,
                    minEvalPoint, -tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, minEvalPoint, maxEvalPoint, 
                    middleEvalPoint, -tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, maxEvalPoint, minEvalPoint,
                    middleEvalPoint, -tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, middleEvalPoint, minEvalPoint, 
                    maxEvalPoint, -tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, middleEvalPoint, maxEvalPoint,
                    minEvalPoint, -tolerance);
            fail("IsnvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        
        
        //Force IllegalArgumentException
        try {
            optimizer = new DerivativeBrentSingleOptimizer(listener,
                    derivativeListener, minEvalPoint, middleEvalPoint,
                    maxEvalPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(optimizer);
    }
    
    @Test
    public void testGetSetTolerance() throws LockedException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer();
        
        assertEquals(optimizer.getTolerance(),
                DerivativeBrentSingleOptimizer.DEFAULT_TOLERANCE, 0.0);
        
        //set new tolerance
        optimizer.setTolerance(tolerance);
        
        //check corectness
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        
        //Force IllegalArgumentException
        try {
            optimizer.setTolerance(-tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetBracketAndAvailability() throws NotAvailableException,
        LockedException, InvalidBracketRangeException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double middleEvalPoint = randomizer.nextDouble(minEvalPoint, 
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(middleEvalPoint, 
                MAX_EVAL_POINT);
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer();
        
        assertTrue(optimizer.isBracketAvailable());
        assertEquals(optimizer.getMinEvaluationPoint(),
                DerivativeBrentSingleOptimizer.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(),
                DerivativeBrentSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(),
                DerivativeBrentSingleOptimizer.DEFAULT_MAX_EVAL_POINT, 0.0);
        
        //set new bracket
        optimizer.setBracket(minEvalPoint, middleEvalPoint, maxEvalPoint);
        
        //check correctness
        assertEquals(optimizer.getMinEvaluationPoint(), minEvalPoint, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(), middleEvalPoint, 
                0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(), maxEvalPoint, 0.0);
        
        //Force InvalidBracketRangeException
        try {
            optimizer.setBracket(maxEvalPoint, middleEvalPoint, minEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer.setBracket(minEvalPoint, maxEvalPoint, middleEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer.setBracket(maxEvalPoint, minEvalPoint, middleEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer.setBracket(middleEvalPoint, minEvalPoint, maxEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        try {
            optimizer.setBracket(middleEvalPoint, maxEvalPoint, minEvalPoint);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
    }
    
    @Test
    public void testGetEvaluationsAndEvaluateBracket() 
            throws InvalidBracketRangeException, LockedException, 
            NotReadyException, OptimizationException, NotAvailableException, 
            Throwable {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set bracket
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double middleEvalPoint = randomizer.nextDouble(minEvalPoint, 
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(middleEvalPoint,
                MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer(listener, derivativeListener,
                minEvalPoint, middleEvalPoint, maxEvalPoint, 
                DerivativeBrentSingleOptimizer.DEFAULT_TOLERANCE);
        
        //attempting to retrieve evaluations fails because although bracket is
        //available, evaluations have not yet been computed
        try {
            optimizer.getEvaluationAtMin();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtMiddle();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtMax();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        assertFalse(optimizer.areBracketEvaluationsAvailable());
        
        //we compute evaluations
        optimizer.evaluateBracket();
        assertTrue(optimizer.areBracketEvaluationsAvailable());
        
        //check correctness
        assertEquals(optimizer.getEvaluationAtMin(), 
                listener.evaluate(minEvalPoint), 0.0);
        assertEquals(optimizer.getEvaluationAtMiddle(),
                listener.evaluate(middleEvalPoint), 0.0);
        assertEquals(optimizer.getEvaluationAtMax(),
                listener.evaluate(maxEvalPoint), 0.0);
    }
    
    @Test
    public void testComputeBracket() throws InvalidBracketRangeException, 
        LockedException, NotReadyException, OptimizationException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer(listener, derivativeListener,
                DerivativeBrentSingleOptimizer.DEFAULT_MIN_EVAL_POINT,
                DerivativeBrentSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT,
                DerivativeBrentSingleOptimizer.DEFAULT_MAX_EVAL_POINT,
                DerivativeBrentSingleOptimizer.DEFAULT_TOLERANCE);
        
        assertTrue(optimizer.isBracketAvailable());
        assertFalse(optimizer.areBracketEvaluationsAvailable());

        optimizer.computeBracket();
        
        assertTrue(optimizer.isBracketAvailable());
        assertTrue(optimizer.areBracketEvaluationsAvailable());
        
        //after computing bracket we can only ensure that ax < bx < cx and also
        //that fa > fb and fc > fb
        assertTrue(optimizer.getMinEvaluationPoint() <= 
                optimizer.getMiddleEvaluationPoint());
        assertTrue(optimizer.getMiddleEvaluationPoint() <=
                optimizer.getMaxEvaluationPoint());
        
        assertTrue(optimizer.getEvaluationAtMin() >=
                optimizer.getEvaluationAtMiddle());
        assertTrue(optimizer.getEvaluationAtMax() >=
                optimizer.getEvaluationAtMiddle());
        
        //also bracket limits must surround the real minimum location, that is
        //ax < minimum and cv > minimum
        assertTrue(optimizer.getMinEvaluationPoint() <= minimum);
        assertTrue(optimizer.getMaxEvaluationPoint() >= minimum);
    }
    
    @Test
    public void testGetSetListenerAndAvailability() throws LockedException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer();
        
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
    public void testGetSetDerivativeListenerAndAvailability() 
            throws LockedException, NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer();
        
        assertFalse(optimizer.isDerivativeListenerAvailable());
        try {
            optimizer.getDerivativeListener();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        
        //set listener
        optimizer.setDerivativeListener(derivativeListener);
        
        //check correctness
        assertTrue(optimizer.isDerivativeListenerAvailable());
        assertEquals(optimizer.getDerivativeListener(), derivativeListener);        
    }
    
    @Test
    public void testIsLocked() {
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer();
        assertFalse(optimizer.isLocked());
    }
    
    @Test
    public void testIsReady() throws LockedException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer();
        assertFalse(optimizer.isReady());
        
        //set listener
        optimizer.setListener(listener);
        
        assertFalse(optimizer.isReady());
        
        //set derivative listener
        optimizer.setDerivativeListener(derivativeListener);
        
        //check correctness
        assertTrue(optimizer.isReady());
    }
    
    @Test
    public void testMinimizeGetResultAndAvailability() 
            throws InvalidBracketRangeException, LockedException, 
            NotReadyException, OptimizationException, NotAvailableException, 
            Throwable {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        DerivativeBrentSingleOptimizer optimizer = 
                new DerivativeBrentSingleOptimizer(listener, derivativeListener,
                DerivativeBrentSingleOptimizer.DEFAULT_MIN_EVAL_POINT,
                DerivativeBrentSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT,
                DerivativeBrentSingleOptimizer.DEFAULT_MAX_EVAL_POINT,
                DerivativeBrentSingleOptimizer.DEFAULT_TOLERANCE);
        
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
        } catch (NotAvailableException e) { }
        try {
            optimizer.getEvaluationAtResult();
            fail("NotAvailabelException expected but not thrown");
        } catch (NotAvailableException e) { }
        
        //minimize
        if (minimum > 0.0) {
            optimizer.setBracket(0.8 * minimum, 1.1 * minimum, 1.5 * minimum);
        } else {
            optimizer.setBracket(1.5 * minimum, 1.1 * minimum, 0.8 * minimum);
        }
        
        optimizer.minimize();
        
        //test correctness
        assertTrue(optimizer.isResultAvailable());
        assertEquals(optimizer.getResult(), minimum, optimizer.getTolerance());
        assertEquals(optimizer.getEvaluationAtResult(), listener.evaluate(
                optimizer.getResult()), 0.0);
    }
}
