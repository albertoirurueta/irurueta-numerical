/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.optimization.BrentSingleOptimizer
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 17, 2012
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

public class BrentSingleOptimizerTest 
    implements SingleDimensionFunctionEvaluatorListener {
    
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
    
    public BrentSingleOptimizerTest() { }

    @BeforeClass
    public static void setUpClass() throws Exception { }

    @AfterClass
    public static void tearDownClass() throws Exception { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Override
    public double evaluate(double point) throws Throwable {
        return (point - minimum) * (point - minimum) / width + offset;
    }
    
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
        
        BrentSingleOptimizer optimizer;
        
        
        //test 1st constructor
        optimizer = new BrentSingleOptimizer();
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), 
                BrentSingleOptimizer.DEFAULT_TOLERANCE, 0.0);
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
        
        
        //test 2nd constructor
        optimizer = new BrentSingleOptimizer(minEvalPoint, middleEvalPoint,
                maxEvalPoint, tolerance);
        assertNotNull(optimizer);
        assertEquals(optimizer.getTolerance(), tolerance, 0.0);
        assertFalse(optimizer.isReady());
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
        
        //Force InvalidBracketRangeException
        optimizer = null;
        try {
            optimizer = new BrentSingleOptimizer(maxEvalPoint, middleEvalPoint,
                    minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(minEvalPoint, maxEvalPoint,
                    middleEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(maxEvalPoint, minEvalPoint, 
                    middleEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(middleEvalPoint, minEvalPoint,
                    maxEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(middleEvalPoint, maxEvalPoint,
                    minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        
        //Force IllegalArgumentException
        try {
            optimizer = new BrentSingleOptimizer(minEvalPoint, middleEvalPoint,
                    maxEvalPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(optimizer);
        
        
        
        //test 3rd constructor
        optimizer = new BrentSingleOptimizer(this, minEvalPoint, 
                middleEvalPoint, maxEvalPoint, tolerance);
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
        assertEquals(optimizer.getListener(), this);
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
        
        //Force InvalidBracketRangeException
        optimizer = null;
        try {
            optimizer = new BrentSingleOptimizer(this, maxEvalPoint, 
                    middleEvalPoint, minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(this, minEvalPoint, 
                    maxEvalPoint, middleEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(this, maxEvalPoint, 
                    minEvalPoint, middleEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(this, middleEvalPoint, 
                    minEvalPoint, maxEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        try {
            optimizer = new BrentSingleOptimizer(this, middleEvalPoint, 
                    maxEvalPoint, minEvalPoint, tolerance);
            fail("InvalidBracketRangeException expected but not thrown");
        } catch (InvalidBracketRangeException e) { }
        assertNull(optimizer);
        
        //Force IllegalArgumentException
        try {
            optimizer = new BrentSingleOptimizer(this, minEvalPoint, 
                    middleEvalPoint, maxEvalPoint, -tolerance);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(optimizer);        
    }
    
    @Test
    public void testGetSetTolerance() throws LockedException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double tolerance = randomizer.nextDouble(MIN_TOLERANCE, MAX_TOLERANCE);
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
        
        assertEquals(optimizer.getTolerance(), 
                BrentSingleOptimizer.DEFAULT_TOLERANCE, 0.0);
        
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
    public void testGetSetBracketAndAvailability() throws NotAvailableException,
        LockedException, InvalidBracketRangeException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double minEvalPoint = randomizer.nextDouble(MIN_EVAL_POINT, 
                MAX_EVAL_POINT);
        double middleEvalPoint = randomizer.nextDouble(minEvalPoint,
                MAX_EVAL_POINT);
        double maxEvalPoint = randomizer.nextDouble(middleEvalPoint,
                MAX_EVAL_POINT);
        
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
        
        assertTrue(optimizer.isBracketAvailable());
        
        assertEquals(optimizer.getMinEvaluationPoint(), 
                BracketedSingleOptimizer.DEFAULT_MIN_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMiddleEvaluationPoint(),
                BracketedSingleOptimizer.DEFAULT_MIDDLE_EVAL_POINT, 0.0);
        assertEquals(optimizer.getMaxEvaluationPoint(),
                BracketedSingleOptimizer.DEFAULT_MAX_EVAL_POINT, 0.0);
        
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
        
        //set minimum to be estimated
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
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer(this,
                minEvalPoint, middleEvalPoint, maxEvalPoint, 
                BrentSingleOptimizer.DEFAULT_TOLERANCE);
        
        //attempting to retrieve evaluation fails because altough bracket is
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
        assertEquals(optimizer.getEvaluationAtMin(), evaluate(minEvalPoint), 
                0.0);
        assertEquals(optimizer.getEvaluationAtMiddle(), evaluate(
                middleEvalPoint), 0.0);
        assertEquals(optimizer.getEvaluationAtMax(), evaluate(maxEvalPoint), 
                0.0);
    }
    
    @Test
    public void testComputeBracket() throws LockedException, NotReadyException, 
        OptimizationException, InvalidBracketRangeException, 
        NotAvailableException {
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
        optimizer.setListener(this);
        
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
        //ax < minimum and cx > minimum
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
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
        
        assertFalse(optimizer.isListenerAvailable());
        try{
            optimizer.getListener();
        }catch(NotAvailableException e){}
        assertFalse(optimizer.isReady());
        
        //set listener
        optimizer.setListener(this);
        
        //check correctness
        assertTrue(optimizer.isListenerAvailable());
        assertEquals(optimizer.getListener(), this);
        assertTrue(optimizer.isReady());
    }
    
    @Test
    public void testIsLocked() {
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
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
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
        
        assertFalse(optimizer.isReady());
        
        //set listener
        optimizer.setListener(this);
        
        //check correctness
        assertTrue(optimizer.isReady());
    }
    
    @Test
    public void testMinimizeGetResultAndAvailability() throws LockedException, 
            InvalidBracketRangeException, NotReadyException, 
            OptimizationException, NotAvailableException, Throwable {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        //set minimum to be estimated
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        //set offset
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        //set width
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        BrentSingleOptimizer optimizer = new BrentSingleOptimizer();
        optimizer.setListener(this);
        
        assertFalse(optimizer.isResultAvailable());
        try {
            optimizer.getResult();
            fail("NotAvailableException expected but not thrown");
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
        assertEquals(optimizer.getEvaluationAtResult(), evaluate(
                optimizer.getResult()), 0.0);
    }
}
