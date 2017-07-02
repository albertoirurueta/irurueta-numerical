/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.DerivativeEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 20, 2012
 */
package com.irurueta.numerical;

import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import org.junit.*;

public class DerivativeEstimatorTest 
    implements SingleDimensionFunctionEvaluatorListener{
    
    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;
    
    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    public static final double ABSOLUTE_ERROR = 1e-2;
    
    private double minimum;
    private double width;
    private double offset;
    
    public DerivativeEstimatorTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Override
    public double evaluate(double point) throws Throwable {
        return (point - minimum) * (point - minimum) / width + offset;
    }
    
    public double derivative(double x){
        return 2.0 * (x - minimum) / width;
    }
    
    @Test
    public void testConstructor(){                
        DerivativeEstimator estimator = new DerivativeEstimator(this);
        assertNotNull(estimator);
    }
    
    @Test
    public void testDerivative() throws EvaluationException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        double x = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);

        DerivativeEstimator estimator = new DerivativeEstimator(this);
        
        //estimate derivative
        double estDerivative = estimator.derivative(x);
        
        //real derivative
        double realDerivative = derivative(x);
        
        //compare both results
        assertEquals(estDerivative, realDerivative, ABSOLUTE_ERROR);
    }
}
