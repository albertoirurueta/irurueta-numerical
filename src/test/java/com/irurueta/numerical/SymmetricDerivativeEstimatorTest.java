/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.SymmetricDerivativeEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 28, 2012
 */
package com.irurueta.numerical;

import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import org.junit.*;

public class SymmetricDerivativeEstimatorTest 
    implements SingleDimensionFunctionEvaluatorListener{
    
    public static final double MIN_EVAL_POINT = -10.0;
    public static final double MAX_EVAL_POINT = 10.0;
    
    public static final double MIN_OFFSET = -10.0;
    public static final double MAX_OFFSET = 10.0;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    public static final double ABSOLUTE_ERROR = 1e-3;
    
    public static final int TIMES = 100;
    
    private double minimum;
    private double width;
    private double offset;
    
    public SymmetricDerivativeEstimatorTest() {
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
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);        
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);        
        width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
        SymmetricDerivativeEstimator estimator = 
                new SymmetricDerivativeEstimator(this);
        assertNotNull(estimator);
    }
    
    @Test
    public void testDerivative() throws EvaluationException{
        
        for(int i = 0; i < TIMES; i++){
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            minimum = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
            width = randomizer.nextDouble(MIN_WIDTH, MAX_WIDTH);
        
            double x = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
            SymmetricDerivativeEstimator estimator =
                    new SymmetricDerivativeEstimator(this);
        
            //estimate derivative
            double estDerivative = estimator.derivative(x);
        
            //real derivative
            double realDerivative = derivative(x);
        
            //compare both results
            assertEquals(estDerivative, realDerivative, ABSOLUTE_ERROR);
        }
    }
}
