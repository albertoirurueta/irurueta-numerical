/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.DirectionalEvaluator
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

public class DirectionalEvaluatorTest {
    
    public static final int MIN_DIMS = 2;
    public static final int MAX_DIMS = 5;
    
    public static final double MIN_EVAL_POINT = -1e3;
    public static final double MAX_EVAL_POINT = 1e3;
    
    public static final double MIN_DIRECTION = -1.0;
    public static final double MAX_DIRECTION = 1.0;
    
    public static final double MIN_OFFSET = -1e3;
    public static final double MAX_OFFSET = 1e3;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;
    
    private MultiDimensionFunctionEvaluatorListener listener;
    
    
    public DirectionalEvaluatorTest() {
        listener = new MultiDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double[] point) throws Throwable {
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

    
    @Test
    public void testConstructor(){
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        minimum = new double[ndims];
        double[] point = new double[ndims];
        double[] direction = new double[ndims];
        width = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        DirectionalEvaluator evaluator = new DirectionalEvaluator(listener,
                point, direction);
        assertNotNull(evaluator);
    }
    
    @Test
    public void testEvaluateAt() throws EvaluationException, Throwable{

        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        minimum = new double[ndims];
        double[] point = new double[ndims];
        double[] direction = new double[ndims];
        width = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        double x = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);
        
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        
        DirectionalEvaluator evaluator = new DirectionalEvaluator(listener,
                point, direction);
        
        double value = evaluator.evaluateAt(x);
        
        //check correctness
        double[] xt = new double[ndims];
        for(int i = 0; i < ndims; i++){
            xt[i] = point[i] + x * direction[i];
        }
        
        double value2 = listener.evaluate(xt);
        
        assertEquals(value, value2, 0.0);
    }
}
