/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.SymmetricGradientEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 28, 2012
 */
package com.irurueta.numerical;

import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import static org.junit.Assert.*;
import org.junit.*;

public class SymmetricGradientEstimatorTest 
    implements MultiDimensionFunctionEvaluatorListener{
    
    public static final int MIN_DIMS = 1;
    public static final int MAX_DIMS = 3;
    
    public static final double MIN_EVAL_POINT = -10.0;
    public static final double MAX_EVAL_POINT = 10.0;
    
    public static final double MIN_OFFSET = -10.0;
    public static final double MAX_OFFSET = 10.0;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 2.0;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    
    public static final int TIMES = 100;
    
    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;
    
    public SymmetricGradientEstimatorTest() {
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
    
    double[] gradient(double[] params){        
        int dims = Math.min(Math.min(params.length, minimum.length), 
                width.length);
        
        double[] gradient = new double[dims];
        
        double value;
        for(int j = 0; j < dims; j++){
            value = 1.0;
            for(int i = 0; i < dims; i++){
                if(i != j){
                    value *= Math.pow(params[i] - minimum[i], 2.0) / width[i];
                }else{
                    value *= 2.0 * (params[i] - minimum[i]) / width[i];
                }
            }
            
            gradient[j] = value;
        }
        
        return gradient;
    }
    
    @Test
    public void testConstructor(){
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        minimum = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        double[] point = new double[ndims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        width = new double[ndims];
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        
        SymmetricGradientEstimator estimator = new SymmetricGradientEstimator(
                this);
        assertNotNull(estimator);
    }
    
    @Test
    public void testGradient() throws EvaluationException{
        
        for(int t = 0; t < TIMES; t++){
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
            minimum = new double[ndims];
            randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
            double[] point = new double[ndims];
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
            offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
            width = new double[ndims];
            randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        
            SymmetricGradientEstimator estimator = new SymmetricGradientEstimator(
                    this);

            double[] gradient1 = estimator.gradient(point);
            double[] gradient2 = new double[ndims];
            estimator.gradient(point, gradient2);
        
            //check correctness
            double[] gradient3 = gradient(point);
        
            assertEquals(gradient1.length, ndims);
            assertEquals(gradient2.length, ndims);
            assertEquals(gradient3.length, ndims);
            for(int i = 0; i < ndims; i++){
                assertEquals(gradient1[i], gradient3[i], ABSOLUTE_ERROR);
                assertEquals(gradient2[i], gradient3[i], ABSOLUTE_ERROR);
            }
        }
    }
}
