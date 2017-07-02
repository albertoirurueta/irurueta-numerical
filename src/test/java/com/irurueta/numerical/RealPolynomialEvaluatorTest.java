/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.RealPolynomialEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 9, 2016.
 */
package com.irurueta.numerical;

import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class RealPolynomialEvaluatorTest {
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    public static final int MIN_LENGTH = 1;
    public static final int MAX_LENGTH = 5;
    
    public RealPolynomialEvaluatorTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}

    @Test
    public void testConstructor(){
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
        
        double[] polyParams = new double[length];
        
        RealPolynomialEvaluator evaluator = new RealPolynomialEvaluator(
                polyParams);
        
        //check correctness
        assertSame(evaluator.getPolyParams(), polyParams);        
    }
    
    @Test
    public void testEvaluate(){
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
        
        double[] polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        RealPolynomialEvaluator evaluator = new RealPolynomialEvaluator(
                polyParams);
        
        assertEquals(evaluator.evaluate(x), 
                PolynomialEvaluator.evaluate(polyParams, x), 0.0);
    }
}
