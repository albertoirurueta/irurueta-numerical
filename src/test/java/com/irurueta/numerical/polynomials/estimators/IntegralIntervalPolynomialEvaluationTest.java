/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.IntegralIntervalPolynomialEvaluation
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class IntegralIntervalPolynomialEvaluationTest {

    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    
    public static final int MIN_ORDER = 1;
    public static final int MAX_ORDER = 5;
    
    public IntegralIntervalPolynomialEvaluationTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() {
        IntegralIntervalPolynomialEvaluation eval = 
                new IntegralIntervalPolynomialEvaluation();
        
        //check default values
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        assertEquals(eval.getStartX(), 0.0, 0.0);
        assertEquals(eval.getEndX(), 0.0, 0.0);
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_INTERVAL);
        
        //test constructor with values
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);        
        double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        
        eval = new IntegralIntervalPolynomialEvaluation(startX, endX, 
                evaluation, order);
        
        //check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getStartX(), startX, 0.0);
        assertEquals(eval.getEndX(), endX, 0.0);
        assertEquals(eval.getIntegralOrder(), order);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_INTERVAL);
        
        
        eval = new IntegralIntervalPolynomialEvaluation(startX, endX, 
                evaluation);
        
        //check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getStartX(), startX, 0.0);
        assertEquals(eval.getEndX(), endX, 0.0);
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_INTERVAL);        
    }
    
    @Test
    public void testGetSetStartX() {
        IntegralIntervalPolynomialEvaluation eval = 
                new IntegralIntervalPolynomialEvaluation();
        
        //check default value
        assertEquals(eval.getStartX(), 0.0, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double startX = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        eval.setStartX(startX);
        
        //check correctness
        assertEquals(eval.getStartX(), startX, 0.0);
    }
    
    @Test
    public void testGetSetEndX() {
        IntegralIntervalPolynomialEvaluation eval =
                new IntegralIntervalPolynomialEvaluation();
        
        //check default value
        assertEquals(eval.getEndX(), 0.0, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double endX = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setEndX(endX);
        
        //check correctness
        assertEquals(eval.getEndX(), endX, 0.0);
    }
    
    @Test
    public void testGetSetIntegralOrder() {
        IntegralIntervalPolynomialEvaluation eval =
                new IntegralIntervalPolynomialEvaluation();
        
        //check default value
        assertEquals(eval.getIntegralOrder(), 1);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setIntegralOrder(order);
        
        //check correctness
        assertEquals(eval.getIntegralOrder(), order);
    }
}
