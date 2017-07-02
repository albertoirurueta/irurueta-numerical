/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.PolynomialEvaluation
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

public class DirectPolynomialEvaluationTest {
    
    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    
    public DirectPolynomialEvaluationTest() { }
    
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
        //test empty constructor
        DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation();
        
        //check initial values
        assertEquals(eval.getX(), 0.0, 0.0);
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.DIRECT_EVALUATION);
        
        //test constructor with values
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        
        eval = new DirectPolynomialEvaluation(x, evaluation);
        
        //check correctness
        assertEquals(eval.getX(), x, 0.0);
        assertEquals(eval.getEvaluation(), evaluation, 0.0);        
    }
    
    @Test
    public void testGetSetX() {
        DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation();
        
        //check initial value
        assertEquals(eval.getX(), 0.0, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setX(x);
        
        //check correctness
        assertEquals(eval.getX(), x, 0.0);
    }
    
    @Test
    public void testGetSetEvaluation() {
        DirectPolynomialEvaluation eval = new DirectPolynomialEvaluation();
        
        //check initial value
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        eval.setEvaluation(evaluation);
        
        //check correctness
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
    }
}
