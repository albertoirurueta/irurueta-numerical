/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.IntegralPolynomialEvaluation
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

public class IntegralPolynomialEvaluationTest {

    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    
    public static final int MIN_ORDER = 1;
    public static final int MAX_ORDER = 5;
    
    public IntegralPolynomialEvaluationTest() { }
    
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
        IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();
        
        //check default values
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        assertEquals(eval.getX(), 0.0, 0.0);
        assertNull(eval.getConstants());
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_EVALUATION);
        
        //test constructor with values
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        double[] constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        
        eval = new IntegralPolynomialEvaluation(x, evaluation, constants, 
                order);
        
        //check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertSame(eval.getConstants(), constants);
        assertEquals(eval.getIntegralOrder(), order);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_EVALUATION);
                
        
        eval = new IntegralPolynomialEvaluation(x, evaluation, constants);
        
        //check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertSame(eval.getConstants(), constants);
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_EVALUATION);
        
        
        eval = new IntegralPolynomialEvaluation(x, evaluation, order);
        
        //check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertNull(eval.getConstants());
        assertEquals(eval.getIntegralOrder(), order);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_EVALUATION);
        
        
        eval = new IntegralPolynomialEvaluation(x, evaluation);
        
        //check default values
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
        assertEquals(eval.getX(), x, 0.0);
        assertNull(eval.getConstants());
        assertEquals(eval.getIntegralOrder(), 1);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.INTEGRAL_EVALUATION);        
    }
    
    @Test
    public void testGetSetX() {
        IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();
        
        //check default value
        assertEquals(eval.getX(), 0.0, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setX(x);
        
        //check correctness
        assertEquals(eval.getX(), x, 0.0);
    }
    
    @Test
    public void testGetSetConstants() {
        IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();
        
        //check default value
        assertNull(eval.getConstants());
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        double[] constants = new double[1];
        randomizer.fill(constants, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        eval.setConstant(constants);
        
        //check correctness
        assertSame(eval.getConstants(), constants);
    }
    
    @Test
    public void testGetSetIntegralOrder() {
        IntegralPolynomialEvaluation eval = new IntegralPolynomialEvaluation();
        
        //check default value
        assertEquals(eval.getIntegralOrder(), 1);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setIntegralOrder(order);
        
        //check correctness
        assertEquals(eval.getIntegralOrder(), order);
        
        //Force IllegalArgumentException
        try {
            eval.setIntegralOrder(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
}
