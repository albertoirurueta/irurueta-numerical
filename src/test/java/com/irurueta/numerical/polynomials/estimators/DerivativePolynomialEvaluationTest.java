/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.DerivativePolynomialEvaluation
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

public class DerivativePolynomialEvaluationTest {
    
    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;
    public static final int MIN_ORDER = 1;
    public static final int MAX_ORDER = 4;
    
    public DerivativePolynomialEvaluationTest() { }
    
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
        DerivativePolynomialEvaluation eval = 
                new DerivativePolynomialEvaluation();
        
        //check initial values
        assertEquals(eval.getDerivativeOrder(), 1);
        assertEquals(eval.getX(), 0.0, 0.0);
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.DERIVATIVE_EVALUATION);
        
        //test constructor with values
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        
        eval = new DerivativePolynomialEvaluation(x, evaluation, order);
        
        //check correctness
        assertEquals(eval.getDerivativeOrder(), order);
        assertEquals(eval.getX(), x, 0.0);
        assertEquals(eval.getEvaluation(), evaluation, 0.0);        
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.DERIVATIVE_EVALUATION);

        //test constructor with values (without order)
        eval = new DerivativePolynomialEvaluation(x, evaluation);
        
        //check correctness
        assertEquals(eval.getDerivativeOrder(), 1);
        assertEquals(eval.getX(), x, 0.0);
        assertEquals(eval.getEvaluation(), evaluation, 0.0);        
        assertEquals(eval.getType(), 
                PolynomialEvaluationType.DERIVATIVE_EVALUATION);        
    }
    
    @Test
    public void testGetSetEvaluation() {
        DerivativePolynomialEvaluation eval = 
                new DerivativePolynomialEvaluation();

        //check default value
        assertEquals(eval.getEvaluation(), 0.0, 0.0);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        double evaluation = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                MAX_RANDOM_VALUE);
        eval.setEvaluation(evaluation);
        
        //check correctness
        assertEquals(eval.getEvaluation(), evaluation, 0.0);
    }
    
    @Test
    public void testGetSetX() {
        DerivativePolynomialEvaluation eval =
                new DerivativePolynomialEvaluation();
        
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
    public void testGetSetDerivativeOrder() {
        DerivativePolynomialEvaluation eval =
                new DerivativePolynomialEvaluation();
        
        //check default value
        assertEquals(eval.getDerivativeOrder(), 1);
        
        //set new value
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int order = randomizer.nextInt(MIN_ORDER, MAX_ORDER);
        eval.setDerivativeOrder(order);
        
        //check correctness
        assertEquals(eval.getDerivativeOrder(), order);
        
        //Force IllegalArgumentException
        try {
            eval.setDerivativeOrder(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
}
