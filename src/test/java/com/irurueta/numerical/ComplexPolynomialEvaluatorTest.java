/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.ComplexPolynomialEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 9, 2016.
 */
package com.irurueta.numerical;

import com.irurueta.algebra.Complex;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class ComplexPolynomialEvaluatorTest {
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    public static final int MIN_LENGTH = 1;
    public static final int MAX_LENGTH = 5;    
    
    public ComplexPolynomialEvaluatorTest() {}
    
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
        
        Complex[] polyParams = new Complex[length];
        
        ComplexPolynomialEvaluator evaluator = 
                new ComplexPolynomialEvaluator(polyParams);
        
        //check correctness
        assertSame(evaluator.getPolyParams(), polyParams);        
    }
    
    @Test
    public void testEvaluate(){
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
        
        Complex[] polyParams = new Complex[length];
        for(int i = 0; i < polyParams.length; i++){
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }
        
        Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        
        ComplexPolynomialEvaluator evaluator = 
                new ComplexPolynomialEvaluator(polyParams);
        
        assertTrue(evaluator.evaluate(x).equals(
                PolynomialEvaluator.evaluate(polyParams, x)));
    }
}
