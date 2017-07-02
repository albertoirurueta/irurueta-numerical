/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.PolynomialEvaluator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 9, 2016
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

public class PolynomialEvaluatorTest {
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    public static final int MIN_LENGTH = 1;
    public static final int MAX_LENGTH = 5;
    
    public static final double ABSOLUTE_ERROR = 1e-9;

    public static final int TIMES = 10;
    
    public PolynomialEvaluatorTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}
    
    @Test
    public void testEvaluateRealConstant(){
        
        double[] polyParams = new double[1];
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        assertEquals(PolynomialEvaluator.evaluate(polyParams, x), 
                polyParams[0], 0.0);
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((double[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEvaluateRealFirstDegree(){
        
        double[] polyParams = new double[2];
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0]*x + polyParams[1], ABSOLUTE_ERROR);
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((double[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }

    @Test
    public void testEvaluateRealSecondDegree(){
        
        double[] polyParams = new double[3];
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0]*x*x + polyParams[1]*x + polyParams[2], 
                ABSOLUTE_ERROR);
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((double[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testEvaluateReal(){
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
        int degree = length - 1;
        
        double[] polyParams = new double[length];
        randomizer.fill(polyParams, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        double x = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        double value = 0.0;
        for(int i = 0; i < length; i++){
            value += polyParams[i]*Math.pow(x, degree - i);
        }
        
        assertEquals(PolynomialEvaluator.evaluate(polyParams, x), value, 
                ABSOLUTE_ERROR);
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((double[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new double[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testEvaluateComplexConstant(){
        
        Complex[] polyParams = new Complex[1];
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for(int i = 0; i < polyParams.length; i++){
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }
        
        Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        
        assertEquals(PolynomialEvaluator.evaluate(polyParams, x), 
                polyParams[0]);
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((Complex[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new Complex[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEvaluateComplexFirstDegree(){
        
        Complex[] polyParams = new Complex[2];
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for(int i = 0; i < polyParams.length; i++){
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }
        
        Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        
        assertEquals(PolynomialEvaluator.evaluate(polyParams, x),
                polyParams[0].multiplyAndReturnNew(x).
                addAndReturnNew(polyParams[1]));
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((Complex[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new Complex[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }

    @Test
    public void testEvaluateComplexSecondDegree(){
        
        Complex[] polyParams = new Complex[3];
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        for(int i = 0; i < polyParams.length; i++){
            polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
        }
        
        Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));

        assertTrue(PolynomialEvaluator.evaluate(polyParams, x).equals(
                polyParams[0].multiplyAndReturnNew(x.multiplyAndReturnNew(x)).
                addAndReturnNew(polyParams[1].multiplyAndReturnNew(x)).
                addAndReturnNew(polyParams[2]), ABSOLUTE_ERROR));
        
        //Force IllegalArgumentException
        try{
            PolynomialEvaluator.evaluate((Complex[])null, x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            PolynomialEvaluator.evaluate(new Complex[0], x);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testEvaluateComplex(){
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
            int degree = length - 1;

            Complex[] polyParams = new Complex[length];
            for (int i = 0; i < polyParams.length; i++) {
                polyParams[i] = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));
            }

            Complex x = new Complex(randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE), randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));

            Complex value = new Complex();
            for (int i = 0; i < length; i++) {
                value.add(polyParams[i].multiplyAndReturnNew(
                        x.powAndReturnNew(degree - i)));
            }

            if (!PolynomialEvaluator.evaluate(polyParams, x).equals(value, ABSOLUTE_ERROR)) {
                continue;
            }
            assertTrue(PolynomialEvaluator.evaluate(polyParams, x).
                    equals(value, ABSOLUTE_ERROR));

            //Force IllegalArgumentException
            try {
                PolynomialEvaluator.evaluate((Complex[]) null, x);
                fail("IllegalArgumentException expected but not thrown");
            } catch (IllegalArgumentException e) { }
            try {
                PolynomialEvaluator.evaluate(new Complex[0], x);
                fail("IllegalArgumentException expected but not thrown");
            } catch (IllegalArgumentException e) { }

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }    
}
