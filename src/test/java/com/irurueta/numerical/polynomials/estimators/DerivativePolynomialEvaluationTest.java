/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

public class DerivativePolynomialEvaluationTest {
    
    private static final double MIN_RANDOM_VALUE = -10.0;
    private static final double MAX_RANDOM_VALUE = 10.0;
    private static final int MIN_ORDER = 1;
    private static final int MAX_ORDER = 4;
    
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
        } catch (IllegalArgumentException ignore) { }
    }
}
