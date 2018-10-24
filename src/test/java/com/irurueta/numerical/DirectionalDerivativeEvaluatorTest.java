/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical;

import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

@SuppressWarnings("Duplicates")
public class DirectionalDerivativeEvaluatorTest {
    
    private static final int MIN_DIMS = 2;
    private static final int MAX_DIMS = 5;
    
    private static final double MIN_EVAL_POINT = -1e3;
    private static final double MAX_EVAL_POINT = 1e3;
    
    private static final double MIN_DIRECTION = -1.0;
    private static final double MAX_DIRECTION = 1.0;
    
    private static final double MIN_OFFSET = -1e3;
    private static final double MAX_OFFSET = 1e3;
    
    private static final double MIN_WIDTH = 1.0;
    private static final double MAX_WIDTH = 2.0;
    
    private int ndims;
    private double[] minimum;
    private double[] width;
    private double offset;
    
    private MultiDimensionFunctionEvaluatorListener listener; 
    private GradientFunctionEvaluatorListener gradientListener;
    
    public DirectionalDerivativeEvaluatorTest() {
        listener = new MultiDimensionFunctionEvaluatorListener() {

            @Override
            public double evaluate(double[] point) {
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
        
        gradientListener = new GradientFunctionEvaluatorListener() {

            @Override
            public void evaluateGradient(double[] params, double[] result) {
                int dims = Math.min(Math.min(params.length, minimum.length), 
                    width.length);        
        
                double value;
                for(int j = 0; j < dims; j++){
                    value = 1.0;
                    for(int i = 0; i < dims; i++){
                        if(i != j){
                            value *= Math.pow(params[i] - minimum[i], 2.0) / 
                                    width[i];
                        }else{
                            value *= 2.0 * (params[i] - minimum[i]) / width[i];
                        }
                    }
            
                    result[j] = value;
                }        
            }
        };
    }

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
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        
        minimum = new double[ndims];
        double[] point = new double[ndims];
        double[] direction = new double[ndims];
        width = new double[ndims];
        randomizer.fill(minimum, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);
        randomizer.fill(direction, MIN_DIRECTION, MAX_DIRECTION);
        offset = randomizer.nextDouble(MIN_OFFSET, MAX_OFFSET);
        randomizer.fill(width, MIN_WIDTH, MAX_WIDTH);
        
        DirectionalDerivativeEvaluator evaluator = 
                new DirectionalDerivativeEvaluator(listener, gradientListener, 
                point, direction);
        assertNotNull(evaluator);
    }
    
    @Test
    public void testEvaluateAt() throws Throwable {
        
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
        
        DirectionalDerivativeEvaluator evaluator = 
                new DirectionalDerivativeEvaluator(listener, gradientListener,
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
    
    @Test
    public void testDifferentiateAt() throws Throwable {
        
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
        
        DirectionalDerivativeEvaluator evaluator = 
                new DirectionalDerivativeEvaluator(listener, gradientListener,
                point, direction);
        
        double diff = evaluator.differentiateAt(x);
        
        //check correctness
        double[] xt = new double[ndims];
        for (int i = 0; i < ndims; i++) {
            xt[i] = point[i] + x * direction[i];
        }
        
        double[] gradient = new double[ndims];
        gradientListener.evaluateGradient(xt, gradient);
        
        double diff2 = 0.0;
        for (int i = 0; i < ndims; i++) {
            diff2 += gradient[i] * direction[i];
        }
        
        assertEquals(diff, diff2, 0.0);
    }
}
