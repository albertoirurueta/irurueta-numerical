/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.JacobianEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 31, 2015
 */
package com.irurueta.numerical;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class JacobianEstimatorTest 
        implements MultiVariateFunctionEvaluatorListener{
    
    public static final int MIN_DIMS = 2;
    public static final int MAX_DIMS = 3;
    
    public static final int MIN_VARS = 1;
    public static final int MAX_VARS = 4;
    
    public static final double MIN_EVAL_POINT = -10.0;
    public static final double MAX_EVAL_POINT = 10.0;
    
    public static final double MIN_OFFSET = 0.0;
    public static final double MAX_OFFSET = 1.0;
    
    public static final double MIN_WIDTH = 1.0;
    public static final double MAX_WIDTH = 10.0;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    
    public static final int TIMES = 100;
    
    private int ndims;
    private int nvars;
    private Matrix minimums;
    private Matrix widths;
    private double[] offsets;
    
    public JacobianEstimatorTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}

    @Override
    public void evaluate(double[] point, double[] result) throws Throwable {
        int dims = Math.min(Math.min(point.length, minimums.getColumns()),
                widths.getColumns());
        int vars = result.length;
        
        for(int j = 0; j < vars; j++){
            
            double value = 1.0;
            
            for(int i = 0; i < dims; i++){
                value *= Math.pow(point[i] - minimums.getElementAt(j, i), 2.0) / 
                        widths.getElementAt(j, i);
            }
            
            value += offsets[j];
            
            result[j] = value;
        }
    }

    @Override
    public int getNumberOfVariables() {
        return nvars;
    }
    
    public Matrix jacobian(double[] params){
        int dims = Math.min(Math.min(params.length, minimums.getColumns()),
                widths.getColumns());
        int vars = nvars;
        
        try{
            Matrix jacobian = new Matrix(vars, dims);
            
            for(int k = 0; k < vars; k++){
                
                double value;
                for(int j = 0; j < dims; j++){
                    value = 1.0;
                    for(int i = 0; i < dims; i++){
                        if(i != j){
                            value *= Math.pow(params[i] - minimums.getElementAt(k, i), 2.0) /
                                widths.getElementAt(k, i);
                        }else{
                            value *= 2.0 * (params[i] - minimums.getElementAt(k, i)) / 
                                    widths.getElementAt(k, i);
                        }
                    }
            
                    jacobian.setElementAt(k, j, value);
                }
            }
            
            return jacobian;
        }catch(WrongSizeException e){
            return null;
        }
    }
    
    @Test
    public void testConstructor() throws WrongSizeException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
        nvars = randomizer.nextInt(MIN_VARS, MAX_VARS);
        
        minimums = Matrix.createWithUniformRandomValues(nvars, ndims, 
                MIN_EVAL_POINT, MAX_EVAL_POINT, new Random());
        double[] point = new double[ndims];
        randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);  
        offsets = new double[nvars];
        randomizer.fill(offsets, MIN_OFFSET, MAX_OFFSET);        
        
        widths = Matrix.createWithUniformRandomValues(nvars, ndims, 
                MIN_WIDTH, MAX_WIDTH, new Random());
        
        JacobianEstimator estimator = new JacobianEstimator(this);
        assertNotNull(estimator);        
    }
    
    @Test
    public void testJacobian() throws EvaluationException, WrongSizeException{
        
        for(int t = 0; t < TIMES; t++){
            UniformRandomizer randomizer = new UniformRandomizer(new Random());
            
            ndims = randomizer.nextInt(MIN_DIMS, MAX_DIMS);
            nvars = randomizer.nextInt(MIN_VARS, MAX_VARS);
        
            minimums = Matrix.createWithUniformRandomValues(nvars, ndims, 
                    MIN_EVAL_POINT, MAX_EVAL_POINT, new Random());
            double[] point = new double[ndims];
            randomizer.fill(point, MIN_EVAL_POINT, MAX_EVAL_POINT);        
            offsets = new double[nvars];
            randomizer.fill(offsets, MIN_OFFSET, MAX_OFFSET);        
        
            widths = Matrix.createWithUniformRandomValues(nvars, ndims, 
                    MIN_WIDTH, MAX_WIDTH, new Random());
            
            JacobianEstimator estimator = new JacobianEstimator(this);
            
            Matrix jacobian1 = estimator.jacobian(point);
            Matrix jacobian2 = new Matrix(nvars, ndims);
            estimator.jacobian(point, jacobian2);
            
            //check correctness
            Matrix jacobian3 = jacobian(point);
            for(int j = 0; j < nvars; j++){
                for(int i = 0; i < ndims; i++){
                    assertEquals(jacobian1.getElementAt(j, i), 
                            jacobian3.getElementAt(j, i), ABSOLUTE_ERROR);
                    assertEquals(jacobian2.getElementAt(j, i), 
                            jacobian3.getElementAt(j, i), ABSOLUTE_ERROR);                    
                }
            }
        }
    }
}
