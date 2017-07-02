/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.fitting.SvdSingleDimensionLinearFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 26, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class SvdSingleDimensionLinearFitterTest {
    public static final int MIN_POLY_PARAMS = 1;
    public static final int MAX_POLY_PARAMS = 5;
    
    public static final int MIN_POINTS = 500;
    public static final int MAX_POINTS = 1000;
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    
    public static final double MIN_SIGMA_VALUE = 1e-4;
    public static final double MAX_SIGMA_VALUE = 1.0;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    
    public static final int TRIGO_PARAMS = 2;

    
    public SvdSingleDimensionLinearFitterTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}

    @Test
    public void testConstructor() throws FittingException{
        //test empty constructor
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter();
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //test constructor with input data
        double[] x = new double[2];
        double[] y = new double[2];
        double[] sig = new double[2];
        
        fitter = new SvdSingleDimensionLinearFitter(x, y, sig);
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        double[] shortX = new double[1];
        double[] shortSig = new double[1];
        
        fitter = null;
        try{
            fitter = new SvdSingleDimensionLinearFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdSingleDimensionLinearFitter(x, shortX, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdSingleDimensionLinearFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with input data (constant sigma)
        fitter = new SvdSingleDimensionLinearFitter(x, y, 1.0);
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++){
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new SvdSingleDimensionLinearFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdSingleDimensionLinearFitter(x, shortX, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with evaluator
        LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[2];
            }

            @Override
            public void evaluate(double point, double[] result) 
                    throws Throwable {
            }
        };
        
        fitter = new SvdSingleDimensionLinearFitter(evaluator);
        
        //check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //test constructor with evaluator and input data
        fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, sig);
        
        //check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new SvdSingleDimensionLinearFitter(evaluator, shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdSingleDimensionLinearFitter(evaluator, x, shortX, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with evaluator and input data (constant sigma)
        fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, 1.0);
        
        //check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertTrue(fitter.isReady());
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new SvdSingleDimensionLinearFitter(evaluator, shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new SvdSingleDimensionLinearFitter(evaluator, x, shortX, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
    }
    
    @Test
    public void testGetSetFunctionEvaluator() throws FittingException{
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter();
        
        //check default values
        assertNull(fitter.getFunctionEvaluator());
        
        //set new value
        LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[2];
            }

            @Override
            public void evaluate(double point, double[] result) 
                    throws Throwable {}
        };
        
        //set new value
        fitter.setFunctionEvaluator(evaluator);
        
        //check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
    }
    
    @Test
    public void testGetSetInputData(){
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter();
        
        //check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        double[] x = new double[2];
        double[] y = new double[2];
        double[] sig = new double[2];
        
        //set input data
        fitter.setInputData(x, y, sig);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        
        //Force IllegalArgumentException
        double[] wrong = new double[1];
        
        try{
            fitter.setInputData(wrong, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputData(x, wrong, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputData(x, y, wrong);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testGetSetInputDataWithConstantSigma(){
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter();
        
        //check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        double[] x = new double[2];
        double[] y = new double[2];
        
        //set input data
        fitter.setInputData(x, y, 1.0);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++){
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        
        //Force IllegalArgumentException
        double[] wrong = new double[1];
        
        try{
            fitter.setInputData(wrong, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter.setInputData(x, wrong, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testIsReady() throws FittingException{
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter();
        
        //check default value
        assertFalse(fitter.isReady());
        
        //set new values
        double[] x = new double[2];
        double[] y = new double[2];
        double[] sig = new double[2];
        
        fitter.setInputData(x, y, sig);
        
        assertFalse(fitter.isReady());
        
        LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator(){

            @Override
            public double[] createResultArray() {
                return new double[2];
            }

            @Override
            public void evaluate(double point, double[] result) 
                    throws Throwable {}
        };
        
        fitter.setFunctionEvaluator(evaluator);
        
        assertTrue(fitter.isReady());
    }
    
    @Test
    public void testGetSetTol(){
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter();
        
        //check default values
        assertEquals(fitter.getTol(), SvdSingleDimensionLinearFitter.TOL, 0.0);
        
        //set new value
        fitter.setTol(1e-3);
        
        //check correctness
        assertEquals(fitter.getTol(), 1e-3, 0.0);
    }
    
    @Test
    public void testFitPolynomial() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        final int nparams = randomizer.nextInt(MIN_POLY_PARAMS, 
                MAX_POLY_PARAMS);
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        double[] params = new double[nparams];
        for(int i = 0; i < nparams; i++){
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
        
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = 0.0;
            for(int j = 0; j < nparams; j++){
                y[i] += params[j] * Math.pow(x[i], j);
            }
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[nparams];
            }

            @Override
            public void evaluate(double point, double[] result) 
                    throws Throwable {
                for(int i = 0; i < result.length; i++){
                    result[i] = Math.pow(point, i);
                }
            }
        };

        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, 1.0);
        
        //check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        //fit
        fitter.fit();
        
        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, nparams);
        for(int i = 0; i < nparams; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);          
    }
    
    @Test
    public void testTrigo() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        double[] params = new double[TRIGO_PARAMS];
        for(int i = 0; i < TRIGO_PARAMS; i++){
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
        
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        //function: y(x) = a * cos(x) + b * sin(x)
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            error = errorRandomizer.nextDouble();
            y[i] = params[0] * Math.sin(x[i]) + params[1] * Math.cos(x[i]) + 
                    error;
        }
        
        LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[TRIGO_PARAMS];
            }

            @Override
            public void evaluate(double point, double[] result) throws Throwable {
                result[0] = Math.sin(point);
                result[1] = Math.cos(point);
            }
        };
        
        SvdSingleDimensionLinearFitter fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, 1.0);
        
        //check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        
        //fit
        fitter.fit();
        
        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, TRIGO_PARAMS);
        for(int i = 0; i < TRIGO_PARAMS; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);
    }    
}
