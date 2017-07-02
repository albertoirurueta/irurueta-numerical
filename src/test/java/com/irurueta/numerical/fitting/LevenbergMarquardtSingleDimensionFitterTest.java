/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.fitting.LevenbergMarquardtSingleDimensionFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 30, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.numerical.GradientEstimator;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
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

public class LevenbergMarquardtSingleDimensionFitterTest {
    public static final int MIN_POINTS = 500;
    public static final int MAX_POINTS = 1000;
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    
    public static final double MIN_SIGMA_VALUE = 1e-4;
    public static final double MAX_SIGMA_VALUE = 1e-3;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    
    public static final int GAUSS_PARAMS = 3;
    
    public static final int MIN_GAUSSIANS = 1;
    public static final int MAX_GAUSSIANS = 3;
    
    public static final double MIN_SINE_AMPLITUDE = 0.5;
    public static final double MAX_SINE_AMPLITUDE = 10.0;
    
    public static final double MIN_SINE_FREQ = 0.5;
    public static final double MAX_SINE_FREQ = 100.0;
    
    public static final double MIN_SINE_PHASE = -Math.PI;
    public static final double MAX_SINE_PHASE = Math.PI;
    
    public static final int SINE_PARAMS = 3;
    
    public LevenbergMarquardtSingleDimensionFitterTest() {}
    
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
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();
        
        //check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        
        //test constructor with input data
        double[] x = new double[2];
        double[] y = new double[2];
        double[] sig = new double[2];
        
        fitter = new LevenbergMarquardtSingleDimensionFitter(x, y, sig);
        
        //check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());        
        
        //Force IllegalArgumentException
        double[] shortX = new double[1];
        double[] shortSig = new double[1];
        
        fitter = null;
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(shortX, y, 
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(x, shortX, 
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(x, y, 
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtSingleDimensionFitter(x, y, 1.0);
        
        //check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());        
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(shortX, y, 
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(x, shortX, 
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with evaluator
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_PARAMS];
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                return 0.0;
            }
        };
        
        fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator);
        
        //check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_PARAMS);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_PARAMS);
        
        //test constructor with input data
        fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                sig);
        
        //check default value
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_PARAMS);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_PARAMS);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, 
                    shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, 
                    shortX, sig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, 
                    y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);
        
        //test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                1.0);
        
        //check default value
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for(int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_PARAMS);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_PARAMS);
        
        //Force IllegalArgumentException
        fitter = null;
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator,
                    shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        try{
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, 
                    shortX, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(fitter);        
    }
    
    @Test
    public void testGetSetNdone(){
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();
        
        //check default values
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtSingleDimensionFitter.NDONE);
        
        //set new value
        fitter.setNdone(5);
        
        //check correctness
        assertEquals(fitter.getNdone(), 5);
        
        //force IllegalArgumentException
        try{
            fitter.setNdone(0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testGetSetItmax(){
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();

        //check default values
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtSingleDimensionFitter.ITMAX);
        
        //set new value
        fitter.setItmax(10);
        
        //check correctness
        assertEquals(fitter.getItmax(), 10);
        
        //force IllegalArgumentException
        try{
            fitter.setItmax(0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testGetSetTol(){
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();

        //check default values
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtSingleDimensionFitter.TOL, 0.0);
        
        //set new value
        fitter.setTol(1e-1);
        
        //check correctness
        assertEquals(fitter.getTol(), 1e-1, 0.0);
        
        //force IllegalArgumentException
        try{
            fitter.setTol(0.0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
    }
    
    @Test
    public void testGetSetFunctionEvaluator() throws FittingException{
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();

        //check default value
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());        
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_PARAMS];
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                return 0.0;
            }
        };
        
        //set new value
        fitter.setFunctionEvaluator(evaluator);
        
        //check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_PARAMS);        
    }
    
    @Test
    public void testGetSetInputData(){
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();
        
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
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();
        
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
        for(int i = 0; i < fitter.getSig().length; i++) {
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
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter();
        
        //check default value
        assertFalse(fitter.isReady());
        
        //set new values
        double[] x = new double[2];
        double[] y = new double[2];
        double[] sig = new double[2];        

        fitter.setInputData(x, y, sig);
        
        assertFalse(fitter.isReady());
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_PARAMS];
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                return 0.0;
            }
        };

        fitter.setFunctionEvaluator(evaluator);
        
        assertTrue(fitter.isReady());
    }
    
    @Test
    public void testFitGaussian() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());        
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
        final int numParams = numgaussians * GAUSS_PARAMS;
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[numParams];
        for(int i = 0; i < numParams; i++){
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
                
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = 0.0;
            for(int k = 0; k < numgaussians; k++){
                double b = params[k * GAUSS_PARAMS];
                double e = params[k * GAUSS_PARAMS + 1];
                double g = params[k * GAUSS_PARAMS + 2];
                y[i] += b * Math.exp(-Math.pow((x[i] - e) / g, 2.0));
            }      
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[numParams];
                double error;
                for(int i = 0; i < numParams; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(int pos, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                int i, na = params.length;
                double fac, ex, arg;
                double y = 0.0;
                for(i = 0; i < na - 1; i += 3){
                    arg = (point - params[i + 1]) / params[i + 2];
                    ex = Math.exp(-Math.pow(arg, 2.0));
                    fac = params[i] * ex * 2. * arg;
                    y += params[i] * ex;
                    derivatives[i] = ex;
                    derivatives[i + 1] = fac / params[i + 2];
                    derivatives[i + 2] = fac * arg / params[i + 2];
                }
                
                return y;
            }
        };
        
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                1.0);
        
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
        assertEquals(fitter.getA().length, numParams);
        for(int i = 0; i < numParams; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
    @Test
    public void testFitSine() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());        
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;
                
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = amplitude * Math.sin(freq * x[i] + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[SINE_PARAMS];
                double error;
                for(int i = 0; i < SINE_PARAMS; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                double y = amplitude * Math.sin(freq * point + phase);
                
                //derivative respect amplitude
                derivatives[0] = Math.sin(freq * point + phase);
                //derivative respect frequency
                derivatives[1] = amplitude * Math.cos(freq * point + phase) * point;
                //derivative respect phase
                derivatives[2] = amplitude * Math.cos(freq * point + phase);
                
                return y;
            }
        };
        
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                1.0);
        
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
        assertEquals(fitter.getA().length, SINE_PARAMS);
        for(int i = 0; i < SINE_PARAMS; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
    @Test
    public void testFitSineWithHoldAndFree() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());        
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;
                
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = amplitude * Math.sin(freq * x[i] + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[SINE_PARAMS];
                initParams[0] = params[0];
                double error;
                for(int i = 1; i < SINE_PARAMS; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                double y = amplitude * Math.sin(freq * point + phase);
                
                //derivative respect amplitude
                derivatives[0] = Math.sin(freq * point + phase);
                //derivative respect frequency
                derivatives[1] = amplitude * Math.cos(freq * point + phase) * point;
                //derivative respect phase
                derivatives[2] = amplitude * Math.cos(freq * point + phase);
                
                return y;
            }
        };
        
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                1.0);
        
        //check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        
        //hold first parameter
        fitter.hold(0, params[0]);
        
        //fit
        fitter.fit();
        
        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, SINE_PARAMS);
        //first parameter is hold and matches exactly
        assertEquals(fitter.getA()[0], params[0], 0.0);
        for(int i = 0; i < SINE_PARAMS; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);   
        
        //release first parameter
        fitter.free(0);
        
        //fit and check correctness
        fitter.fit();
        
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, SINE_PARAMS);
        for(int i = 0; i < SINE_PARAMS; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);           
    }    
    
    @Test
    public void testFitGaussianWithGradientEstimator() throws FittingException, NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());        
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
        final int numParams = numgaussians * GAUSS_PARAMS;
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[numParams];
        for(int i = 0; i < numParams; i++){
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
                
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = 0.0;
            for(int k = 0; k < numgaussians; k++){
                double b = params[k * GAUSS_PARAMS];
                double e = params[k * GAUSS_PARAMS + 1];
                double g = params[k * GAUSS_PARAMS + 2];
                y[i] += b * Math.exp(-Math.pow((x[i] - e) / g, 2.0));
            }      
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
                    
            private double point;

            private GradientEstimator gradientEstimator = 
                    new GradientEstimator(
                            new MultiDimensionFunctionEvaluatorListener() {

                @Override
                public double evaluate(double[] params) throws Throwable {
                    return evaluateParams(point, params);
                }
            });
                
            public double evaluateParams(double point, double[] params){
                int i, na = params.length;
                double ex, arg;
                double y = 0.0;
                for(i = 0; i < na - 1; i += 3){
                    arg = (point - params[i + 1]) / params[i + 2];
                    ex = Math.exp(-Math.pow(arg, 2.0));
                    y += params[i] * ex;
                }

                return y;                    
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[numParams];
                double error;
                for(int i = 0; i < numParams; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                this.point = point;
                double y = evaluateParams(point, params);
                gradientEstimator.gradient(params, derivatives);
                
                return y;
            }
        };
        
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                1.0);
        
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
        assertEquals(fitter.getA().length, numParams);
        for(int i = 0; i < numParams; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
    @Test
    public void testFitSineWithGradientEstimator() throws FittingException, 
            NotReadyException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());        
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;
                
        double[] y = new double[npoints];
        double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for(int i = 0; i < npoints; i++){
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = amplitude * Math.sin(freq * x[i] + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator = 
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
                    
            private double point;
            
            private GradientEstimator gradientEstimator =
                    new GradientEstimator(
                            new MultiDimensionFunctionEvaluatorListener(){

                @Override
                public double evaluate(double[] params) throws Throwable {
                    return evaluateParams(point, params);
                }
            });
            
            public double evaluateParams(double point, double [] params){
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                return amplitude * Math.sin(freq * point + phase);                
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[SINE_PARAMS];
                double error;
                for(int i = 0; i < SINE_PARAMS; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(int i, double point, double[] params, 
                    double[] derivatives) throws Throwable {
                this.point = point;
                double y = evaluateParams(point, params);
                gradientEstimator.gradient(params, derivatives);
                
                return y;
            }
        };
        
        LevenbergMarquardtSingleDimensionFitter fitter = 
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y, 
                1.0);
        
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
        assertEquals(fitter.getA().length, SINE_PARAMS);
        for(int i = 0; i < SINE_PARAMS; i++){
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
}
