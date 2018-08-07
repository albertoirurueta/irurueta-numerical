/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.fitting.LevenbergMarquardtMultiVariateFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date June 1, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.JacobianEstimator;
import com.irurueta.numerical.MultiVariateFunctionEvaluatorListener;
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

public class LevenbergMarquardtMultiVariateFitterTest {
    
    public static final int MIN_POINTS = 500;
    public static final int MAX_POINTS = 1000;
    
    public static final double MIN_RANDOM_VALUE = -100.0;
    public static final double MAX_RANDOM_VALUE = 100.0;
    
    public static final double MIN_SIGMA_VALUE = 1e-4;
    public static final double MAX_SIGMA_VALUE = 1e-3;
    
    public static final double ABSOLUTE_ERROR = 1e-1;
    public static final double LARGE_ABSOLUTE_ERROR = 1.0;
    
    public static final int GAUSS_UNI_PARAMS = 3;
    public static final int GAUSS_MULTI_PARAMS = 5; //B, Ex, Ey, Gx, Gy
    
    public static final int MIN_GAUSSIANS = 1;
    public static final int MAX_GAUSSIANS = 3;
    
    public static final double MIN_SINE_AMPLITUDE = 0.5;
    public static final double MAX_SINE_AMPLITUDE = 10.0;
    
    public static final double MIN_SINE_FREQ = 0.5;
    public static final double MAX_SINE_FREQ = 100.0;
    
    public static final double MIN_SINE_PHASE = -Math.PI;
    public static final double MAX_SINE_PHASE = Math.PI;
    
    public static final int SINE_UNI_PARAMS = 3;
    public static final int SINE_MULTI_PARAMS = 5; //amplitude, freqx, freqy, phasex, phasey
    
    public static final int NUM_DIMENSIONS = 2;
    
    public static final int NUM_VARIABLES = 2;

    public static final int TIMES = 10;
    
    public LevenbergMarquardtMultiVariateFitterTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() throws FittingException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        
        //test constructor with input data
        Matrix x = new Matrix(nPoints, 2);
        Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        double[] sig = new double[nPoints];
        
        fitter = new LevenbergMarquardtMultiVariateFitter(x, y, sig);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        
        //Force IllegalArgumentException
        Matrix shortX = new Matrix(nPoints - 1, 2);
        Matrix shortY = new Matrix(nPoints - 1, NUM_VARIABLES);
        double[] shortSig = new double[nPoints - 1];
        
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(fitter);
        
        //test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtMultiVariateFitter(x, y, 1.0);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        
        //Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(fitter);
        
        //test constructor with evaluator
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public int getNumberOfVariables() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[nPoints];
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
            }
        };
        
        fitter = new LevenbergMarquardtMultiVariateFitter(evaluator);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, nPoints);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), nPoints);
        assertEquals(fitter.getCovar().getColumns(), nPoints);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        
        //test constructor with evaluator and input data
        fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, sig);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, nPoints);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), nPoints);
        assertEquals(fitter.getCovar().getColumns(), nPoints);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        
        //Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, shortX,
                    y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x,
                    shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(fitter);
        
        //test constructor with evaluator and input data (constant sigma)
        fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
        //check default values
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, nPoints);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), nPoints);
        assertEquals(fitter.getCovar().getColumns(), nPoints);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        
        //Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator,
                    shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator,
                    x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        assertNull(fitter);
    }

    @Test
    public void testGetSetNdone() {
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default value
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtMultiVariateFitter.NDONE);
        
        //new value
        fitter.setNdone(5);
        
        //check correctness
        assertEquals(fitter.getNdone(), 5);
        
        //force IllegalArgumentException
        try {
            fitter.setNdone(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetItmax() {
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default value
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiVariateFitter.ITMAX);
        
        //new value
        fitter.setItmax(10);
        
        //check correctness
        assertEquals(fitter.getItmax(), 10);
        
        //force IllegalArgumentException
        try {
            fitter.setItmax(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetTol() {
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default value
        assertEquals(fitter.getTol(),
                LevenbergMarquardtMultiVariateFitter.TOL, 0.0);
        
        //new value
        fitter.setTol(1e-1);
        
        //check correctness
        assertEquals(fitter.getTol(), 1e-1, 0.0);
        
        //force IllegalArgumentException
        try {
            fitter.setTol(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetFunctionEvaluator() throws FittingException {
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default value
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertNull(fitter.getAlpha());
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public int getNumberOfVariables() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_UNI_PARAMS];
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
            }
        };
        
        //new value
        fitter.setFunctionEvaluator(evaluator);
        
        //check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_UNI_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_UNI_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_UNI_PARAMS);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_UNI_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_UNI_PARAMS);
    }
    
    @Test
    public void testGetSetInputData() throws WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        Matrix x = new Matrix(nPoints, 2);
        Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        double[] sig = new double[nPoints];
        
        //set input data
        fitter.setInputData(x, y, sig);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        
        //Force IllegalArgumentException
        Matrix shortX = new Matrix(nPoints - 1, 2);
        Matrix shortY = new Matrix(nPoints - 1, NUM_VARIABLES);
        double[] shortSig = new double[nPoints - 1];
        
        try {
            fitter.setInputData(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter.setInputData(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter.setInputData(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testGetSetInputDataWithConstantSigma() 
            throws WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        Matrix x = new Matrix(nPoints, 2);
        Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        
        //set input data
        fitter.setInputData(x, y, 1.0);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        
        
        //Force IllegalArgumentException
        Matrix shortX = new Matrix(nPoints - 1, 2);
        Matrix shortY = new Matrix(nPoints - 1, NUM_VARIABLES);
        
        try {
            fitter.setInputData(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
        try {
            fitter.setInputData(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException e) { }
    }
    
    @Test
    public void testIsReady() throws WrongSizeException, FittingException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();
        
        //check default value
        assertFalse(fitter.isReady());
        
        //set new values
        Matrix x = new Matrix(nPoints, 2);
        Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        double[] sig = new double[nPoints];
        
        fitter.setInputData(x, y, sig);
        
        assertFalse(fitter.isReady());
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public int getNumberOfVariables() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_UNI_PARAMS];
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
            }
        };

        fitter.setFunctionEvaluator(evaluator);
        
        assertTrue(fitter.isReady());
    }
    
    @Test
    public void testFitUnidimensionalGaussian() throws FittingException,
            NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
        final int numParams = numgaussians * GAUSS_UNI_PARAMS;
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
        
        Matrix y = new Matrix(npoints, 1);
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            double value = 0.0;
            for(int k = 0; k < numgaussians; k++){
                double b = params[k * GAUSS_UNI_PARAMS];
                double e = params[k * GAUSS_UNI_PARAMS + 1];
                double g = params[k * GAUSS_UNI_PARAMS + 2];
                value += b * Math.exp(-Math.pow((x.getElementAt(i, 0) - e) / g, 2.0));
            }
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public int getNumberOfVariables() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams = new double[numParams];
                double error;
                for (int i = 0; i < numParams; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int pos, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                int i, na = params.length;
                double fac, ex, arg;
                double y = 0.0;
                for (i = 0; i < na - 1; i += 3) {
                    arg = (point[0] - params[i + 1]) / params[i + 2];
                    ex = Math.exp(-Math.pow(arg, 2.0));
                    fac = params[i] * ex * 2. * arg;
                    y += params[i] * ex;
                    jacobian.setElementAt(0, i, ex);
                    jacobian.setElementAt(0, i + 1, fac / params[i + 2]);
                    jacobian.setElementAt(0, i + 2, fac * arg / params[i + 2]);
                }

                result[0] = y;
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        for (int i = 0; i < numParams; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);                        
    }
    
    @Test
    public void testFitUnidimensionalSine() throws FittingException,
            NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;
        
        Matrix y = new Matrix(npoints, 1);
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            double value = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public int getNumberOfVariables() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams = new double[SINE_UNI_PARAMS];
                double error;
                for(int i = 0; i < SINE_UNI_PARAMS; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                double y = amplitude * Math.sin(freq * point[0] + phase);
                
                //derivative respect amplitude
                jacobian.setElementAt(0, 0, Math.sin(freq * point[0] + phase));
                jacobian.setElementAt(0, 1, amplitude * Math.cos(
                        freq * point[0] + phase) * point[0]);
                jacobian.setElementAt(0, 2, amplitude * Math.cos(
                        freq * point[0] + phase));

                result[0] = y;
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
    @Test
    public void testFitUnidimensionalSineWithHoldAndFree()
            throws FittingException, NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;
        
        Matrix y = new Matrix(npoints, 1);
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            double value = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public int getNumberOfVariables() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams = new double[SINE_UNI_PARAMS];
                double error;
                for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                double y = amplitude * Math.sin(freq * point[0] + phase);
                
                //derivative respect amplitude
                jacobian.setElementAt(0, 0, Math.sin(freq * point[0] + phase));
                jacobian.setElementAt(0, 1, amplitude * Math.cos(
                        freq * point[0] + phase) * point[0]);
                jacobian.setElementAt(0, 2, amplitude * Math.cos(
                        freq * point[0] + phase));

                result[0] = y;
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
        //first parameter is hold and matches exactly
        assertEquals(fitter.getA()[0], params[0], 0.0);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
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
        assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);                           
    }
    
    @Test
    public void testFitUnidimensionalGaussianWithJacobianEstimator()
            throws FittingException, NotReadyException, WrongSizeException {

        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_UNI_PARAMS;

            double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            Matrix y = new Matrix(npoints, 1);
            Matrix x = new Matrix(npoints, 1);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));
                double value = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    double b = params[k * GAUSS_UNI_PARAMS];
                    double e = params[k * GAUSS_UNI_PARAMS + 1];
                    double g = params[k * GAUSS_UNI_PARAMS + 2];
                    value += b * Math.exp(-Math.pow((x.getElementAt(i, 0) - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(double[] params, double[] result)
                                                    throws Throwable {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 1;
                                            }
                                        });

                        public void evaluateParams(double[] point, double[] params,
                                                   double[] result) {
                            int i, na = params.length;
                            double ex, arg;
                            double y = 0.0;
                            for (i = 0; i < na - 1; i += 3) {
                                arg = (point[0] - params[i + 1]) / params[i + 2];
                                ex = Math.exp(-Math.pow(arg, 2.0));
                                y += params[i] * ex;
                            }

                            result[0] = y;
                        }

                        @Override
                        public int getNumberOfDimensions() {
                            return 1;
                        }

                        @Override
                        public int getNumberOfVariables() {
                            return 1;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(int i, double[] point, double[] result,
                                             double[] params, Matrix jacobian) throws Throwable {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }
                    };

            LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

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
            boolean allValid = true;
            for (int i = 0; i < numParams; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    allValid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], 10 * ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            if (allValid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }
    
    @Test
    public void testFitUnidimensionalSineWithJacobianEstimator() 
            throws FittingException, NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;
        
        Matrix y = new Matrix(npoints, 1);
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            double value = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {
                    
            private double[] point;
            
            private JacobianEstimator jacobianEstimator =
                    new JacobianEstimator(
                            new MultiVariateFunctionEvaluatorListener() {

                @Override
                public void evaluate(double[] params, double[] result) 
                        throws Throwable {
                    evaluateParams(point, params, result);
                }

                @Override
                public int getNumberOfVariables() {
                    return 1;
                }
            });
            
            public void evaluateParams(double[] point, double[] params, 
                    double[] result) {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                result[0] = amplitude * Math.sin(freq * point[0] + phase);
            }

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public int getNumberOfVariables() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams = new double[SINE_UNI_PARAMS];
                double error;
                for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                this.point = point;
                evaluateParams(point, params, result);
                jacobianEstimator.jacobian(params, jacobian);
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
    @Test
    public void testFitMultidimensionalGaussian() throws FittingException,
            NotReadyException, WrongSizeException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS,
                    MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_MULTI_PARAMS;

            double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            Matrix y = new Matrix(npoints, 1);
            Matrix x = new Matrix(npoints, 2);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                double value = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    double b = params[k * GAUSS_MULTI_PARAMS];
                    double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    value += b * Math.exp(-(Math.pow(x.getElementAt(i, 0) - ex, 2.0) +
                            Math.pow(x.getElementAt(i, 1) - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                }
                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(double[] params, double[] result)
                                                    throws Throwable {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 1;
                                            }
                                        });

                        public void evaluateParams(double[] point, double[] params,
                                                   double[] result) {
                            double y = 0.0;
                            for (int k = 0; k < numgaussians; k++) {
                                double b = params[k * GAUSS_MULTI_PARAMS];
                                double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                                double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                                double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                                double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                                y += b * Math.exp(-(Math.pow(point[0] - ex, 2.0) +
                                        Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                            }

                            result[0] = y;
                        }

                        @Override
                        public int getNumberOfDimensions() {
                            return NUM_DIMENSIONS;
                        }

                        @Override
                        public int getNumberOfVariables() {
                            return 1;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(int i, double[] point, double[] result,
                                             double[] params, Matrix jacobian) throws Throwable {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }
                    };

            LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

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
            boolean failed = false;
            for (int i = 0; i < numParams; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > LARGE_ABSOLUTE_ERROR) {
                    failed = true;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], LARGE_ABSOLUTE_ERROR);
            }

            if (failed) {
                continue;
            }
            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);

        assertTrue(numValid > 0);
    }
    
    @Test
    public void testFitMultidimensionalSine() throws FittingException,
            NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double freqy = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        double phasey = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_MULTI_PARAMS];
        params[0] = amplitude;
        params[1] = freqx;
        params[2] = freqy;
        params[3] = phasex;
        params[4] = phasey;
        
        Matrix y = new Matrix(npoints, 1);
        Matrix x = new Matrix(npoints, 2);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            
            double value = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                    Math.sin(freqy * x.getElementAt(i, 1) + phasey);
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

            private double[] point;
            
            private JacobianEstimator jacobianEstimator =
                    new JacobianEstimator(
                            new MultiVariateFunctionEvaluatorListener() {

                @Override
                public void evaluate(double[] params, double[] result) 
                        throws Throwable {
                    evaluateParams(point, params, result);
                }

                @Override
                public int getNumberOfVariables() {
                    return 1;
                }
            });
            
            public void evaluateParams(double[] point, double[] params, 
                    double[] result) {
                double amplitude = params[0];
                double freqx = params[1];
                double freqy = params[2];
                double phasex = params[3];
                double phasey = params[4];
                
                double value = amplitude * Math.sin(freqx * point[0] + phasex) *
                    Math.sin(freqy * point[1] + phasey);
                
                result[0] = value;
            }
                    
            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public int getNumberOfVariables() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[SINE_MULTI_PARAMS];
                double error;
                for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                this.point = point;
                evaluateParams(point, params, result);
                jacobianEstimator.jacobian(params, jacobian);
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        assertEquals(fitter.getA().length, SINE_MULTI_PARAMS);
        for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);                
    }
    
    @Test
    public void testFitMultidimensionalSineRepeatInOneDimension()
            throws FittingException, NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, 
                MAX_SINE_AMPLITUDE);
        double freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        double phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freqx;
        params[2] = phasex;
        
        Matrix y = new Matrix(npoints, 1);
        Matrix x = new Matrix(npoints, 2);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            
            double value = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex);
            error = errorRandomizer.nextDouble();
            value += error;
            
            y.setElementAt(i, 0, value);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {
                    
            private double[] point;
            
            private JacobianEstimator jacobianEstimator =
                    new JacobianEstimator(
                            new MultiVariateFunctionEvaluatorListener() {

                @Override
                public void evaluate(double[] params, double[] result) 
                        throws Throwable {
                    evaluateParams(point, params, result);
                }

                @Override
                public int getNumberOfVariables() {
                    return 1;
                }
            });
            
            public void evaluateParams(double[] point, double[] params, 
                    double[] result) {
                double amplitude = params[0];
                double freqx = params[1];
                double phasex = params[2];
                
                result[0] = amplitude * Math.sin(freqx * point[0] + phasex);                
            }

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public int getNumberOfVariables() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[SINE_UNI_PARAMS];
                double error;
                for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                this.point = point;
                evaluateParams(point, params, result);
                jacobianEstimator.jacobian(params, jacobian);
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
    
    @Test
    public void testFitMultiVariateGaussianAndSine() throws FittingException,
            NotReadyException, WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        
        int npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
        final int numParams = SINE_MULTI_PARAMS;
        
        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
        
        final double[] params = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE);
        }
        
        Matrix y = new Matrix(npoints, 2);
        Matrix x = new Matrix(npoints, 2);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            double valueGaussian, valueSine;
            
            double b = params[0];
            double ex = params[1];
            double ey = params[2];
            double gx = params[3];
            double gy = params[4];
            error = errorRandomizer.nextDouble();
            
            valueGaussian = b * Math.exp(-(Math.pow(x.getElementAt(i, 0) - ex, 2.0) + 
                        Math.pow(x.getElementAt(i, 1) - ey, 2.0)) / Math.pow(gx * gy, 2.0));
            valueGaussian += error;
            
            double amplitude = params[0];
            double freqx = params[1];
            double freqy = params[2];
            double phasex = params[3];
            double phasey = params[4];
            
            valueSine = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                    Math.sin(freqy * x.getElementAt(i, 1) + phasey);
            valueSine += error;
            
            y.setElementAt(i, 0, valueGaussian);
            y.setElementAt(i, 1, valueSine);
        }
        
        LevenbergMarquardtMultiVariateFunctionEvaluator evaluator = 
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {
                    
            private double[] point;
            
            private JacobianEstimator jacobianEstimator =
                    new JacobianEstimator(
                            new MultiVariateFunctionEvaluatorListener() {

                @Override
                public void evaluate(double[] params, double[] result) 
                        throws Throwable {
                    evaluateParams(point, params, result);
                }

                @Override
                public int getNumberOfVariables() {
                    return 2;
                }
            });
            
            public void evaluateParams(double[] point, double[] params, 
                    double[] result) {
                double b = params[0];
                double ex = params[1];
                double ey = params[2];
                double gx = params[3];
                double gy = params[4];
                
                double amplitude = params[0];
                double freqx = params[1];
                double freqy = params[2];
                double phasex = params[3];
                double phasey = params[4];
            
                double valueGaussian = b * Math.exp(-(Math.pow(point[0] - ex, 2.0) + 
                        Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));
            
                double valueSine = amplitude * Math.sin(freqx * point[0] + phasex) *
                    Math.sin(freqy * point[1] + phasey);
            
                result[0] = valueGaussian;
                result[1] = valueSine;
            }

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public int getNumberOfVariables() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[numParams];
                double error;
                for (int i = 0; i < numParams; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public void evaluate(int i, double[] point, double[] result, 
                    double[] params, Matrix jacobian) throws Throwable {
                this.point = point;
                evaluateParams(point, params, result);
                jacobianEstimator.jacobian(params, jacobian);
            }
        };
        
        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        
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
        assertEquals(fitter.getA().length, SINE_MULTI_PARAMS);
        for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);        
    }
}
