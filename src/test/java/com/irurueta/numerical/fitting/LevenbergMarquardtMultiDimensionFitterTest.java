/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.Utils;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.GradientEstimator;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.*;
import org.junit.*;

import java.util.Arrays;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.Assert.*;

@SuppressWarnings("Duplicates")
public class LevenbergMarquardtMultiDimensionFitterTest {

    private static final Logger LOGGER = Logger.getLogger(
            LevenbergMarquardtMultiDimensionFitterTest.class.getName());
    
    private static final int MIN_POINTS = 500;
    private static final int MAX_POINTS = 1000;
    
    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;
    
    private static final double MIN_SIGMA_VALUE = 1e-4;
    private static final double MAX_SIGMA_VALUE = 1e-3;
    
    private static final double ABSOLUTE_ERROR = 1e-1;
    private static final double LARGE_ABSOLUTE_ERROR = 1.0;
    private static final double SMALL_ABSOLUTE_ERROR = 1e-5;

    private static final int CONSTANT_PARAMS = 1;
    private static final double MIN_CONSTANT = -100.0;
    private static final double MAX_CONSTANT = 100.0;

    private static final int LINE1_PARAMS = 1;
    private static final double MIN_LINE1_A = -3.0;
    private static final double MAX_LINE1_A = 3.0;

    private static final int LINE2_PARAMS = 2;
    private static final double MIN_LINE2_A = -3.0;
    private static final double MAX_LINE2_A = 3.0;
    private static final double MIN_LINE2_B = -10.0;
    private static final double MAX_LINE2_B = 10.0;

    private static final int SINE_UNI_PARAMS = 3;
    private static final int SINE_MULTI_PARAMS = 5; //amplitude, freqx, freqy, phasex, phasey
    private static final double MIN_SINE_AMPLITUDE = 0.5;
    private static final double MAX_SINE_AMPLITUDE = 10.0;
    private static final double MIN_SINE_FREQ = 0.5;
    private static final double MAX_SINE_FREQ = 100.0;
    private static final double MIN_SINE_PHASE = -Math.PI;
    private static final double MAX_SINE_PHASE = Math.PI;

    private static final int GAUSS_UNI_PARAMS = 3;
    private static final int GAUSS_MULTI_PARAMS = 5; //B, Ex, Ey, Gx, Gy
    private static final int MIN_GAUSSIANS = 1;
    private static final int MAX_GAUSSIANS = 3;

    private static final int NUM_DIMENSIONS = 2;

    public static final int TIMES = 10;
    private static final int N_SAMPLES = 1000000;
    
    public LevenbergMarquardtMultiDimensionFitterTest() { }
    
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
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
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
        assertTrue(fitter.isCovarianceAdjusted());
        
        //test constructor with input data
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        
        fitter = new LevenbergMarquardtMultiDimensionFitter(x, y, sig);
        
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
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        //Force IllegalArgumentException
        Matrix shortX = new Matrix(nPoints - 1, 2);
        double[] shortY = new double[nPoints - 1];
        double[] shortSig = new double[nPoints - 1];
        
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(fitter);
        
        //test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtMultiDimensionFitter(x, y, 1.0);
        
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
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        //Force IllegalArgumentException        
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(fitter);
        
        //test constructor with evaluator
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[nPoints];
            }

            @Override
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) {
                return 0.0;
            }
        };
        
        fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator);
        
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
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());
        
        //test constructor with evaluator and input data
        fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 
                sig);
        
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
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        //Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, 
                    shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, 
                    shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(fitter);
        
        //test constructor with evaluator and input data (constant sigma)
        fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 
                1.0);
        
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
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        //Force IllegalArgumentException        
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, 
                    shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, 
                    shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(fitter);        
    }
    
    @Test
    public void testGetSetNdone() {
        LevenbergMarquardtMultiDimensionFitter fitter = 
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default value
        assertEquals(fitter.getNdone(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE);
        
        //new value
        fitter.setNdone(5);
        
        //check correctness
        assertEquals(fitter.getNdone(), 5);
        
        //force IllegalArgumentException
        try {
            fitter.setNdone(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testGetSetItmax() {
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default value
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX);
        
        //new value
        fitter.setItmax(10);
        
        //check correctness
        assertEquals(fitter.getItmax(), 10);
        
        //force IllegalArgumentException
        try {
            fitter.setItmax(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testGetSetTol() {
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default value
        assertEquals(fitter.getTol(), 
                LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, 0.0);
        
        //new value
        fitter.setTol(1e-1);
        
        //check correctness
        assertEquals(fitter.getTol(), 1e-1, 0.0);
        
        //force IllegalArgumentException
        try {
            fitter.setTol(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testGetSetFunctionEvaluator() throws FittingException {
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default value
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertNull(fitter.getAlpha());
        
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_UNI_PARAMS];
            }

            @Override
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) {
                return 0.0;
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
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        
        //set input data
        fitter.setInputData(x, y, sig);
        
        //check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        
        //Force IllegalArgumentException
        Matrix shortX = new Matrix(nPoints - 1, 2);
        double[] shortY = new double[nPoints - 1];
        double[] shortSig = new double[nPoints - 1];
        
        try {
            fitter.setInputData(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter.setInputData(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter.setInputData(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testGetSetInputDataWithConstantSigma() 
            throws WrongSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        
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
        double[] shortY = new double[nPoints - 1];
        
        try {
            fitter.setInputData(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        try {
            fitter.setInputData(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testIsReady() throws WrongSizeException, FittingException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();
        
        //check default value
        assertFalse(fitter.isReady());
        
        //set new values
        Matrix x = new Matrix(nPoints, 2);
        double[] y = new double[nPoints];
        double[] sig = new double[nPoints];
        
        fitter.setInputData(x, y, sig);
        
        assertFalse(fitter.isReady());
        
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_UNI_PARAMS];
            }

            @Override
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) {
                return 0.0;
            }
        };
        
        fitter.setFunctionEvaluator(evaluator);
        
        assertTrue(fitter.isReady());        
    }

    @Test
    public void testIsSetCovarianceAdjusted() {
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter();

        //check default value
        assertTrue(fitter.isCovarianceAdjusted());

        //set new value
        fitter.setCovarianceAdjusted(false);

        //check
        assertFalse(fitter.isCovarianceAdjusted());
    }

    @Test
    public void testFitUnidimensionalConstant() throws WrongSizeException,
            FittingException, NotReadyException, MaxIterationsExceededException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());

        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[CONSTANT_PARAMS];
        params[0] = constant;

        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0,
                    randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            y[i] = constant;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                    @Override
                    public int getNumberOfDimensions() {
                        return 1;
                    }

                    @Override
                    public double[] createInitialParametersArray() {
                        double[] initParams = new double[CONSTANT_PARAMS];
                        double error;
                        for (int i = 0; i < CONSTANT_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public double evaluate(int i, double[] point, double[] params,
                                           double[] derivatives) {
                        double constant = params[0];

                        //derivative of evaluated function respect constant parameter
                        derivatives[0] = 1.0;

                        //evaluated function
                        return constant;
                    }
                };

        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

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
        assertEquals(fitter.getA().length, CONSTANT_PARAMS);
        for (int i = 0; i < CONSTANT_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        double chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalLine1() throws WrongSizeException,
            FittingException, NotReadyException, MaxIterationsExceededException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());

        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[LINE1_PARAMS];
        params[0] = a;

        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            y[i] = a * xi;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                    @Override
                    public int getNumberOfDimensions() {
                        return 1;
                    }

                    @Override
                    public double[] createInitialParametersArray() {
                        double[] initParams = new double[LINE1_PARAMS];
                        double error;
                        for (int i = 0; i < LINE1_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public double evaluate(int i, double[] point, double[] params,
                                           double[] derivatives) {
                        double a = params[0];

                        //derivative of evaluated function respect constant parameter
                        derivatives[0] = a;

                        //evaluated function
                        return a * point[0];
                    }
                };

        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

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
        assertEquals(fitter.getA().length, LINE1_PARAMS);
        for (int i = 0; i < LINE1_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        double chiSqrDegreesOfFreedom = npoints - LINE1_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalLine2() throws WrongSizeException,
            FittingException, NotReadyException, MaxIterationsExceededException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());

        int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        double a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
        double b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

        double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[LINE2_PARAMS];
        params[0] = a;
        params[1] = b;

        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            y[i] = a * xi + b;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                    @Override
                    public int getNumberOfDimensions() {
                        return 1;
                    }

                    @Override
                    public double[] createInitialParametersArray() {
                        double[] initParams = new double[LINE2_PARAMS];
                        double error;
                        for (int i = 0; i < LINE2_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public double evaluate(int i, double[] point, double[] params,
                                           double[] derivatives) {
                        double a = params[0];
                        double b = params[1];

                        //derivatives of function f(x) respect parameters a and b
                        derivatives[0] = point[0];
                        derivatives[1] = 1.0;

                        //evaluated function f(x) = a * x + b
                        return a * point[0] + b;
                    }
                };

        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

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
        assertEquals(fitter.getA().length, LINE2_PARAMS);
        for (int i = 0; i < LINE2_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        double chiSqrDegreesOfFreedom = npoints - LINE2_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalSine() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
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
        
        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            y[i] = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
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
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                double y = amplitude * Math.sin(freq * point[0] + phase);
                
                //derivative respect amplitude
                derivatives[0] = Math.sin(freq * point[0] + phase);
                //derivative respect frequency
                derivatives[1] = amplitude * Math.cos(freq * point[0] + phase) * point[0];
                //derivative respect phase
                derivatives[2] = amplitude * Math.cos(freq * point[0] + phase);
                
                return y;
            }
        };
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 
                1.0);
        fitter.setCovarianceAdjusted(false);
        
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

        double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalGaussian() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
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

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 1);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);
                y[i] = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    double b = params[k * GAUSS_UNI_PARAMS];
                    double e = params[k * GAUSS_UNI_PARAMS + 1];
                    double g = params[k * GAUSS_UNI_PARAMS + 2];
                    y[i] += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        @Override
                        public int getNumberOfDimensions() {
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
                        public double evaluate(int pos, double[] point, double[] params,
                                               double[] derivatives) {
                            int i, na = params.length;
                            double fac, ex, arg;
                            double y = 0.0;
                            for (i = 0; i < na - 1; i += 3) {
                                arg = (point[0] - params[i + 1]) / params[i + 2];
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

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            1.0);
            fitter.setCovarianceAdjusted(false);

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
            boolean valid = true;
            for (int i = 0; i < numParams; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }

            if (!valid) {
                continue;
            }

            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - numParams;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            numValid++;
            break;
        }
        assertTrue(numValid > 0);
    }

    @Test
    public void testFitUnidimensionalSineWithHoldAndFree()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
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
        
        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            y[i] = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                double[] initParams =  new double[SINE_UNI_PARAMS];
                double error;
                for(int i = 0; i < SINE_UNI_PARAMS; i++){
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) {
                double amplitude = params[0];
                double freq = params[1];
                double phase = params[2];
                double y = amplitude * Math.sin(freq * point[0] + phase);
                
                //derivative respect amplitude
                derivatives[0] = Math.sin(freq * point[0] + phase);
                //derivative respect frequency
                derivatives[1] = amplitude * Math.cos(freq * point[0] + phase) * point[0];
                //derivative respect phase
                derivatives[2] = amplitude * Math.cos(freq * point[0] + phase);
                
                return y;
            }
        };
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 
                1.0);
        fitter.setCovarianceAdjusted(false);
        
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

        double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalSineWithGradientEstimator()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
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

        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            double xi = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            y[i] = amplitude * Math.sin(freq * xi + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                    private double[] point;

                    private GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

                    @Override
                    public int getNumberOfDimensions() {
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
                    public double evaluate(int i, double[] point, double[] params,
                                           double[] derivatives) throws EvaluationException {
                        this.point = point;
                        double y = evaluateParams(point, params);
                        gradientEstimator.gradient(params, derivatives);

                        return y;
                    }

                    double evaluateParams(double[] point, double[] params) {
                        double amplitude = params[0];
                        double freq = params[1];
                        double phase = params[2];
                        return amplitude * Math.sin(freq * point[0] + phase);
                    }
                };

        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
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
        assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalGaussianWithGradientEstimator()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
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
        
        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            y[i] = 0.0;
            for (int k = 0; k < numgaussians; k++) {
                double b = params[k * GAUSS_UNI_PARAMS];
                double e = params[k * GAUSS_UNI_PARAMS + 1];
                double g = params[k * GAUSS_UNI_PARAMS + 2];
                y[i] += b * Math.exp(-Math.pow((x.getElementAt(i, 0) - e) / g, 2.0));
            }      
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                    
            private double[] point;
            
            private GradientEstimator gradientEstimator =
                    new GradientEstimator(
                            new MultiDimensionFunctionEvaluatorListener() {

                @Override
                public double evaluate(double[] params) {
                    return evaluateParams(point, params);
                }
            });

            @Override
            public int getNumberOfDimensions() {
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
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) throws EvaluationException {
                this.point = point;
                double y = evaluateParams(point, params);
                gradientEstimator.gradient(params, derivatives);
                
                return y;
            }

            double evaluateParams(double[] point, double[] params) {
                int i, na = params.length;
                double ex, arg;
                double y = 0.0;
                for (i = 0; i < na - 1; i += 3) {
                    arg = (point[0] - params[i + 1]) / params[i + 2];
                    ex = Math.exp(-Math.pow(arg, 2.0));
                    y += params[i] * ex;
                }

                return y;
            }
        };
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 
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
        for (int i = 0; i < numParams; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        double chiSqrDegreesOfFreedom = npoints - numParams;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitMultidimensionalSine() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
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

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 2);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));
                x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));

                y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                        Math.sin(freqy * x.getElementAt(i, 1) + phasey);
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        private double[] point;

                        private GradientEstimator gradientEstimator =
                                new GradientEstimator(
                                        new MultiDimensionFunctionEvaluatorListener() {

                                            @Override
                                            public double evaluate(double[] params) {
                                                return evaluateParams(point, params);
                                            }
                                        });

                        @Override
                        public int getNumberOfDimensions() {
                            return 2;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[SINE_MULTI_PARAMS];
                            double error;
                            for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) throws EvaluationException {
                            this.point = point;
                            double y = evaluateParams(point, params);
                            gradientEstimator.gradient(params, derivatives);

                            return y;
                        }

                        double evaluateParams(double[] point, double[] params) {
                            double amplitude = params[0];
                            double freqx = params[1];
                            double freqy = params[2];
                            double phasex = params[3];
                            double phasey = params[4];

                            return amplitude * Math.sin(freqx * point[0] + phasex) *
                                    Math.sin(freqy * point[1] + phasey);
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            1.0);
            fitter.setCovarianceAdjusted(false);

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
            boolean failed = false;
            for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    failed = true;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }

            if (failed) {
                continue;
            }

            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - SINE_MULTI_PARAMS;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitMultidimensionalGaussian() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
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

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 2);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                y[i] = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    double b = params[k * GAUSS_MULTI_PARAMS];
                    double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    y[i] += b * Math.exp(-(Math.pow(x.getElementAt(i, 0) -
                            ex, 2.0) + Math.pow(x.getElementAt(i, 1) - ey, 2.0)) /
                            Math.pow(gx * gy, 2.0));
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        private double[] point;

                        private GradientEstimator gradientEstimator =
                                new GradientEstimator(
                                        new MultiDimensionFunctionEvaluatorListener() {

                                            @Override
                                            public double evaluate(double[] params) {
                                                return evaluateParams(point, params);
                                            }
                                        });

                        @Override
                        public int getNumberOfDimensions() {
                            return NUM_DIMENSIONS;
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
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) throws EvaluationException {
                            this.point = point;
                            double y = evaluateParams(point, params);
                            gradientEstimator.gradient(params, derivatives);

                            return y;
                        }

                        double evaluateParams(double[] point, double[] params) {
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

                            return y;
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            1.0);
            fitter.setCovarianceAdjusted(false);

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

            double chiSqrDegreesOfFreedom = npoints - numParams;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitMultidimensionalSineRepeatInOneDimension()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
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
        
        double[] y = new double[npoints];
        Matrix x = new Matrix(npoints, 2);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, 
                    MAX_RANDOM_VALUE));
            
            y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }
        
        LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                    
            private double[] point;
            
            private GradientEstimator gradientEstimator = 
                    new GradientEstimator(
                            new MultiDimensionFunctionEvaluatorListener() {

                @Override
                public double evaluate(double[] params) {
                    return evaluateParams(point, params);
                }
            });

            @Override
            public int getNumberOfDimensions() {
                return 2;
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
            public double evaluate(int i, double[] point, double[] params, 
                    double[] derivatives) throws EvaluationException {
                this.point = point;
                double y = evaluateParams(point, params);
                gradientEstimator.gradient(params, derivatives);
                
                return y;
            }

            double evaluateParams(double[] point, double[] params) {
                double amplitude = params[0];
                double freqx = params[1];
                double phasex = params[2];

                return amplitude * Math.sin(freqx * point[0] + phasex);
            }
        };
        
        LevenbergMarquardtMultiDimensionFitter fitter =
                new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 
                1.0);
        fitter.setCovarianceAdjusted(false);
        
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

        double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        double chiSqr = fitter.getChisq();

        //probability that chi square can be smaller
        //(the smaller is p the better)
        double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        //measure of quality (1.0 indicates maximum quality)
        double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalConstantCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            double constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

            double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[CONSTANT_PARAMS];
            params[0] = constant;

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 1);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            NormalDist dist = new NormalDist();
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0,
                        randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

                //propagate standard deviation of a (sigma) for x[i] value:
                NormalDist.propagate(new NormalDist.DerivativeEvaluator() {
                    @Override
                    public double evaluate(double param) {
                        //y[i] = constant
                        return param;
                    }

                    @Override
                    public double evaluateDerivative(double param) {
                        //derivative respect param
                        return 1.0;
                    }
                }, params[0], sigma, dist);

                //expression below is equal to y[i] = constant;
                y[i] = dist.getMean();
                assertEquals(y[i], constant, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(sigmas[i]);

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                        @Override
                        public int getNumberOfDimensions() {
                            return 1;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[CONSTANT_PARAMS];
                            double error;
                            for (int i = 0; i < CONSTANT_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) {
                            double constant = params[0];

                            //derivative of evaluated function respect constant parameter
                            derivatives[0] = 1.0;

                            //evaluated function
                            return constant;
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }

            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, CONSTANT_PARAMS);
            for (int i = 0; i < CONSTANT_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], SMALL_ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), CONSTANT_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), CONSTANT_PARAMS);
            assertTrue(fitter.getChisq() > 0.0);

            double chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] point = new double[evaluator.getNumberOfDimensions()];
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            double standardDeviation3 = Math.sqrt(
                    cov.getElementAt(0, 0));

            assertEquals(standardDeviation3, sigma, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "real parameter: " + params[0] +
                    ", estimated parameter: " + fitter.getA()[0]);
            LOGGER.log(Level.INFO, "real parameter sigma: " + sigma +
                    ", estimated parameter sigma: " + standardDeviation3);

            numValid++;

            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitUnidimensionalLine1Covariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            double a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

            double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[LINE1_PARAMS];
            params[0] = a;

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 1);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            NormalDist dist = new NormalDist();
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                //propagate standard deviation of a (sigma) for x[i] value:
                NormalDist.propagate(new NormalDist.DerivativeEvaluator() {
                    @Override
                    public double evaluate(double param) {
                        //y[i] = a * x[i]
                        return param * xi;
                    }

                    @Override
                    public double evaluateDerivative(double param) {
                        //derivative respect param
                        return xi;
                    }
                }, params[0], sigma, dist);

                //expression below is equal to y[i] = a * x[i];
                y[i] = dist.getMean();
                assertEquals(y[i], a * xi, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(sigmas[i]);

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                        @Override
                        public int getNumberOfDimensions() {
                            return 1;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[LINE1_PARAMS];
                            double error;
                            for (int i = 0; i < LINE1_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) {
                            double a = params[0];

                            //derivative of evaluated function f(x) respect parameter a
                            derivatives[0] = point[0];

                            //evaluated function
                            return a * point[0];
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }

            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, LINE1_PARAMS);
            for (int i = 0; i < LINE1_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), LINE1_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), LINE1_PARAMS);
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - LINE1_PARAMS;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] point = new double[evaluator.getNumberOfDimensions()];
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            double standardDeviation3 = Math.sqrt(
                    cov.getElementAt(0, 0));

            assertEquals(standardDeviation3, sigma, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "real parameter: " + params[0] +
                    ", estimated parameter: " + fitter.getA()[0]);
            LOGGER.log(Level.INFO, "real parameter sigma: " + sigma +
                    ", estimated parameter sigma: " + standardDeviation3);

            numValid++;

            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitUnidimensionalLine2Covariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            double a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
            double b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

            double sigmaA = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaB = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double varianceA = sigmaA * sigmaA;
            double varianceB = sigmaB * sigmaB;

            final double[] params = new double[LINE2_PARAMS];
            params[0] = a;
            params[1] = b;

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 1);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            MultivariateNormalDist dist = new MultivariateNormalDist();
            Matrix covariance = Matrix.diagonal(new double[]{varianceA, varianceB});
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                //propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(double[] params, double[] y, Matrix jacobian) {
                        //y[i] = a * x[i] + b
                        y[0] = params[0] * xi + params[1];

                        //derivatives respect parameters
                        jacobian.setElementAt(0, 0, xi);
                        jacobian.setElementAt(0, 1, 1.0);
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }
                }, params, covariance, dist);

                //expression below is equal to y[i] = a * x[i] + b
                y[i] = dist.getMean()[0];
                assertEquals(y[i], a * xi + b, SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(sigmas[i]);

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                        @Override
                        public int getNumberOfDimensions() {
                            return 1;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[LINE2_PARAMS];
                            double error;
                            for (int i = 0; i < LINE2_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) {
                            double a = params[0];
                            double b = params[1];

                            //derivatives of function f(x) respect parameters a and b
                            derivatives[0] = point[0];
                            derivatives[1] = 1.0;

                            //evaluated function f(x) = a * x + b
                            return a * point[0] + b;
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, LINE2_PARAMS);
            for (int i = 0; i < LINE2_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), LINE2_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), LINE2_PARAMS);
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - LINE2_PARAMS;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] point = new double[evaluator.getNumberOfDimensions()];
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            double standardDeviationA = Math.sqrt(
                    cov.getElementAt(0, 0));
            double standardDeviationB = Math.sqrt(
                    cov.getElementAt(1, 1));

            LOGGER.log(Level.INFO, "real parameter A: " + params[0] +
                    ", estimated parameter A: " + fitter.getA()[0]);
            LOGGER.log(Level.INFO, "real parameter sigma A: " + sigmaA +
                    ", estimated parameter sigma A: " + standardDeviationA);


            LOGGER.log(Level.INFO, "real parameter B: " + params[1] +
                    ", estimated parameter B: " + fitter.getA()[1]);
            LOGGER.log(Level.INFO, "real parameter sigma B: " + sigmaB +
                    ", estimated parameter sigma B: " + standardDeviationB);

            assertEquals(standardDeviationA, sigmaA, SMALL_ABSOLUTE_ERROR);

            numValid++;

            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitUnidimensionalSineCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                    MAX_SINE_AMPLITUDE);
            double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            double sigmaAmplitude = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaFreq = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaPhase = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            double varianceAmplitude = sigmaAmplitude * sigmaAmplitude;
            double varianceFreq = sigmaFreq * sigmaFreq;
            double variancePhase = sigmaPhase * sigmaPhase;

            final double[] params = new double[SINE_UNI_PARAMS];
            params[0] = amplitude;
            params[1] = freq;
            params[2] = phase;

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 1);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            MultivariateNormalDist dist = new MultivariateNormalDist();
            Matrix covariance = Matrix.diagonal(new double[]{
                    varianceAmplitude, varianceFreq, variancePhase});
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                //propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(double[] params, double[] y, Matrix jacobian) {
                        //y[i] = amplitude * Math.sin(freq * xi + phase)
                        y[0] = params[0] * Math.sin(params[1] * xi + params[2]);

                        //derivatives respect parameters
                        jacobian.setElementAt(0, 0,
                                Math.sin(params[1] * xi + params[2]));
                        jacobian.setElementAt(0, 1,
                                params[0] * Math.cos(params[1] * xi + params[2]) * xi);
                        jacobian.setElementAt(0, 2,
                                params[0] * Math.cos(params[1] * xi + params[2]));
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }
                }, params, covariance, dist);

                //expression below is equal to y[i] = amplitude * Math.sin(freq * x[i] + phase)
                y[i] = dist.getMean()[0];
                assertEquals(y[i], amplitude * Math.sin(freq * xi + phase),
                        SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(sigmas[i]);

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        @Override
                        public int getNumberOfDimensions() {
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
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) {
                            double amplitude = params[0];
                            double freq = params[1];
                            double phase = params[2];
                            double y = amplitude * Math.sin(freq * point[0] + phase);

                            //derivative respect amplitude
                            derivatives[0] = Math.sin(freq * point[0] + phase);
                            //derivative respect frequency
                            derivatives[1] = amplitude * Math.cos(freq * point[0] + phase) * point[0];
                            //derivative respect phase
                            derivatives[2] = amplitude * Math.cos(freq * point[0] + phase);

                            return y;
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, SINE_UNI_PARAMS);
            boolean valid = true;
            for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), SINE_UNI_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), SINE_UNI_PARAMS);
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] point = new double[evaluator.getNumberOfDimensions()];
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            double standardDeviationAmplitude = Math.sqrt(
                    cov.getElementAt(0, 0));
            double standardDeviationFreq = Math.sqrt(
                    cov.getElementAt(1, 1));
            double standardDeviationPhase = Math.sqrt(
                    cov.getElementAt(2, 2));

            LOGGER.log(Level.INFO, "real parameter Amplitude: " + params[0] +
                    ", estimated parameter Amplitude: " + fitter.getA()[0]);
            LOGGER.log(Level.INFO, "real parameter sigma Amplitude: " + sigmaAmplitude +
                    ", estimated parameter sigma Amplitude: " + standardDeviationAmplitude);

            LOGGER.log(Level.INFO, "real parameter Freq: " + params[1] +
                    ", estimated parameter Freq: " + fitter.getA()[1]);
            LOGGER.log(Level.INFO, "real parameter sigma Freq: " + sigmaFreq +
                    ", estimated parameter sigma Freq: " + standardDeviationFreq);

            LOGGER.log(Level.INFO, "real parameter Phase: " + params[2] +
                    ", estimated parameter Phase: " + fitter.getA()[2]);
            LOGGER.log(Level.INFO, "real parameter sigma Phase: " + sigmaPhase +
                    ", estimated parameter sigma Phase: " + standardDeviationPhase);

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitUnidimensionalGaussianCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_UNI_PARAMS;

            double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            double[] varianceParams = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 1);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            MultivariateNormalDist dist = new MultivariateNormalDist();
            Matrix covariance = Matrix.diagonal(varianceParams);
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                //propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(double[] params, double[] y, Matrix jacobian) {
                        y[0] = 0.0;
                        for (int k = 0; k < numgaussians; k++) {
                            double b = params[k * GAUSS_UNI_PARAMS];
                            double e = params[k * GAUSS_UNI_PARAMS + 1];
                            double g = params[k * GAUSS_UNI_PARAMS + 2];
                            y[0] += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                        }

                        int i, na = params.length;
                        double fac, ex, arg;
                        for (i = 0; i < na - 1; i += 3) {
                            arg = (xi - params[i + 1]) / params[i + 2];
                            ex = Math.exp(-Math.pow(arg, 2.0));
                            fac = params[i] * ex * 2. * arg;

                            jacobian.setElementAt(0, i, ex);
                            jacobian.setElementAt(0, i + 1,
                                    fac / params[i + 2]);
                            jacobian.setElementAt(0, i + 2,
                                    fac * arg / params[i + 2]);
                        }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }
                }, params, covariance, dist);

                double yi = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    double b = params[k * GAUSS_UNI_PARAMS];
                    double e = params[k * GAUSS_UNI_PARAMS + 1];
                    double g = params[k * GAUSS_UNI_PARAMS + 2];
                    yi += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                }

                y[i] = dist.getMean()[0];
                assertEquals(y[i], yi, SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(sigmas[i]);

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        @Override
                        public int getNumberOfDimensions() {
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
                        public double evaluate(int pos, double[] point, double[] params,
                                               double[] derivatives) {
                            int i, na = params.length;
                            double fac, ex, arg;
                            double y = 0.0;
                            for (i = 0; i < na - 1; i += 3) {
                                arg = (point[0] - params[i + 1]) / params[i + 2];
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

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, numParams);
            boolean valid = true;
            for (int i = 0; i < numParams; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), numParams);
            assertEquals(fitter.getCovar().getColumns(), numParams);
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - numParams;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] point = new double[evaluator.getNumberOfDimensions()];
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            for (int i = 0; i < numParams; i++) {
                double standardDeviation = Math.sqrt(cov.getElementAt(i, i));

                LOGGER.log(Level.INFO, "real parameter: " + params[i] +
                        ", estimated parameter: " + fitter.getA()[i]);
                LOGGER.log(Level.INFO, "real parameter sigma: " + sigma +
                        ", estimated parameter sigma: " + standardDeviation);
            }

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitMultidimensionalSineCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                    MAX_SINE_AMPLITUDE);
            double freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            double freqy = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            double phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
            double phasey = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            double sigmaAmplitude = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaFreqx = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaFreqy = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaPhasex = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            double sigmaPhasey = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            double varianceAmplitude = sigmaAmplitude * sigmaAmplitude;
            double varianceFreqx = sigmaFreqx * sigmaFreqx;
            double varianceFreqy = sigmaFreqy * sigmaFreqy;
            double variancePhasex = sigmaPhasex * sigmaPhasex;
            double variancePhasey = sigmaPhasey * sigmaPhasey;

            final double[] params = new double[SINE_MULTI_PARAMS];
            params[0] = amplitude;
            params[1] = freqx;
            params[2] = freqy;
            params[3] = phasex;
            params[4] = phasey;

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 2);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final double[] point = new double[2];
            MultivariateNormalDist dist = new MultivariateNormalDist();
            Matrix covariance = Matrix.diagonal(new double[]{
                    varianceAmplitude, varianceFreqx, varianceFreqy, variancePhasex, variancePhasey});
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi0 = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double xi1 = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi0);
                x.setElementAt(i, 1, xi1);
                point[0] = xi0;
                point[1] = xi1;

                //propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {

                    double[] derivatives = new double[params.length];

                    private GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

                    @Override
                    public void evaluate(double[] params, double[] y, Matrix jacobian)  {
                        try {
                            gradientEstimator.gradient(params, derivatives);
                            y[0] = evaluateParams(point, params);
                            jacobian.fromArray(derivatives);
                        } catch (Exception ignore) { }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }

                    double evaluateParams(double[] point, double[] params) {
                        double amplitude = params[0];
                        double freqx = params[1];
                        double freqy = params[2];
                        double phasex = params[3];
                        double phasey = params[4];

                        return amplitude * Math.sin(freqx * point[0] + phasex) *
                                Math.sin(freqy * point[1] + phasey);
                    }

                }, params, covariance, dist);

                //expression below is equal to:
                //y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                //      Math.sin(freqy * x.getElementAt(i, 1) + phasey);
                y[i] = dist.getMean()[0];
                assertEquals(y[i], amplitude * Math.sin(freqx * xi0 + phasex) *
                        Math.sin(freqy * xi1 + phasey), SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(sigmas[i]);

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        private double[] point;

                        private GradientEstimator gradientEstimator =
                                new GradientEstimator(
                                        new MultiDimensionFunctionEvaluatorListener() {

                                            @Override
                                            public double evaluate(double[] params) {
                                                return evaluateParams(point, params);
                                            }
                                        });

                        @Override
                        public int getNumberOfDimensions() {
                            return 2;
                        }

                        @Override
                        public double[] createInitialParametersArray() {
                            double[] initParams = new double[SINE_MULTI_PARAMS];
                            double error;
                            for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) throws EvaluationException {
                            this.point = point;
                            double y = evaluateParams(point, params);
                            gradientEstimator.gradient(params, derivatives);

                            return y;
                        }

                        double evaluateParams(double[] point, double[] params) {
                            double amplitude = params[0];
                            double freqx = params[1];
                            double freqy = params[2];
                            double phasex = params[3];
                            double phasey = params[4];

                            return amplitude * Math.sin(freqx * point[0] + phasex) *
                                    Math.sin(freqy * point[1] + phasey);
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, SINE_MULTI_PARAMS);
            boolean valid = true;
            for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }

            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), SINE_MULTI_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), SINE_MULTI_PARAMS);
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - SINE_MULTI_PARAMS;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            double standardDeviationAmplitude = Math.sqrt(
                    cov.getElementAt(0, 0));
            double standardDeviationFreqx = Math.sqrt(
                    cov.getElementAt(1, 1));
            double standardDeviationFreqy = Math.sqrt(
                    cov.getElementAt(2, 2));
            double standardDeviationPhasex = Math.sqrt(
                    cov.getElementAt(3, 3));
            double standardDeviationPhasey = Math.sqrt(
                    cov.getElementAt(4,4));

            LOGGER.log(Level.INFO, "real parameter Amplitude: " + params[0] +
                    ", estimated parameter Amplitude: " + fitter.getA()[0]);
            LOGGER.log(Level.INFO, "real parameter sigma Amplitude: " + sigmaAmplitude +
                    ", estimated parameter sigma Amplitude: " + standardDeviationAmplitude);

            LOGGER.log(Level.INFO, "real parameter Freqx: " + params[1] +
                    ", estimated parameter Freqx: " + fitter.getA()[1]);
            LOGGER.log(Level.INFO, "real parameter sigma Freqx: " + sigmaFreqx +
                    ", estimated parameter sigma Freqx: " + standardDeviationFreqx);

            LOGGER.log(Level.INFO, "real parameter Freqy: " + params[2] +
                    ", estimated parameter Freqy: " + fitter.getA()[2]);
            LOGGER.log(Level.INFO, "real parameter sigma Freqy: " + sigmaFreqy +
                    ", estimated parameter sigma Freqy: " + standardDeviationFreqy);

            LOGGER.log(Level.INFO, "real parameter Phasex: " + params[3] +
                    ", estimated parameter Phasex: " + fitter.getA()[3]);
            LOGGER.log(Level.INFO, "real parameter sigma Phasex: " + sigmaPhasex +
                    ", estimated parameter sigma Phasex: " + standardDeviationPhasex);

            LOGGER.log(Level.INFO, "real parameter Phasey: " + params[4] +
                    ", estimated parameter Phasey: " + fitter.getA()[4]);
            LOGGER.log(Level.INFO, "real parameter sigma Phasey: " + sigmaPhasey +
                    ", estimated parameter sigma Phasey: " + standardDeviationPhasey);

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitMultidimensionalGaussianCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            UniformRandomizer randomizer = new UniformRandomizer(new Random());

            int npoints = N_SAMPLES;
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS,
                    MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_MULTI_PARAMS;

            double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            double[] varianceParams = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            double[] y = new double[npoints];
            Matrix x = new Matrix(npoints, 2);
            double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final double[] point = new double[2];
            MultivariateNormalDist dist = new MultivariateNormalDist();
            Matrix covariance = Matrix.diagonal(varianceParams);
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi0 = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double xi1 = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi0);
                x.setElementAt(i, 1, xi1);
                point[0] = xi0;
                point[1] = xi1;

                //propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {

                    double[] derivatives = new double[params.length];

                    private GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

                    @Override
                    public void evaluate(double[] params, double[] y, Matrix jacobian)  {
                        try {
                            gradientEstimator.gradient(params, derivatives);
                            y[0] = evaluateParams(point, params);
                            jacobian.fromArray(derivatives);
                        } catch (Exception ignore) { }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }

                    double evaluateParams(double[] point, double[] params) {
                        double y = 0.0;
                        for (int k = 0; k < numgaussians; k++) {
                            double b = params[k * GAUSS_MULTI_PARAMS];
                            double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                            double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                            double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                            double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                            y += b * Math.exp(-(Math.pow(point[0] -
                                    ex, 2.0) + Math.pow(point[1] - ey, 2.0)) /
                                    Math.pow(gx * gy, 2.0));
                        }

                        return y;
                    }

                }, params, covariance, dist);

                double yi = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    double b = params[k * GAUSS_MULTI_PARAMS];
                    double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    yi += b * Math.exp(-(Math.pow(x.getElementAt(i, 0) -
                            ex, 2.0) + Math.pow(x.getElementAt(i, 1) - ey, 2.0)) /
                            Math.pow(gx * gy, 2.0));
                }

                y[i] = dist.getMean()[0];
                assertEquals(y[i], yi, SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            LevenbergMarquardtMultiDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                        private double[] point;

                        private GradientEstimator gradientEstimator =
                                new GradientEstimator(
                                        new MultiDimensionFunctionEvaluatorListener() {

                                            @Override
                                            public double evaluate(double[] params) {
                                                return evaluateParams(point, params);
                                            }
                                        });

                        @Override
                        public int getNumberOfDimensions() {
                            return NUM_DIMENSIONS;
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
                        public double evaluate(int i, double[] point, double[] params,
                                               double[] derivatives) throws EvaluationException {
                            this.point = point;
                            double y = evaluateParams(point, params);
                            gradientEstimator.gradient(params, derivatives);

                            return y;
                        }

                        double evaluateParams(double[] point, double[] params) {
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

                            return y;
                        }
                    };

            LevenbergMarquardtMultiDimensionFitter fitter =
                    new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            //check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            //fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            //check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, numParams);
            boolean valid = true;
            for (int i = 0; i < numParams; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > LARGE_ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], LARGE_ABSOLUTE_ERROR);
            }

            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), numParams);
            assertEquals(fitter.getCovar().getColumns(), numParams);
            assertTrue(fitter.getChisq() > 0);

            double chiSqrDegreesOfFreedom = npoints - numParams;
            double chiSqr = fitter.getChisq();

            //probability that chi square can be smaller
            //(the smaller is p the better)
            double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            //measure of quality (1.0 indicates maximum quality)
            double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            //Covariance of parameters is
            //Vp = (J’ * W * J)^-1

            //where J is the jacobian of measures, each row contains derivatives
            //for one sample, each column is the derivative for one parameter
            //W is the inverse of the measurement error covariance. Assuming
            //independent samples, W is diagonal with diagonal terms 1/sigma^2
            //where sigma is the standard deviation of each sample
            //More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            Matrix invCov = new Matrix(params.length, params.length);
            Matrix tmp1 = new Matrix(params.length, 1);
            Matrix tmp2 = new Matrix(1, params.length);
            Matrix tmpInvCov = new Matrix(params.length, params.length);
            double[] estimatedParams = fitter.getA();
            double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                double yi = evaluator.evaluate(i, point, estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            for (int i = 0; i < numParams; i++) {
                double standardDeviation = Math.sqrt(cov.getElementAt(i, i));

                LOGGER.log(Level.INFO, "real parameter: " + params[i] +
                        ", estimated parameter: " + fitter.getA()[i]);
                LOGGER.log(Level.INFO, "real parameter sigma: " + sigma +
                        ", estimated parameter sigma: " + standardDeviation);
            }

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }
}
