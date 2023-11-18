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
import com.irurueta.numerical.JacobianEstimator;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
import com.irurueta.numerical.MultiVariateFunctionEvaluatorListener;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.ChiSqDist;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.MaxIterationsExceededException;
import com.irurueta.statistics.MultivariateNormalDist;
import com.irurueta.statistics.NormalDist;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.Assert.*;

public class LevenbergMarquardtMultiVariateFitterTest {

    private static final Logger LOGGER = Logger.getLogger(
            LevenbergMarquardtMultiVariateFitterTest.class.getName());

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
    // amplitude, freqx, freqy, phasex, phasey
    private static final int SINE_MULTI_PARAMS = 5;
    private static final double MIN_SINE_AMPLITUDE = 0.5;
    private static final double MAX_SINE_AMPLITUDE = 10.0;
    private static final double MIN_SINE_FREQ = 0.5;
    private static final double MAX_SINE_FREQ = 100.0;
    private static final double MIN_SINE_PHASE = -Math.PI;
    private static final double MAX_SINE_PHASE = Math.PI;

    private static final int GAUSS_UNI_PARAMS = 3;
    // B, Ex, Ey, Gx, Gy
    private static final int GAUSS_MULTI_PARAMS = 5;
    private static final int MIN_GAUSSIANS = 1;
    private static final int MAX_GAUSSIANS = 3;

    private static final int NUM_DIMENSIONS = 2;
    private static final int NUM_VARIABLES = 2;

    private static final int TIMES = 10;
    private static final int N_SAMPLES = 1000000;

    @Test
    public void testConstructor() throws FittingException, WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // test constructor with input data
        final Matrix x = new Matrix(nPoints, 2);
        final Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        final double[] sig = new double[nPoints];

        fitter = new LevenbergMarquardtMultiVariateFitter(x, y, sig);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        final Matrix shortX = new Matrix(nPoints - 1, 2);
        final Matrix shortY = new Matrix(nPoints - 1, NUM_VARIABLES);
        final double[] shortSig = new double[nPoints - 1];

        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtMultiVariateFitter(x, y, 1.0);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator
        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                    }
                };

        fitter = new LevenbergMarquardtMultiVariateFitter(evaluator);

        // check default values
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
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        // test constructor with evaluator and input data
        fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, sig);

        // check default values
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
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, shortX,
                    y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x,
                    shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator and input data (constant sigma)
        fitter = new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, nPoints);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), nPoints);
        assertEquals(fitter.getCovar().getColumns(), nPoints);
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator,
                    shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtMultiVariateFitter(evaluator,
                    x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);
    }

    @Test
    public void testGetSetNdone() {
        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default value
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_NDONE,
                fitter.getNdone());

        // new value
        fitter.setNdone(5);

        // check correctness
        assertEquals(5, fitter.getNdone());

        // force IllegalArgumentException
        try {
            fitter.setNdone(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetItmax() {
        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default value
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_ITMAX,
                fitter.getItmax());

        // new value
        fitter.setItmax(10);

        // check correctness
        assertEquals(10, fitter.getItmax());

        // force IllegalArgumentException
        try {
            fitter.setItmax(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetTol() {
        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default value
        assertEquals(LevenbergMarquardtMultiVariateFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);

        // new value
        fitter.setTol(1e-1);

        // check correctness
        assertEquals(1e-1, fitter.getTol(), 0.0);

        // force IllegalArgumentException
        try {
            fitter.setTol(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetFunctionEvaluator() throws FittingException {
        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default value
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertNull(fitter.getAlpha());

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                    }
                };

        // new value
        fitter.setFunctionEvaluator(evaluator);

        // check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getA());
        assertEquals(GAUSS_UNI_PARAMS, fitter.getA().length);
        assertNotNull(fitter.getCovar());
        assertEquals(GAUSS_UNI_PARAMS, fitter.getCovar().getRows());
        assertEquals(GAUSS_UNI_PARAMS, fitter.getCovar().getColumns());
        assertNotNull(fitter.getAlpha());
        assertEquals(GAUSS_UNI_PARAMS, fitter.getAlpha().getRows());
        assertEquals(GAUSS_UNI_PARAMS, fitter.getAlpha().getColumns());
    }

    @Test
    public void testGetSetInputData() throws WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final Matrix x = new Matrix(nPoints, 2);
        final Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        final double[] sig = new double[nPoints];

        // set input data
        fitter.setInputData(x, y, sig);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);

        // Force IllegalArgumentException
        final Matrix shortX = new Matrix(nPoints - 1, 2);
        final Matrix shortY = new Matrix(nPoints - 1, NUM_VARIABLES);
        final double[] shortSig = new double[nPoints - 1];

        try {
            fitter.setInputData(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetInputDataWithConstantSigma()
            throws WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final Matrix x = new Matrix(nPoints, 2);
        final Matrix y = new Matrix(nPoints, NUM_VARIABLES);

        // set input data
        fitter.setInputData(x, y, 1.0);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }


        // Force IllegalArgumentException
        final Matrix shortX = new Matrix(nPoints - 1, 2);
        final Matrix shortY = new Matrix(nPoints - 1, NUM_VARIABLES);

        try {
            fitter.setInputData(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsReady() throws WrongSizeException, FittingException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final Matrix x = new Matrix(nPoints, 2);
        final Matrix y = new Matrix(nPoints, NUM_VARIABLES);
        final double[] sig = new double[nPoints];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                    }
                };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());
    }

    @Test
    public void testIsSetCovarianceAdjusted() {
        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter();

        // check default value
        assertTrue(fitter.isCovarianceAdjusted());

        // set new value
        fitter.setCovarianceAdjusted(false);

        // check
        assertFalse(fitter.isCovarianceAdjusted());
    }

    @Test
    public void testFitUnidimensionalConstant() throws WrongSizeException,
            FittingException, NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[CONSTANT_PARAMS];
        params[0] = constant;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0,
                    randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            double value = constant;
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                        final double[] initParams = new double[CONSTANT_PARAMS];
                        double error;
                        for (int i = 0; i < CONSTANT_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                        final double constant = params[0];

                        // derivative of evaluated function respect constant parameter
                        jacobian.setElementAt(0, 0, 1.0);

                        // evaluated function
                        result[0] = constant;
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(CONSTANT_PARAMS, fitter.getA().length);
        for (int i = 0; i < CONSTANT_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalLine1() throws WrongSizeException,
            FittingException, NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[LINE1_PARAMS];
        params[0] = a;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            double value = a * xi;
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                        final double[] initParams = new double[LINE1_PARAMS];
                        double error;
                        for (int i = 0; i < LINE1_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                        final double a = params[0];

                        // derivative of evaluated function respect constant parameter
                        jacobian.setElementAt(0, 0, a);

                        // evaluated function
                        result[0] = a * point[0];
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(LINE1_PARAMS, fitter.getA().length);
        for (int i = 0; i < LINE1_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - LINE1_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalLine2() throws WrongSizeException,
            FittingException, NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
        final double b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[LINE2_PARAMS];
        params[0] = a;
        params[1] = b;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 1);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            double value = a * xi + b;
            error = errorRandomizer.nextDouble();
            value += error;
            y.setElementAt(i, 0, value);
        }

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                        final double[] initParams = new double[LINE2_PARAMS];
                        double error;
                        for (int i = 0; i < LINE2_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                        final double a = params[0];
                        final double b = params[1];

                        // derivatives of function f(x) respect parameters a and b
                        jacobian.setElementAt(0, 0, point[0]);
                        jacobian.setElementAt(0, 1, 1.0);

                        // evaluated function f(x) = a * x + b
                        result[0] = a * point[0] + b;
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(LINE2_PARAMS, fitter.getA().length);
        for (int i = 0; i < LINE2_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - LINE2_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalSine() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 1);
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

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                        final double[] initParams = new double[SINE_UNI_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                        final double amplitude = params[0];
                        final double freq = params[1];
                        final double phase = params[2];
                        final double y = amplitude * Math.sin(freq * point[0] + phase);

                        // derivative respect amplitude
                        jacobian.setElementAt(0, 0, Math.sin(freq * point[0] + phase));
                        jacobian.setElementAt(0, 1, amplitude * Math.cos(
                                freq * point[0] + phase) * point[0]);
                        jacobian.setElementAt(0, 2, amplitude * Math.cos(
                                freq * point[0] + phase));

                        result[0] = y;
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_UNI_PARAMS, fitter.getA().length);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_UNI_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));
                double value = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    final double b = params[k * GAUSS_UNI_PARAMS];
                    final double e = params[k * GAUSS_UNI_PARAMS + 1];
                    final double g = params[k * GAUSS_UNI_PARAMS + 2];
                    value += b * Math.exp(-Math.pow((x.getElementAt(i, 0) - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                            final double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int pos, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) {
                            int i;
                            final int na = params.length;
                            double fac;
                            double ex;
                            double arg;
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

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            1.0);
            fitter.setCovarianceAdjusted(false);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            fitter.fit();

            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(fitter.getA().length, numParams);
            boolean failed = false;
            for (int i = 0; i < numParams; i++) {
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

            final double chiSqrDegreesOfFreedom = npoints - numParams;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
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
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 1);
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

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                        final double[] initParams = new double[SINE_UNI_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(final int i, final double[] point, final double[] result,
                                         final double[] params, final Matrix jacobian) {
                        final double amplitude = params[0];
                        final double freq = params[1];
                        final double phase = params[2];
                        final double y = amplitude * Math.sin(freq * point[0] + phase);

                        // derivative respect amplitude
                        jacobian.setElementAt(0, 0, Math.sin(freq * point[0] + phase));
                        jacobian.setElementAt(0, 1, amplitude * Math.cos(
                                freq * point[0] + phase) * point[0]);
                        jacobian.setElementAt(0, 2, amplitude * Math.cos(
                                freq * point[0] + phase));

                        result[0] = y;
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // hold first parameter
        fitter.hold(0, params[0]);

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_UNI_PARAMS, fitter.getA().length);
        // first parameter is hold and matches exactly
        assertEquals(fitter.getA()[0], params[0], 0.0);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        // release first parameter
        fitter.free(0);

        // fit and check correctness
        fitter.fit();

        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_UNI_PARAMS, fitter.getA().length);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalGaussianWithJacobianEstimator()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {

        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_UNI_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE));
                double value = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    final double b = params[k * GAUSS_UNI_PARAMS];
                    final double e = params[k * GAUSS_UNI_PARAMS + 1];
                    final double g = params[k * GAUSS_UNI_PARAMS + 2];
                    value += b * Math.exp(-Math.pow((x.getElementAt(i, 0) - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private final JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(
                                                    final double[] params, final double[] result) {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 1;
                                            }
                                        });

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
                            final double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian)
                                throws EvaluationException {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }

                        void evaluateParams(final double[] point, final double[] params,
                                            final double[] result) {
                            int i;
                            final int na = params.length;
                            double ex;
                            double arg;
                            double y = 0.0;
                            for (i = 0; i < na - 1; i += 3) {
                                arg = (point[0] - params[i + 1]) / params[i + 2];
                                ex = Math.exp(-Math.pow(arg, 2.0));
                                y += params[i] * ex;
                            }

                            result[0] = y;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            fitter.fit();

            // check correctness
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

            final double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            if (allValid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitUnidimensionalSineWithJacobianEstimator()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 1);
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

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                    private double[] point;

                    private final JacobianEstimator jacobianEstimator =
                            new JacobianEstimator(
                                    new MultiVariateFunctionEvaluatorListener() {

                                        @Override
                                        public void evaluate(
                                                final double[] params, final double[] result) {
                                            evaluateParams(point, params, result);
                                        }

                                        @Override
                                        public int getNumberOfVariables() {
                                            return 1;
                                        }
                                    });

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
                        final double[] initParams = new double[SINE_UNI_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(
                            final int i, final double[] point, final double[] result,
                            final double[] params, final Matrix jacobian) throws EvaluationException {
                        this.point = point;
                        evaluateParams(point, params, result);
                        jacobianEstimator.jacobian(params, jacobian);
                    }

                    void evaluateParams(final double[] point, final double[] params,
                                        final double[] result) {
                        final double amplitude = params[0];
                        final double freq = params[1];
                        final double phase = params[2];
                        result[0] = amplitude * Math.sin(freq * point[0] + phase);
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        //check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_UNI_PARAMS, fitter.getA().length);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitMultidimensionalSine() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double freqy = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
        final double phasey = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_MULTI_PARAMS];
        params[0] = amplitude;
        params[1] = freqx;
        params[2] = freqy;
        params[3] = phasex;
        params[4] = phasey;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 2);
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

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                    private double[] point;

                    private final JacobianEstimator jacobianEstimator =
                            new JacobianEstimator(
                                    new MultiVariateFunctionEvaluatorListener() {

                                        @Override
                                        public void evaluate(
                                                final double[] params, final double[] result) {
                                            evaluateParams(point, params, result);
                                        }

                                        @Override
                                        public int getNumberOfVariables() {
                                            return 1;
                                        }
                                    });

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
                        final double[] initParams = new double[SINE_MULTI_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(
                            final int i, final double[] point, final double[] result,
                            final double[] params, final Matrix jacobian) throws EvaluationException {
                        this.point = point;
                        evaluateParams(point, params, result);
                        jacobianEstimator.jacobian(params, jacobian);
                    }

                    void evaluateParams(final double[] point, final double[] params,
                                        final double[] result) {
                        final double amplitude = params[0];
                        final double freqx = params[1];
                        final double freqy = params[2];
                        final double phasex = params[3];
                        final double phasey = params[4];

                        final double value = amplitude * Math.sin(freqx * point[0] + phasex) *
                                Math.sin(freqy * point[1] + phasey);

                        result[0] = value;
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_MULTI_PARAMS, fitter.getA().length);
        for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - SINE_MULTI_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitMultidimensionalGaussian() throws FittingException,
            NotReadyException, WrongSizeException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS,
                    MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_MULTI_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 2);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                double value = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    final double b = params[k * GAUSS_MULTI_PARAMS];
                    final double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    final double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    final double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    final double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    value += b * Math.exp(-(Math.pow(x.getElementAt(i, 0) - ex, 2.0) +
                            Math.pow(x.getElementAt(i, 1) - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                }
                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private final JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(
                                                    final double[] params, final double[] result) {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 1;
                                            }
                                        });

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
                            final double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) throws EvaluationException {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }

                        void evaluateParams(final double[] point, final double[] params,
                                            final double[] result) {
                            double y = 0.0;
                            for (int k = 0; k < numgaussians; k++) {
                                final double b = params[k * GAUSS_MULTI_PARAMS];
                                final double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                                final double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                                final double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                                final double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                                y += b * Math.exp(-(Math.pow(point[0] - ex, 2.0) +
                                        Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                            }

                            result[0] = y;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            fitter.fit();

            // check correctness
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
    }

    @Test
    public void testFitMultidimensionalSineRepeatInOneDimension()
            throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freqx;
        params[2] = phasex;

        final Matrix y = new Matrix(npoints, 1);
        final Matrix x = new Matrix(npoints, 2);
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

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                    private double[] point;

                    private final JacobianEstimator jacobianEstimator =
                            new JacobianEstimator(
                                    new MultiVariateFunctionEvaluatorListener() {

                                        @Override
                                        public void evaluate(
                                                final double[] params, final double[] result) {
                                            evaluateParams(point, params, result);
                                        }

                                        @Override
                                        public int getNumberOfVariables() {
                                            return 1;
                                        }
                                    });

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
                        final double[] initParams = new double[SINE_UNI_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(
                            final int i, final double[] point, final double[] result,
                            final double[] params, final Matrix jacobian) throws EvaluationException {
                        this.point = point;
                        evaluateParams(point, params, result);
                        jacobianEstimator.jacobian(params, jacobian);
                    }

                    void evaluateParams(final double[] point, final double[] params,
                                        final double[] result) {
                        final double amplitude = params[0];
                        final double freqx = params[1];
                        final double phasex = params[2];

                        result[0] = amplitude * Math.sin(freqx * point[0] + phasex);
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_UNI_PARAMS, fitter.getA().length);
        for (int i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitMultiVariateGaussianAndSine() throws FittingException,
            NotReadyException, WrongSizeException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
        final int numParams = SINE_MULTI_PARAMS;

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
        }

        final Matrix y = new Matrix(npoints, 2);
        final Matrix x = new Matrix(npoints, 2);
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
            x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE));
            double valueGaussian;
            double valueSine;

            final double b = params[0];
            final double ex = params[1];
            final double ey = params[2];
            final double gx = params[3];
            final double gy = params[4];
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

        final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                    private double[] point;

                    private final JacobianEstimator jacobianEstimator =
                            new JacobianEstimator(
                                    new MultiVariateFunctionEvaluatorListener() {

                                        @Override
                                        public void evaluate(
                                                final double[] params, final double[] result) {
                                            evaluateParams(point, params, result);
                                        }

                                        @Override
                                        public int getNumberOfVariables() {
                                            return 2;
                                        }
                                    });

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
                        final double[] initParams = new double[numParams];
                        double error;
                        for (int i = 0; i < numParams; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public void evaluate(
                            final int i, final double[] point, final double[] result,
                            final double[] params, final Matrix jacobian) throws EvaluationException {
                        this.point = point;
                        evaluateParams(point, params, result);
                        jacobianEstimator.jacobian(params, jacobian);
                    }

                    void evaluateParams(final double[] point, final double[] params,
                                        final double[] result) {
                        final double b = params[0];
                        final double ex = params[1];
                        final double ey = params[2];
                        final double gx = params[3];
                        final double gy = params[4];

                        final double amplitude = params[0];
                        final double freqx = params[1];
                        final double freqy = params[2];
                        final double phasex = params[3];
                        final double phasey = params[4];

                        final double valueGaussian = b * Math.exp(-(Math.pow(point[0] - ex, 2.0) +
                                Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));

                        final double valueSine = amplitude * Math.sin(freqx * point[0] + phasex) *
                                Math.sin(freqy * point[1] + phasey);

                        result[0] = valueGaussian;
                        result[1] = valueSine;
                    }
                };

        final LevenbergMarquardtMultiVariateFitter fitter =
                new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, 1.0);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(SINE_MULTI_PARAMS, fitter.getA().length);
        for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final double chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final double chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final double q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                ", probability smaller chi sqr: " + p * 100.0 +
                "%, quality: " + q * 100.0 + "%");
    }

    @Test
    public void testFitUnidimensionalConstantCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[CONSTANT_PARAMS];
            params[0] = constant;

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            final NormalDist dist = new NormalDist();
            double error;
            for (int i = 0; i < npoints; i++) {
                x.setElementAt(i, 0,
                        randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

                // propagate standard deviation of a (sigma) for x[i] value:
                NormalDist.propagate(new NormalDist.DerivativeEvaluator() {
                    @Override
                    public double evaluate(final double param) {
                        // y[i] = constant
                        return param;
                    }

                    @Override
                    public double evaluateDerivative(final double param) {
                        // derivative respect param
                        return 1.0;
                    }
                }, params[0], sigma, dist);

                // expression below is equal to y[i] = constant;
                double value = dist.getMean();
                assertEquals(value, constant, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                            final double[] initParams = new double[CONSTANT_PARAMS];
                            double error;
                            for (int i = 0; i < CONSTANT_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) {
                            final double constant = params[0];

                            // derivative of evaluated function respect constant parameter
                            jacobian.setElementAt(0, 0, 1.0);

                            // evaluated function
                            result[0] = constant;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }

            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(CONSTANT_PARAMS, fitter.getA().length);
            for (int i = 0; i < CONSTANT_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(CONSTANT_PARAMS, fitter.getCovar().getRows());
            assertEquals(CONSTANT_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final double chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] point = new double[evaluator.getNumberOfDimensions()];
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            final double standardDeviation3 = Math.sqrt(
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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[LINE1_PARAMS];
            params[0] = a;

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            NormalDist dist = new NormalDist();
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                NormalDist.propagate(new NormalDist.DerivativeEvaluator() {
                    @Override
                    public double evaluate(final double param) {
                        // y[i] = a * x[i]
                        return param * xi;
                    }

                    @Override
                    public double evaluateDerivative(final double param) {
                        // derivative respect param
                        return xi;
                    }
                }, params[0], sigma, dist);

                // expression below is equal to y[i] = a * x[i];
                double value = dist.getMean();
                assertEquals(value, a * xi, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                            final double[] initParams = new double[LINE1_PARAMS];
                            double error;
                            for (int i = 0; i < LINE1_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) {
                            final double a = params[0];

                            // derivative of evaluated function respect constant parameter
                            jacobian.setElementAt(0, 0, point[0]);

                            // evaluated function
                            result[0] = a * point[0];
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }

            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(LINE1_PARAMS, fitter.getA().length);
            for (int i = 0; i < LINE1_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(CONSTANT_PARAMS, fitter.getCovar().getRows());
            assertEquals(CONSTANT_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final double chiSqrDegreesOfFreedom = npoints - LINE1_PARAMS;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] point = new double[evaluator.getNumberOfDimensions()];
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            final double standardDeviation3 = Math.sqrt(
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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
            final double b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

            final double sigmaA = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaB = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double varianceA = sigmaA * sigmaA;
            final double varianceB = sigmaB * sigmaB;

            final double[] params = new double[LINE2_PARAMS];
            params[0] = a;
            params[1] = b;

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(new double[]{varianceA, varianceB});
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(
                            final double[] params, final double[] y, final Matrix jacobian) {
                        // y[i] = a * x[i] + b
                        y[0] = params[0] * xi + params[1];

                        // derivatives respect parameters
                        jacobian.setElementAt(0, 0, xi);
                        jacobian.setElementAt(0, 1, 1.0);
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }
                }, params, covariance, dist);

                // expression below is equal to y[i] = a * x[i] + b
                double value = dist.getMean()[0];
                assertEquals(value, a * xi + b, SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                            final double[] initParams = new double[LINE2_PARAMS];
                            double error;
                            for (int i = 0; i < LINE2_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) {
                            final double a = params[0];
                            final double b = params[1];

                            // derivatives of function f(x) respect parameters a and b
                            jacobian.setElementAt(0, 0, point[0]);
                            jacobian.setElementAt(0, 1, 1.0);

                            // evaluated function f(x) = a * x + b
                            result[0] = a * point[0] + b;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (final FittingException e) {
                continue;
            }

            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(LINE2_PARAMS, fitter.getA().length);
            for (int i = 0; i < LINE2_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(LINE2_PARAMS, fitter.getCovar().getRows());
            assertEquals(LINE2_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final double chiSqrDegreesOfFreedom = npoints - LINE2_PARAMS;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] point = new double[evaluator.getNumberOfDimensions()];
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            final double standardDeviationA = Math.sqrt(
                    cov.getElementAt(0, 0));
            final double standardDeviationB = Math.sqrt(
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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                    MAX_SINE_AMPLITUDE);
            final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            final double sigmaAmplitude = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaFreq = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaPhase = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double varianceAmplitude = sigmaAmplitude * sigmaAmplitude;
            final double varianceFreq = sigmaFreq * sigmaFreq;
            final double variancePhase = sigmaPhase * sigmaPhase;

            final double[] params = new double[SINE_UNI_PARAMS];
            params[0] = amplitude;
            params[1] = freq;
            params[2] = phase;

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(new double[]{
                    varianceAmplitude, varianceFreq, variancePhase});
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(
                            final double[] params, final double[] y, final Matrix jacobian) {
                        // y[i] = amplitude * Math.sin(freq * xi + phase)
                        y[0] = params[0] * Math.sin(params[1] * xi + params[2]);

                        // derivatives respect parameters
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

                // expression below is equal to y[i] = amplitude * Math.sin(freq * x[i] + phase)
                double value = dist.getMean()[0];
                assertEquals(value, amplitude * Math.sin(freq * xi + phase),
                        SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                            final double[] initParams = new double[SINE_UNI_PARAMS];
                            double error;
                            for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) {
                            final double amplitude = params[0];
                            final double freq = params[1];
                            final double phase = params[2];
                            final double y = amplitude * Math.sin(freq * point[0] + phase);

                            // derivative respect amplitude
                            jacobian.setElementAt(0, 0, Math.sin(freq * point[0] + phase));
                            jacobian.setElementAt(0, 1, amplitude * Math.cos(
                                    freq * point[0] + phase) * point[0]);
                            jacobian.setElementAt(0, 2, amplitude * Math.cos(
                                    freq * point[0] + phase));

                            result[0] = y;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (final FittingException e) {
                continue;
            }


            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(SINE_UNI_PARAMS, fitter.getA().length);
            boolean valid = true;
            for (int i = 0; i < SINE_UNI_PARAMS; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(SINE_UNI_PARAMS, fitter.getCovar().getRows());
            assertEquals(SINE_UNI_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final double chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] point = new double[evaluator.getNumberOfDimensions()];
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            final double standardDeviationAmplitude = Math.sqrt(
                    cov.getElementAt(0, 0));
            final double standardDeviationFreq = Math.sqrt(
                    cov.getElementAt(1, 1));
            final double standardDeviationPhase = Math.sqrt(
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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_UNI_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            final double[] varianceParams = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 1);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(varianceParams);
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(
                            final double[] params, final double[] y, final Matrix jacobian) {
                        y[0] = 0.0;
                        for (int k = 0; k < numgaussians; k++) {
                            final double b = params[k * GAUSS_UNI_PARAMS];
                            final double e = params[k * GAUSS_UNI_PARAMS + 1];
                            final double g = params[k * GAUSS_UNI_PARAMS + 2];
                            y[0] += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                        }

                        int i;
                        final int na = params.length;
                        double fac;
                        double ex;
                        double arg;
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
                    final double b = params[k * GAUSS_UNI_PARAMS];
                    final double e = params[k * GAUSS_UNI_PARAMS + 1];
                    final double g = params[k * GAUSS_UNI_PARAMS + 2];
                    yi += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                }

                double value = dist.getMean()[0];
                assertEquals(value, yi, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));
                if (sigmas[i] == 0.0) {
                    sigmas[i] = Double.MIN_VALUE;
                }

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
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
                            final double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int pos, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian) {
                            int i;
                            final int na = params.length;
                            double fac;
                            double ex;
                            double arg;
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

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (final FittingException e) {
                continue;
            }

            // check correctness
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

            final double chiSqrDegreesOfFreedom = npoints - numParams;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] point = new double[evaluator.getNumberOfDimensions()];
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            for (int i = 0; i < numParams; i++) {
                final double standardDeviation = Math.sqrt(cov.getElementAt(i, i));

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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                    MAX_SINE_AMPLITUDE);
            final double freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final double freqy = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final double phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
            final double phasey = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            final double sigmaAmplitude = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaFreqx = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaFreqy = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaPhasex = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final double sigmaPhasey = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double varianceAmplitude = sigmaAmplitude * sigmaAmplitude;
            final double varianceFreqx = sigmaFreqx * sigmaFreqx;
            final double varianceFreqy = sigmaFreqy * sigmaFreqy;
            final double variancePhasex = sigmaPhasex * sigmaPhasex;
            final double variancePhasey = sigmaPhasey * sigmaPhasey;

            final double[] params = new double[SINE_MULTI_PARAMS];
            params[0] = amplitude;
            params[1] = freqx;
            params[2] = freqy;
            params[3] = phasex;
            params[4] = phasey;

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 2);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final double[] point = new double[2];
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(new double[]{
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

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {

                    final double[] derivatives = new double[params.length];

                    private final GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(final double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

                    @Override
                    public void evaluate(
                            final double[] params, final double[] y, final Matrix jacobian) {
                        try {
                            gradientEstimator.gradient(params, derivatives);
                            y[0] = evaluateParams(point, params);
                            jacobian.fromArray(derivatives);
                        } catch (final Exception ignore) {
                        }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }

                    double evaluateParams(final double[] point, final double[] params) {
                        final double amplitude = params[0];
                        final double freqx = params[1];
                        final double freqy = params[2];
                        final double phasex = params[3];
                        final double phasey = params[4];

                        return amplitude * Math.sin(freqx * point[0] + phasex) *
                                Math.sin(freqy * point[1] + phasey);
                    }

                }, params, covariance, dist);

                // expression below is equal to:
                // y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                //      Math.sin(freqy * x.getElementAt(i, 1) + phasey)
                double value = dist.getMean()[0];
                assertEquals(value, amplitude * Math.sin(freqx * xi0 + phasex) *
                        Math.sin(freqy * xi1 + phasey), SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private final JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(
                                                    final double[] params, final double[] result) {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 1;
                                            }
                                        });

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
                            final double[] initParams = new double[SINE_MULTI_PARAMS];
                            double error;
                            for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian)
                                throws EvaluationException {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }

                        void evaluateParams(final double[] point, final double[] params,
                                            final double[] result) {
                            final double amplitude = params[0];
                            final double freqx = params[1];
                            final double freqy = params[2];
                            final double phasex = params[3];
                            final double phasey = params[4];

                            final double value = amplitude * Math.sin(freqx * point[0] + phasex) *
                                    Math.sin(freqy * point[1] + phasey);

                            result[0] = value;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertNotNull(fitter.getA());
            assertEquals(SINE_MULTI_PARAMS, fitter.getA().length);
            boolean valid = true;
            for (int i = 0; i < SINE_MULTI_PARAMS; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }

            assertNotNull(fitter.getCovar());
            assertEquals(SINE_MULTI_PARAMS, fitter.getCovar().getRows());
            assertEquals(SINE_MULTI_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final double chiSqrDegreesOfFreedom = npoints - SINE_MULTI_PARAMS;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            final double standardDeviationAmplitude = Math.sqrt(
                    cov.getElementAt(0, 0));
            final double standardDeviationFreqx = Math.sqrt(
                    cov.getElementAt(1, 1));
            final double standardDeviationFreqy = Math.sqrt(
                    cov.getElementAt(2, 2));
            final double standardDeviationPhasex = Math.sqrt(
                    cov.getElementAt(3, 3));
            final double standardDeviationPhasey = Math.sqrt(
                    cov.getElementAt(4, 4));

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
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS,
                    MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_MULTI_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            final double[] varianceParams = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            final Matrix y = new Matrix(npoints, 1);
            final Matrix x = new Matrix(npoints, 2);
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final double[] point = new double[2];
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(varianceParams);
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

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {

                    final double[] derivatives = new double[params.length];

                    private final GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(final double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

                    @Override
                    public void evaluate(
                            final double[] params, final double[] y, final Matrix jacobian) {
                        try {
                            gradientEstimator.gradient(params, derivatives);
                            y[0] = evaluateParams(point, params);
                            jacobian.fromArray(derivatives);
                        } catch (final Exception ignore) {
                        }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }

                    double evaluateParams(final double[] point, final double[] params) {
                        double y = 0.0;
                        for (int k = 0; k < numgaussians; k++) {
                            final double b = params[k * GAUSS_MULTI_PARAMS];
                            final double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                            final double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                            final double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                            final double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                            y += b * Math.exp(-(Math.pow(point[0] -
                                    ex, 2.0) + Math.pow(point[1] - ey, 2.0)) /
                                    Math.pow(gx * gy, 2.0));
                        }

                        return y;
                    }

                }, params, covariance, dist);

                double yi = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    final double b = params[k * GAUSS_MULTI_PARAMS];
                    final double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    final double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    final double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    final double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    yi += b * Math.exp(-(Math.pow(x.getElementAt(i, 0) - ex, 2.0) +
                            Math.pow(x.getElementAt(i, 1) - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                }

                double value = dist.getMean()[0];
                assertEquals(value, yi, SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                value += error;
                y.setElementAt(i, 0, value);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private final JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(
                                                    final double[] params, final double[] result) {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 1;
                                            }
                                        });

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
                            final double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian)
                                throws EvaluationException {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }

                        void evaluateParams(final double[] point, final double[] params,
                                            final double[] result) {
                            double y = 0.0;
                            for (int k = 0; k < numgaussians; k++) {
                                final double b = params[k * GAUSS_MULTI_PARAMS];
                                final double ex = params[k * GAUSS_MULTI_PARAMS + 1];
                                final double ey = params[k * GAUSS_MULTI_PARAMS + 2];
                                final double gx = params[k * GAUSS_MULTI_PARAMS + 3];
                                final double gy = params[k * GAUSS_MULTI_PARAMS + 4];
                                y += b * Math.exp(-(Math.pow(point[0] - ex, 2.0) +
                                        Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                            }

                            result[0] = y;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (FittingException e) {
                continue;
            }


            // check correctness
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

            final double chiSqrDegreesOfFreedom = npoints - numParams;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            for (int i = 0; i < numParams; i++) {
                final double standardDeviation = Math.sqrt(cov.getElementAt(i, i));

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
    public void testFitMultiVariateGaussianAndSineCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
            final int numParams = SINE_MULTI_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final Matrix y = new Matrix(npoints, 2);
            final Matrix x = new Matrix(npoints, 2);
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                final double xi0 = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                final double xi1 = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi0);
                x.setElementAt(i, 1, xi1);

                double valueGaussian;
                double valueSine;

                final double b = params[0];
                final double ex = params[1];
                final double ey = params[2];
                final double gx = params[3];
                final double gy = params[4];

                error = errorRandomizer.nextDouble();

                valueGaussian = b * Math.exp(-(Math.pow(x.getElementAt(i, 0) - ex, 2.0) +
                        Math.pow(x.getElementAt(i, 1) - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                valueGaussian += error;

                final double amplitude = params[0];
                final double freqx = params[1];
                final double freqy = params[2];
                final double phasex = params[3];
                final double phasey = params[4];

                valueSine = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                        Math.sin(freqy * x.getElementAt(i, 1) + phasey);
                valueSine += error;

                y.setElementAt(i, 0, valueGaussian);
                y.setElementAt(i, 1, valueSine);
            }

            final LevenbergMarquardtMultiVariateFunctionEvaluator evaluator =
                    new LevenbergMarquardtMultiVariateFunctionEvaluator() {

                        private double[] point;

                        private final JacobianEstimator jacobianEstimator =
                                new JacobianEstimator(
                                        new MultiVariateFunctionEvaluatorListener() {

                                            @Override
                                            public void evaluate(
                                                    final double[] params, final double[] result) {
                                                evaluateParams(point, params, result);
                                            }

                                            @Override
                                            public int getNumberOfVariables() {
                                                return 2;
                                            }
                                        });

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
                            final double[] initParams = new double[numParams];
                            double error;
                            for (int i = 0; i < numParams; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public void evaluate(
                                final int i, final double[] point, final double[] result,
                                final double[] params, final Matrix jacobian)
                                throws EvaluationException {
                            this.point = point;
                            evaluateParams(point, params, result);
                            jacobianEstimator.jacobian(params, jacobian);
                        }

                        void evaluateParams(final double[] point, final double[] params,
                                            final double[] result) {
                            final double b = params[0];
                            final double ex = params[1];
                            final double ey = params[2];
                            final double gx = params[3];
                            final double gy = params[4];

                            final double amplitude = params[0];
                            final double freqx = params[1];
                            final double freqy = params[2];
                            final double phasex = params[3];
                            final double phasey = params[4];

                            final double valueGaussian = b * Math.exp(-(Math.pow(point[0] - ex, 2.0) +
                                    Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));

                            final double valueSine = amplitude * Math.sin(freqx * point[0] + phasex) *
                                    Math.sin(freqy * point[1] + phasey);

                            result[0] = valueGaussian;
                            result[1] = valueSine;
                        }
                    };

            final LevenbergMarquardtMultiVariateFitter fitter =
                    new LevenbergMarquardtMultiVariateFitter(evaluator, x, y, sigma);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(0.0, fitter.getChisq(), 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            try {
                fitter.fit();
            } catch (final FittingException e) {
                continue;
            }

            // check correctness
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

            final double chiSqrDegreesOfFreedom = npoints - numParams;
            final double chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final double p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final double q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, "chi sqr: " + chiSqr +
                    ", probability smaller chi sqr: " + p * 100.0 +
                    "%, quality: " + q * 100.0 + "%");

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final int nVars = evaluator.getNumberOfVariables();
            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix jacobianTrans = new Matrix(params.length, nVars);
            final Matrix jacobian = new Matrix(nVars, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] point = new double[evaluator.getNumberOfDimensions()];
            final double[] estimatedParams = fitter.getA();
            final double[] result = new double[nVars];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i,
                        evaluator.getNumberOfDimensions() - 1,
                        point);
                evaluator.evaluate(i, point, result, estimatedParams,
                        jacobian);

                jacobian.transpose(jacobianTrans);

                jacobianTrans.multiply(jacobian, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigma * sigma);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                for (int j = 0; j < nVars; j++) {
                    mse2 += Math.pow(y.getElementAt(i, j) - result[j], 2.0) / (npoints - params.length);
                }
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final Matrix cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            for (int i = 0; i < numParams; i++) {
                final double standardDeviation = Math.sqrt(cov.getElementAt(i, i));

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
