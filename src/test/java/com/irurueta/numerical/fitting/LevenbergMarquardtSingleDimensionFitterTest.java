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
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.GradientEstimator;
import com.irurueta.numerical.MultiDimensionFunctionEvaluatorListener;
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

@SuppressWarnings({"Duplicates", "DuplicateExpressions"})
public class LevenbergMarquardtSingleDimensionFitterTest {

    private static final Logger LOGGER = Logger.getLogger(
            LevenbergMarquardtSingleDimensionFitterTest.class.getName());

    private static final int MIN_POINTS = 500;
    private static final int MAX_POINTS = 1000;

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;

    private static final double MIN_SIGMA_VALUE = 1e-4;
    private static final double MAX_SIGMA_VALUE = 1e-3;

    private static final double ABSOLUTE_ERROR = 1e-1;
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

    private static final int SINE_PARAMS = 3;
    private static final double MIN_SINE_AMPLITUDE = 0.5;
    private static final double MAX_SINE_AMPLITUDE = 10.0;
    private static final double MIN_SINE_FREQ = 0.5;
    private static final double MAX_SINE_FREQ = 100.0;
    private static final double MIN_SINE_PHASE = -Math.PI;
    private static final double MAX_SINE_PHASE = Math.PI;

    private static final int GAUSS_PARAMS = 3;
    private static final int MIN_GAUSSIANS = 1;
    private static final int MAX_GAUSSIANS = 3;

    private static final int TIMES = 50;
    private static final int N_SAMPLES = 1000000;

    @Test
    public void testConstructor() throws FittingException {
        LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // test constructor with input data
        final double[] x = new double[2];
        final double[] y = new double[2];
        final double[] sig = new double[2];

        fitter = new LevenbergMarquardtSingleDimensionFitter(x, y, sig);

        // check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        final double[] shortX = new double[1];
        final double[] shortSig = new double[1];

        fitter = null;
        try {
            fitter = new LevenbergMarquardtSingleDimensionFitter(shortX, y,
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            //noinspection SuspiciousNameCombination
            fitter = new LevenbergMarquardtSingleDimensionFitter(x, shortX,
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtSingleDimensionFitter(x, y,
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtSingleDimensionFitter(x, y, 1.0);

        // check default value
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
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtSingleDimensionFitter(shortX, y,
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            //noinspection SuspiciousNameCombination
            fitter = new LevenbergMarquardtSingleDimensionFitter(x, shortX,
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator
        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createInitialParametersArray() {
                        return new double[GAUSS_PARAMS];
                    }

                    @Override
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        return 0.0;
                    }
                };

        fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator);

        // check default value
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
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_PARAMS);
        assertTrue(fitter.isCovarianceAdjusted());

        // test constructor with input data
        fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                sig);

        // check default value
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
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_PARAMS);
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator,
                    shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            //noinspection SuspiciousNameCombination
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x,
                    shortX, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x,
                    y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                1.0);

        // check default value
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_PARAMS);
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getAlpha().getColumns(), GAUSS_PARAMS);
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator,
                    shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            //noinspection SuspiciousNameCombination
            fitter = new LevenbergMarquardtSingleDimensionFitter(evaluator, x,
                    shortX, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);
    }

    @Test
    public void testGetSetNdone() {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default values
        assertEquals(fitter.getNdone(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE);

        // set new value
        fitter.setNdone(5);

        // check correctness
        assertEquals(fitter.getNdone(), 5);

        // force IllegalArgumentException
        try {
            fitter.setNdone(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetItmax() {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default values
        assertEquals(fitter.getItmax(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX);

        // set new value
        fitter.setItmax(10);

        // check correctness
        assertEquals(fitter.getItmax(), 10);

        // force IllegalArgumentException
        try {
            fitter.setItmax(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetTol() {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default values
        assertEquals(fitter.getTol(),
                LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, 0.0);

        // set new value
        fitter.setTol(1e-1);

        // check correctness
        assertEquals(fitter.getTol(), 1e-1, 0.0);

        // force IllegalArgumentException
        try {
            fitter.setTol(0.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetFunctionEvaluator() throws FittingException {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default value
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createInitialParametersArray() {
                        return new double[GAUSS_PARAMS];
                    }

                    @Override
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        return 0.0;
                    }
                };

        // set new value
        fitter.setFunctionEvaluator(evaluator);

        // check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, GAUSS_PARAMS);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), GAUSS_PARAMS);
        assertEquals(fitter.getCovar().getColumns(), GAUSS_PARAMS);
    }

    @Test
    public void testGetSetInputData() {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final double[] x = new double[2];
        final double[] y = new double[2];
        final double[] sig = new double[2];

        // set input data
        fitter.setInputData(x, y, sig);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);

        // Force IllegalArgumentException
        final double[] wrong = new double[1];

        try {
            fitter.setInputData(wrong, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, wrong, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, y, wrong);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetInputDataWithConstantSigma() {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final double[] x = new double[2];
        final double[] y = new double[2];


        // set input data
        fitter.setInputData(x, y, 1.0);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }

        // Force IllegalArgumentException
        final double[] wrong = new double[1];

        try {
            fitter.setInputData(wrong, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, wrong, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsReady() throws FittingException {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final double[] x = new double[2];
        final double[] y = new double[2];
        final double[] sig = new double[2];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createInitialParametersArray() {
                        return new double[GAUSS_PARAMS];
                    }

                    @Override
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        return 0.0;
                    }
                };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());
    }

    @Test
    public void testIsSetCovarianceAdjusted() {
        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter();

        // check default value
        assertTrue(fitter.isCovarianceAdjusted());

        // set new value
        fitter.setCovarianceAdjusted(false);

        // check
        assertFalse(fitter.isCovarianceAdjusted());
    }

    @Test
    public void testFitConstant() throws FittingException, NotReadyException,
            MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[CONSTANT_PARAMS];
        params[0] = constant;

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = constant;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
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
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        final double constant = params[0];

                        // derivative of evaluated function respect constant parameter
                        derivatives[0] = 1.0;

                        // evaluated function
                        return constant;
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, CONSTANT_PARAMS);
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
    public void testFitLine1() throws FittingException, NotReadyException,
            MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[LINE1_PARAMS];
        params[0] = a;

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = a * x[i];
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
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
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        final double a = params[0];
                        derivatives[0] = a;
                        return a * point;
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, LINE1_PARAMS);
        for (int i = 0; i < LINE1_PARAMS; i++) {
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
    public void testFitLine2() throws FittingException, NotReadyException,
            MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
        final double b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[LINE2_PARAMS];
        params[0] = a;
        params[1] = b;

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = a * x[i] + b;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
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
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        double a = params[0];
                        double b = params[1];

                        // derivatives of function f(x) respect parameters a and b
                        derivatives[0] = point;
                        derivatives[1] = 1.0;

                        // evaluated function f(x) = a * x + b
                        return a * point + b;
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, LINE2_PARAMS);
        for (int i = 0; i < LINE2_PARAMS; i++) {
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
    public void testFitSine() throws FittingException, NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = amplitude * Math.sin(freq * x[i] + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createInitialParametersArray() {
                        final double[] initParams = new double[SINE_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        final double amplitude = params[0];
                        final double freq = params[1];
                        final double phase = params[2];
                        final double y = amplitude * Math.sin(freq * point + phase);

                        // derivative respect amplitude
                        derivatives[0] = Math.sin(freq * point + phase);
                        // derivative respect frequency
                        derivatives[1] = amplitude * Math.cos(freq * point + phase) * point;
                        // derivative respect phase
                        derivatives[2] = amplitude * Math.cos(freq * point + phase);

                        return y;
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);
        fitter.setCovarianceAdjusted(false);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, SINE_PARAMS);
        for (int i = 0; i < SINE_PARAMS; i++) {
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
    public void testFitGaussian() throws FittingException, NotReadyException,
            MaxIterationsExceededException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                y[i] = 0.0;
                for (int k = 0; k < numgaussians; k++) {
                    final double b = params[k * GAUSS_PARAMS];
                    final double e = params[k * GAUSS_PARAMS + 1];
                    final double g = params[k * GAUSS_PARAMS + 2];
                    y[i] += b * Math.exp(-Math.pow((x[i] - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

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
                        public double evaluate(
                                final int pos, final double point, final double[] params,
                                final double[] derivatives) {
                            int i;
                            final int na = params.length;
                            double fac;
                            double ex;
                            double arg;
                            double y = 0.0;
                            for (i = 0; i < na - 1; i += 3) {
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

            final LevenbergMarquardtSingleDimensionFitter fitter =
                    new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                            1.0);
            fitter.setCovarianceAdjusted(false);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
            assertFalse(fitter.isResultAvailable());
            assertTrue(fitter.isReady());

            // fit
            fitter.fit();

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

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitSineWithHoldAndFree() throws FittingException,
            NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = amplitude * Math.sin(freq * x[i] + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createInitialParametersArray() {
                        final double[] initParams = new double[SINE_PARAMS];
                        initParams[0] = params[0];
                        double error;
                        for (int i = 1; i < SINE_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) {
                        final double amplitude = params[0];
                        final double freq = params[1];
                        final double phase = params[2];
                        final double y = amplitude * Math.sin(freq * point + phase);

                        // derivative respect amplitude
                        derivatives[0] = Math.sin(freq * point + phase);
                        // derivative respect frequency
                        derivatives[1] = amplitude * Math.cos(freq * point + phase) * point;
                        // derivative respect phase
                        derivatives[2] = amplitude * Math.cos(freq * point + phase);

                        return y;
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // hold first parameter
        fitter.hold(0, params[0]);

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, SINE_PARAMS);
        // first parameter is hold and matches exactly
        assertEquals(fitter.getA()[0], params[0], 0.0);
        for (int i = 0; i < SINE_PARAMS; i++) {
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
        assertEquals(fitter.getA().length, SINE_PARAMS);
        for (int i = 0; i < SINE_PARAMS; i++) {
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
    public void testFitGaussianWithGradientEstimator() throws FittingException,
            NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
        final int numParams = numgaussians * GAUSS_PARAMS;

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                    MAX_RANDOM_VALUE);
        }

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = 0.0;
            for (int k = 0; k < numgaussians; k++) {
                final double b = params[k * GAUSS_PARAMS];
                final double e = params[k * GAUSS_PARAMS + 1];
                final double g = params[k * GAUSS_PARAMS + 2];
                y[i] += b * Math.exp(-Math.pow((x[i] - e) / g, 2.0));
            }
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    private double point;

                    private final GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(final double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

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
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) throws EvaluationException {
                        this.point = point;
                        final double y = evaluateParams(point, params);
                        gradientEstimator.gradient(params, derivatives);

                        return y;
                    }

                    double evaluateParams(final double point, final double[] params) {
                        int i;
                        final int na = params.length;
                        double ex;
                        double arg;
                        double y = 0.0;
                        for (i = 0; i < na - 1; i += 3) {
                            arg = (point - params[i + 1]) / params[i + 2];
                            ex = Math.exp(-Math.pow(arg, 2.0));
                            y += params[i] * ex;
                        }

                        return y;
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, numParams);
        for (int i = 0; i < numParams; i++) {
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
    public void testFitSineWithGradientEstimator() throws FittingException,
            NotReadyException, MaxIterationsExceededException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());

        final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE,
                MAX_SINE_AMPLITUDE);
        final double freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final double phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final double[] params = new double[SINE_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final double[] y = new double[npoints];
        final double[] x = new double[npoints];
        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = amplitude * Math.sin(freq * x[i] + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                    private double point;

                    private final GradientEstimator gradientEstimator =
                            new GradientEstimator(
                                    new MultiDimensionFunctionEvaluatorListener() {

                                        @Override
                                        public double evaluate(final double[] params) {
                                            return evaluateParams(point, params);
                                        }
                                    });

                    @Override
                    public double[] createInitialParametersArray() {
                        final double[] initParams = new double[SINE_PARAMS];
                        double error;
                        for (int i = 0; i < SINE_PARAMS; i++) {
                            error = errorRandomizer.nextDouble();
                            initParams[i] = params[i] + error;
                        }
                        return initParams;
                    }

                    @Override
                    public double evaluate(
                            final int i, final double point, final double[] params,
                            final double[] derivatives) throws EvaluationException {
                        this.point = point;
                        final double y = evaluateParams(point, params);
                        gradientEstimator.gradient(params, derivatives);

                        return y;
                    }

                    double evaluateParams(final double point, final double[] params) {
                        final double amplitude = params[0];
                        final double freq = params[1];
                        final double phase = params[2];
                        return amplitude * Math.sin(freq * point + phase);
                    }
                };

        final LevenbergMarquardtSingleDimensionFitter fitter =
                new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                        1.0);

        // check default values
        assertNotNull(fitter.getA());
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getChisq(), 0.0, 0.0);
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());

        // fit
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, SINE_PARAMS);
        for (int i = 0; i < SINE_PARAMS; i++) {
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
    public void testFitConstantCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[CONSTANT_PARAMS];
            params[0] = constant;

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            final NormalDist dist = new NormalDist();
            double error;
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

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
                y[i] = dist.getMean();
                assertEquals(y[i], constant, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
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
                        public double evaluate(
                                final int i, final double point, final double[] params,
                                final double[] derivatives) {
                            final double constant = params[0];

                            // derivative respect parameter
                            derivatives[0] = 1.0;

                            return constant;
                        }
                    };

            final LevenbergMarquardtSingleDimensionFitter fitter =
                    new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                            sigma);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
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
            assertEquals(fitter.getA().length, CONSTANT_PARAMS);
            for (int i = 0; i < CONSTANT_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], SMALL_ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), CONSTANT_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), CONSTANT_PARAMS);
            assertTrue(fitter.getChisq() > 0.0);

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

            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix tmp1 = new Matrix(params.length, 1);
            final Matrix tmp2 = new Matrix(1, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                final double yi = evaluator.evaluate(i, x[i], estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
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
    public void testFitLine1Covariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = N_SAMPLES;
            final double a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

            // this is the standard deviation that we expect on estimated parameter a
            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[LINE1_PARAMS];
            params[0] = a;

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            final NormalDist dist = new NormalDist();
            double error;
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final double xi = x[i];

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
                y[i] = dist.getMean();
                assertEquals(y[i], a * x[i], SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
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
                        public double evaluate(
                                final int i, final double point, final double[] params,
                                final double[] derivatives) {
                            final double a = params[0];

                            // derivative of function f(x) respect parameter a
                            derivatives[0] = point;

                            return a * point;
                        }
                    };

            final LevenbergMarquardtSingleDimensionFitter fitter =
                    new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
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
            assertEquals(fitter.getA().length, LINE1_PARAMS);
            for (int i = 0; i < LINE1_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), LINE1_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), LINE1_PARAMS);
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

            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix tmp1 = new Matrix(params.length, 1);
            final Matrix tmp2 = new Matrix(1, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                final double yi = evaluator.evaluate(i, x[i], estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
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
    public void testFitLine2Covariance() throws Throwable {
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

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(new double[]{varianceA, varianceB});
            double error;
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final double xi = x[i];

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
                y[i] = dist.getMean()[0];
                assertEquals(y[i], a * x[i] + b, SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtSingleDimensionFunctionEvaluator() {
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
                        public double evaluate(
                                final int i, final double point, final double[] params,
                                final double[] derivatives) {
                            final double a = params[0];
                            final double b = params[1];

                            // derivatives of function f(x) respect parameters a and b
                            derivatives[0] = point;
                            derivatives[1] = 1.0;

                            // evaluated function f(x) = a * x + b
                            return a * point + b;
                        }
                    };

            final LevenbergMarquardtSingleDimensionFitter fitter =
                    new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
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
            assertEquals(fitter.getA().length, LINE2_PARAMS);
            for (int i = 0; i < LINE2_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), LINE2_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), LINE2_PARAMS);
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

            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix tmp1 = new Matrix(params.length, 1);
            final Matrix tmp2 = new Matrix(1, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                final double yi = evaluator.evaluate(i, x[i], estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
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
    public void testFitSineCovariance() throws Throwable {
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

            final double[] params = new double[SINE_PARAMS];
            params[0] = amplitude;
            params[1] = freq;
            params[2] = phase;

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(new double[]{
                    varianceAmplitude, varianceFreq, variancePhase});
            double error;
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final double xi = x[i];

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
                y[i] = dist.getMean()[0];
                assertEquals(y[i], amplitude * Math.sin(freq * x[i] + phase),
                        SMALL_ABSOLUTE_ERROR);

                assertEquals(dist.getCovariance().getRows(), 1);
                assertEquals(dist.getCovariance().getColumns(), 1);

                sigmas[i] = Math.sqrt(dist.getCovariance().
                        getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

                        @Override
                        public double[] createInitialParametersArray() {
                            final double[] initParams = new double[SINE_PARAMS];
                            double error;
                            for (int i = 0; i < SINE_PARAMS; i++) {
                                error = errorRandomizer.nextDouble();
                                initParams[i] = params[i] + error;
                            }
                            return initParams;
                        }

                        @Override
                        public double evaluate(
                                final int i, final double point, final double[] params,
                                final double[] derivatives) {
                            final double amplitude = params[0];
                            final double freq = params[1];
                            final double phase = params[2];
                            final double y = amplitude * Math.sin(freq * point + phase);

                            // derivative respect amplitude
                            derivatives[0] = Math.sin(freq * point + phase);
                            // derivative respect frequency
                            derivatives[1] = amplitude * Math.cos(freq * point + phase) * point;
                            // derivative respect phase
                            derivatives[2] = amplitude * Math.cos(freq * point + phase);

                            return y;
                        }
                    };

            final LevenbergMarquardtSingleDimensionFitter fitter =
                    new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
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
            assertEquals(fitter.getA().length, SINE_PARAMS);
            boolean valid = true;
            for (int i = 0; i < SINE_PARAMS; i++) {
                if (Math.abs(fitter.getA()[i] - params[i]) > ABSOLUTE_ERROR) {
                    valid = false;
                    break;
                }
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getCovar().getRows(), SINE_PARAMS);
            assertEquals(fitter.getCovar().getColumns(), SINE_PARAMS);
            assertTrue(fitter.getChisq() > 0);

            final double chiSqrDegreesOfFreedom = npoints - SINE_PARAMS;
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

            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix tmp1 = new Matrix(params.length, 1);
            final Matrix tmp2 = new Matrix(1, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                final double yi = evaluator.evaluate(i, x[i], estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
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
    public void testFitGaussianCovariance() throws Throwable {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final int numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final int numParams = numgaussians * GAUSS_PARAMS;

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[numParams];
            final double[] varianceParams = new double[numParams];
            for (int i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final double[] sigmas = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, 1.0);
            final MultivariateNormalDist dist = new MultivariateNormalDist();
            final Matrix covariance = Matrix.diagonal(varianceParams);
            double error;
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final double xi = x[i];

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(
                            final double[] params, final double[] y, final Matrix jacobian) {
                        y[0] = 0.0;
                        for (int k = 0; k < numgaussians; k++) {
                            final double b = params[k * GAUSS_PARAMS];
                            final double e = params[k * GAUSS_PARAMS + 1];
                            final double g = params[k * GAUSS_PARAMS + 2];
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
                    final double b = params[k * GAUSS_PARAMS];
                    final double e = params[k * GAUSS_PARAMS + 1];
                    final double g = params[k * GAUSS_PARAMS + 2];
                    yi += b * Math.exp(-Math.pow((x[i] - e) / g, 2.0));
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

            final LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator =
                    new LevenbergMarquardtSingleDimensionFunctionEvaluator() {

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
                        public double evaluate(
                                final int pos, final double point, final double[] params,
                                final double[] derivatives) {
                            int i;
                            final int na = params.length;
                            double fac;
                            double ex;
                            double arg;
                            double y = 0.0;
                            for (i = 0; i < na - 1; i += 3) {
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

            final LevenbergMarquardtSingleDimensionFitter fitter =
                    new LevenbergMarquardtSingleDimensionFitter(evaluator, x, y,
                            sigmas);
            fitter.setCovarianceAdjusted(true);

            // check default values
            assertNotNull(fitter.getA());
            assertNotNull(fitter.getCovar());
            assertEquals(fitter.getChisq(), 0.0, 0.0);
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

            final double mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final Matrix invCov = new Matrix(params.length, params.length);
            final Matrix tmp1 = new Matrix(params.length, 1);
            final Matrix tmp2 = new Matrix(1, params.length);
            final Matrix tmpInvCov = new Matrix(params.length, params.length);
            final double[] estimatedParams = fitter.getA();
            final double[] derivatives = new double[params.length];
            double mse2 = 0.0;
            for (int i = 0; i < npoints; i++) {
                final double yi = evaluator.evaluate(i, x[i], estimatedParams,
                        derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final double w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
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
