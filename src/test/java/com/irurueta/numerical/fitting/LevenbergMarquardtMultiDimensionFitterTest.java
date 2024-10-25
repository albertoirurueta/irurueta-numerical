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
import com.irurueta.statistics.ChiSqDist;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.MaxIterationsExceededException;
import com.irurueta.statistics.MultivariateNormalDist;
import com.irurueta.statistics.NormalDist;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.jupiter.api.Assertions.*;

class LevenbergMarquardtMultiDimensionFitterTest {

    private static final Logger LOGGER = Logger.getLogger(LevenbergMarquardtMultiDimensionFitterTest.class.getName());

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

    public static final int TIMES = 10;
    private static final int N_SAMPLES = 1000000;

    @Test
    void testConstructor() throws FittingException, WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtSingleDimensionFitter.DEFAULT_NDONE, fitter.getNdone());
        assertEquals(LevenbergMarquardtSingleDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());
        assertEquals(LevenbergMarquardtSingleDimensionFitter.DEFAULT_TOL, fitter.getTol(), 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // test constructor with input data
        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];
        final var sig = new double[nPoints];

        fitter = new LevenbergMarquardtMultiDimensionFitter(x, y, sig);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE, fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, fitter.getTol(), 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        final var shortX = new Matrix(nPoints - 1, 2);
        final var shortY = new double[nPoints - 1];
        final var shortSig = new double[nPoints - 1];

        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(shortX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(x, shortY, sig));
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(x, y, shortSig));

        // test constructor with input data (constant sigma)
        fitter = new LevenbergMarquardtMultiDimensionFitter(x, y, 1.0);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertFalse(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE, fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, fitter.getTol(), 0.0);
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getAlpha());
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class,
                () -> new LevenbergMarquardtMultiDimensionFitter(shortX, y, 1.0));
        assertThrows(IllegalArgumentException.class,
                () -> new LevenbergMarquardtMultiDimensionFitter(x, shortY, 1.0));

        // test constructor with evaluator
        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[nPoints];
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                return 0.0;
            }
        };

        fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator);

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
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE, fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, fitter.getTol(), 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        // test constructor with evaluator and input data
        fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sig);

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
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE, fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, fitter.getTol(), 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(evaluator,
                shortX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(evaluator, x,
                shortY, sig));
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y,
                shortSig));

        // test constructor with evaluator and input data (constant sigma)
        fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertTrue(fitter.isReady());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertNotNull(fitter.getA());
        assertEquals(fitter.getA().length, nPoints);
        assertNotNull(fitter.getCovar());
        assertEquals(fitter.getCovar().getRows(), nPoints);
        assertEquals(fitter.getCovar().getColumns(), nPoints);
        assertEquals(0.0, fitter.getChisq(), 0.0);
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE, fitter.getNdone());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL, fitter.getTol(), 0.0);
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNotNull(fitter.getAlpha());
        assertEquals(fitter.getAlpha().getRows(), nPoints);
        assertEquals(fitter.getAlpha().getColumns(), nPoints);
        assertTrue(fitter.isCovarianceAdjusted());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(evaluator,
                shortX, y, 1.0));
        assertThrows(IllegalArgumentException.class, () -> new LevenbergMarquardtMultiDimensionFitter(evaluator, x,
                shortY, 1.0));
    }

    @Test
    void testGetSetNdone() {
        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_NDONE, fitter.getNdone());

        // new value
        fitter.setNdone(5);

        // check correctness
        assertEquals(5, fitter.getNdone());

        // force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> fitter.setNdone(0));
    }

    @Test
    void testGetSetItmax() {
        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_ITMAX, fitter.getItmax());

        // new value
        fitter.setItmax(10);

        // check correctness
        assertEquals(10, fitter.getItmax());

        // force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> fitter.setItmax(0));
    }

    @Test
    void testGetSetTol() {
        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertEquals(LevenbergMarquardtMultiDimensionFitter.DEFAULT_TOL,
                fitter.getTol(), 0.0);

        // new value
        fitter.setTol(1e-1);

        //check correctness
        assertEquals(1e-1, fitter.getTol(), 0.0);

        // force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> fitter.setTol(0.0));
    }

    @Test
    void testGetSetFunctionEvaluator() throws FittingException {
        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getA());
        assertNull(fitter.getCovar());
        assertNull(fitter.getAlpha());

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_UNI_PARAMS];
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                return 0.0;
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
    void testGetSetInputData() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];
        final var sig = new double[nPoints];

        // set input data
        fitter.setInputData(x, y, sig);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);

        // Force IllegalArgumentException
        final var shortX = new Matrix(nPoints - 1, 2);
        final var shortY = new double[nPoints - 1];
        final var shortSig = new double[nPoints - 1];

        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(shortX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, shortY, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, y, shortSig));
    }

    @Test
    void testGetSetInputDataWithConstantSigma() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];

        // set input data
        fitter.setInputData(x, y, 1.0);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }

        // Force IllegalArgumentException
        final var shortX = new Matrix(nPoints - 1, 2);
        final var shortY = new double[nPoints - 1];

        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(shortX, y, 1.0));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, shortY, 1.0));
    }

    @Test
    void testIsReady() throws WrongSizeException, FittingException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];
        final var sig = new double[nPoints];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                return new double[GAUSS_UNI_PARAMS];
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                return 0.0;
            }
        };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());
    }

    @Test
    void testIsSetCovarianceAdjusted() {
        final var fitter = new LevenbergMarquardtMultiDimensionFitter();

        // check default value
        assertTrue(fitter.isCovarianceAdjusted());

        // set new value
        fitter.setCovarianceAdjusted(false);

        // check
        assertFalse(fitter.isCovarianceAdjusted());
    }

    @Test
    void testFitUnidimensionalConstant() throws WrongSizeException, FittingException, NotReadyException,
            MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[CONSTANT_PARAMS];
        params[0] = constant;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 1);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (var i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            y[i] = constant;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[CONSTANT_PARAMS];
                double error;
                for (var i = 0; i < CONSTANT_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                final var constant = params[0];

                // derivative of evaluated function respect constant parameter
                derivatives[0] = 1.0;

                // evaluated function
                return constant;
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
        for (var i = 0; i < CONSTANT_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final var chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalLine1() throws WrongSizeException, FittingException, NotReadyException,
            MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[LINE1_PARAMS];
        params[0] = a;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 1);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (var i = 0; i < npoints; i++) {
            final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            y[i] = a * xi;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[LINE1_PARAMS];
                double error;
                for (var i = 0; i < LINE1_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                final var a = params[0];

                // derivative of evaluated function respect constant parameter
                derivatives[0] = a;

                // evaluated function
                return a * point[0];
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
        for (var i = 0; i < LINE1_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final var chiSqrDegreesOfFreedom = npoints - LINE1_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalLine2() throws WrongSizeException, FittingException, NotReadyException,
            MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
        final var b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[LINE2_PARAMS];
        params[0] = a;
        params[1] = b;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 1);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (var i = 0; i < npoints; i++) {
            final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            y[i] = a * xi + b;
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[LINE2_PARAMS];
                double error;
                for (var i = 0; i < LINE2_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                final var a = params[0];
                final var b = params[1];

                // derivatives of function f(x) respect parameters a and b
                derivatives[0] = point[0];
                derivatives[1] = 1.0;

                // evaluated function f(x) = a * x + b
                return a * point[0] + b;
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
        for (var i = 0; i < LINE2_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final var chiSqrDegreesOfFreedom = npoints - LINE2_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalSine() throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
        final var freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final var phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 1);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (int i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            y[i] = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[SINE_UNI_PARAMS];
                double error;
                for (var i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                final var amplitude = params[0];
                final var freq = params[1];
                final var phase = params[2];
                final var y = amplitude * Math.sin(freq * point[0] + phase);

                // derivative respect amplitude
                derivatives[0] = Math.sin(freq * point[0] + phase);
                // derivative respect frequency
                derivatives[1] = amplitude * Math.cos(freq * point[0] + phase) * point[0];
                // derivative respect phase
                derivatives[2] = amplitude * Math.cos(freq * point[0] + phase);

                return y;
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
        for (var i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final var chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f" + "%%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalGaussian() throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final var numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final var numParams = numgaussians * GAUSS_UNI_PARAMS;

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[numParams];
            for (var i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            }

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);
                y[i] = 0.0;
                for (var k = 0; k < numgaussians; k++) {
                    final var b = params[k * GAUSS_UNI_PARAMS];
                    final var e = params[k * GAUSS_UNI_PARAMS + 1];
                    final var g = params[k * GAUSS_UNI_PARAMS + 2];
                    y[i] += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[numParams];
                    double error;
                    for (var i = 0; i < numParams; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int pos, final double[] point, final double[] params, final double[] derivatives) {
                    int i;
                    final var na = params.length;
                    double fac;
                    double ex;
                    double arg;
                    var y = 0.0;
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

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
            var valid = true;
            for (var i = 0; i < numParams; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - numParams;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            numValid++;
            break;
        }
        assertTrue(numValid > 0);
    }

    @Test
    void testFitUnidimensionalSineWithHoldAndFree() throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
        final var freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final var phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 1);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (var i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            y[i] = amplitude * Math.sin(freq * x.getElementAt(i, 0) + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[SINE_UNI_PARAMS];
                double error;
                for (var i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives) {
                final var amplitude = params[0];
                final var freq = params[1];
                final var phase = params[2];
                final var y = amplitude * Math.sin(freq * point[0] + phase);

                // derivative respect amplitude
                derivatives[0] = Math.sin(freq * point[0] + phase);
                // derivative respect frequency
                derivatives[1] = amplitude * Math.cos(freq * point[0] + phase) * point[0];
                // derivative respect phase
                derivatives[2] = amplitude * Math.cos(freq * point[0] + phase);

                return y;
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
        for (var i = 0; i < SINE_UNI_PARAMS; i++) {
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
        for (var i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final var chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalSineWithGradientEstimator() throws FittingException, NotReadyException,
            WrongSizeException, MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
        final var freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final var phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freq;
        params[2] = phase;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 1);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (var i = 0; i < npoints; i++) {
            final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            x.setElementAt(i, 0, xi);
            y[i] = amplitude * Math.sin(freq * xi + phase);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            private double[] point;

            private final GradientEstimator gradientEstimator =
                    new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                        @Override
                        public double evaluate(final double[] params) {
                            return evaluateParams(point, params);
                        }
                    });

            @Override
            public int getNumberOfDimensions() {
                return 1;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[SINE_UNI_PARAMS];
                double error;
                for (var i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives)
                    throws EvaluationException {
                this.point = point;
                var y = evaluateParams(point, params);
                gradientEstimator.gradient(params, derivatives);

                return y;
            }

            double evaluateParams(final double[] point, final double[] params) {
                final var amplitude = params[0];
                final var freq = params[1];
                final var phase = params[2];
                return amplitude * Math.sin(freq * point[0] + phase);
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);

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
        for (var i = 0; i < SINE_UNI_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);

        final var chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalGaussianWithGradientEstimator() throws FittingException, NotReadyException,
            WrongSizeException, MaxIterationsExceededException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final var numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final var numParams = numgaussians * GAUSS_UNI_PARAMS;

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[numParams];
            for (var i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            }

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            double error;
            for (var i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                y[i] = 0.0;
                for (var k = 0; k < numgaussians; k++) {
                    final var b = params[k * GAUSS_UNI_PARAMS];
                    final var e = params[k * GAUSS_UNI_PARAMS + 1];
                    final var g = params[k * GAUSS_UNI_PARAMS + 2];
                    y[i] += b * Math.exp(-Math.pow((x.getElementAt(i, 0) - e) / g, 2.0));
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                private double[] point;

                private final GradientEstimator gradientEstimator =
                        new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                            @Override
                            public double evaluate(final double[] params) {
                                return evaluateParams(point, params);
                            }
                        });

                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[numParams];
                    double error;
                    for (var i = 0; i < numParams; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives)
                        throws EvaluationException {
                    this.point = point;
                    final var y = evaluateParams(point, params);
                    gradientEstimator.gradient(params, derivatives);

                    return y;
                }

                double evaluateParams(final double[] point, final double[] params) {
                    int i;
                    final var na = params.length;
                    double ex;
                    double arg;
                    var y = 0.0;
                    for (i = 0; i < na - 1; i += 3) {
                        arg = (point[0] - params[i + 1]) / params[i + 2];
                        ex = Math.exp(-Math.pow(arg, 2.0));
                        y += params[i] * ex;
                    }

                    return y;
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);

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
            for (var i = 0; i < numParams; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            final var chiSqrDegreesOfFreedom = npoints - numParams;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitMultidimensionalSine() throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
            final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
            final var freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final var freqy = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final var phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
            final var phasey = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[SINE_MULTI_PARAMS];
            params[0] = amplitude;
            params[1] = freqx;
            params[2] = freqy;
            params[3] = phasex;
            params[4] = phasey;

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 2);
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            double error;
            for (var i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

                y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                        Math.sin(freqy * x.getElementAt(i, 1) + phasey);
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                private double[] point;

                private final GradientEstimator gradientEstimator =
                        new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                            @Override
                            public double evaluate(final double[] params) {
                                return evaluateParams(point, params);
                            }
                        });

                @Override
                public int getNumberOfDimensions() {
                    return 2;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[SINE_MULTI_PARAMS];
                    double error;
                    for (var i = 0; i < SINE_MULTI_PARAMS; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives)
                        throws EvaluationException {
                    this.point = point;
                    final var y = evaluateParams(point, params);
                    gradientEstimator.gradient(params, derivatives);

                    return y;
                }

                double evaluateParams(final double[] point, final double[] params) {
                    final var amplitude = params[0];
                    final var freqx = params[1];
                    final var freqy = params[2];
                    final var phasex = params[3];
                    final var phasey = params[4];

                    return amplitude * Math.sin(freqx * point[0] + phasex) * Math.sin(freqy * point[1] + phasey);
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
            var failed = false;
            for (var i = 0; i < SINE_MULTI_PARAMS; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - SINE_MULTI_PARAMS;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitMultidimensionalGaussian() throws FittingException, NotReadyException, WrongSizeException,
            MaxIterationsExceededException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = randomizer.nextInt(MIN_POINTS * 100, MAX_POINTS * 100);
            final var numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final var numParams = numgaussians * GAUSS_MULTI_PARAMS;

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[numParams];
            for (var i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            }

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 2);
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            double error;
            for (var i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
                y[i] = 0.0;
                for (var k = 0; k < numgaussians; k++) {
                    final var b = params[k * GAUSS_MULTI_PARAMS];
                    final var ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    final var ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    final var gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    final var gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    y[i] += b * Math.exp(-(Math.pow(x.getElementAt(i, 0) -
                            ex, 2.0) + Math.pow(x.getElementAt(i, 1) - ey, 2.0)) /
                            Math.pow(gx * gy, 2.0));
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                private double[] point;

                private final GradientEstimator gradientEstimator =
                        new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                            @Override
                            public double evaluate(final double[] params) {
                                return evaluateParams(point, params);
                            }
                        });

                @Override
                public int getNumberOfDimensions() {
                    return NUM_DIMENSIONS;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[numParams];
                    double error;
                    for (var i = 0; i < numParams; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives)
                        throws EvaluationException {
                    this.point = point;
                    final var y = evaluateParams(point, params);
                    gradientEstimator.gradient(params, derivatives);

                    return y;
                }

                double evaluateParams(final double[] point, final double[] params) {
                    var y = 0.0;
                    for (var k = 0; k < numgaussians; k++) {
                        final var b = params[k * GAUSS_MULTI_PARAMS];
                        final var ex = params[k * GAUSS_MULTI_PARAMS + 1];
                        final var ey = params[k * GAUSS_MULTI_PARAMS + 2];
                        final var gx = params[k * GAUSS_MULTI_PARAMS + 3];
                        final var gy = params[k * GAUSS_MULTI_PARAMS + 4];
                        y += b * Math.exp(-(Math.pow(point[0] - ex, 2.0) +
                                Math.pow(point[1] - ey, 2.0)) / Math.pow(gx * gy, 2.0));
                    }

                    return y;
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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
            var failed = false;
            for (var i = 0; i < numParams; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - numParams;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitMultidimensionalSineRepeatInOneDimension() throws FittingException, NotReadyException,
            WrongSizeException, MaxIterationsExceededException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
        final var freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
        final var phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[SINE_UNI_PARAMS];
        params[0] = amplitude;
        params[1] = freqx;
        params[2] = phasex;

        final var y = new double[npoints];
        final var x = new Matrix(npoints, 2);
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        double error;
        for (var i = 0; i < npoints; i++) {
            x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));
            x.setElementAt(i, 1, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

            y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex);
            error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

            private double[] point;

            private final GradientEstimator gradientEstimator =
                    new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                        @Override
                        public double evaluate(final double[] params) {
                            return evaluateParams(point, params);
                        }
                    });

            @Override
            public int getNumberOfDimensions() {
                return 2;
            }

            @Override
            public double[] createInitialParametersArray() {
                final var initParams = new double[SINE_UNI_PARAMS];
                double error;
                for (var i = 0; i < SINE_UNI_PARAMS; i++) {
                    error = errorRandomizer.nextDouble();
                    initParams[i] = params[i] + error;
                }
                return initParams;
            }

            @Override
            public double evaluate(
                    final int i, final double[] point, final double[] params, final double[] derivatives)
                    throws EvaluationException {
                this.point = point;
                final var y = evaluateParams(point, params);
                gradientEstimator.gradient(params, derivatives);

                return y;
            }

            double evaluateParams(final double[] point, final double[] params) {
                final var amplitude = params[0];
                final var freqx = params[1];
                final var phasex = params[2];

                return amplitude * Math.sin(freqx * point[0] + phasex);
            }
        };

        final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, 1.0);
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

        final var chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
        final var chiSqr = fitter.getChisq();

        // probability that chi square can be smaller
        // (the smaller is p the better)
        final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
        assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

        // measure of quality (1.0 indicates maximum quality)
        final var q = 1.0 - p;
        assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

        LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                chiSqr, p * 100.0, q * 100.0));
    }

    @Test
    void testFitUnidimensionalConstantCovariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var constant = randomizer.nextDouble(MIN_CONSTANT, MAX_CONSTANT);

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[CONSTANT_PARAMS];
            params[0] = constant;

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            final var dist = new NormalDist();
            double error;
            for (var i = 0; i < npoints; i++) {
                x.setElementAt(i, 0, randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

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

                // expression below is equal to y[i] = constant
                y[i] = dist.getMean();
                assertEquals(y[i], constant, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[CONSTANT_PARAMS];
                    double error;
                    for (var i = 0; i < CONSTANT_PARAMS; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives) {
                    final var constant = params[0];

                    // derivative of evaluated function respect constant parameter
                    derivatives[0] = 1.0;

                    // evaluated function
                    return constant;
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            assertEquals(CONSTANT_PARAMS, fitter.getA().length);
            for (var i = 0; i < CONSTANT_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], SMALL_ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(CONSTANT_PARAMS, fitter.getCovar().getRows());
            assertEquals(CONSTANT_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0.0);

            final var chiSqrDegreesOfFreedom = npoints - CONSTANT_PARAMS;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var point = new double[evaluator.getNumberOfDimensions()];
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            final var standardDeviation3 = Math.sqrt(cov.getElementAt(0, 0));

            assertEquals(standardDeviation3, sigma, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("real parameter: %f, estimated parameter: %f", params[0],
                    fitter.getA()[0]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma: %f, estimated parameter sigma: %f", sigma,
                    standardDeviation3));

            numValid++;

            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitUnidimensionalLine1Covariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var a = randomizer.nextDouble(MIN_LINE1_A, MAX_LINE1_A);

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[LINE1_PARAMS];
            params[0] = a;

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            final var dist = new NormalDist();
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
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

                // expression below is equal to y[i] = a * x[i]
                y[i] = dist.getMean();
                assertEquals(y[i], a * xi, SMALL_ABSOLUTE_ERROR);

                sigmas[i] = dist.getStandardDeviation();

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[LINE1_PARAMS];
                    double error;
                    for (var i = 0; i < LINE1_PARAMS; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives) {
                    final var a = params[0];

                    // derivative of evaluated function f(x) respect parameter a
                    derivatives[0] = point[0];

                    // evaluated function
                    return a * point[0];
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            assertEquals(LINE1_PARAMS, fitter.getA().length);
            for (var i = 0; i < LINE1_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(LINE1_PARAMS, fitter.getCovar().getRows());
            assertEquals(LINE1_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final var chiSqrDegreesOfFreedom = npoints - LINE1_PARAMS;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var point = new double[evaluator.getNumberOfDimensions()];
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);

                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            final var standardDeviation3 = Math.sqrt(cov.getElementAt(0, 0));

            assertEquals(standardDeviation3, sigma, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("real parameter: %f, estimated parameter: %f",
                    params[0], fitter.getA()[0]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma: %f, estimated parameter sigma: %f",
                    sigma, standardDeviation3));

            numValid++;

            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitUnidimensionalLine2Covariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var a = randomizer.nextDouble(MIN_LINE2_A, MAX_LINE2_A);
            final var b = randomizer.nextDouble(MIN_LINE2_B, MAX_LINE2_B);

            final var sigmaA = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaB = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var varianceA = sigmaA * sigmaA;
            final var varianceB = sigmaB * sigmaB;

            final var params = new double[LINE2_PARAMS];
            params[0] = a;
            params[1] = b;

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, 1.0);
            final var dist = new MultivariateNormalDist();
            final var covariance = Matrix.diagonal(new double[]{varianceA, varianceB});
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(final double[] params, final double[] y, final Matrix jacobian) {
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
                assertEquals(y[i], a * xi + b, SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {
                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[LINE2_PARAMS];
                    for (var i = 0; i < LINE2_PARAMS; i++) {
                        final var error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives) {
                    final var a = params[0];
                    final var b = params[1];

                    // derivatives of function f(x) respect parameters a and b
                    derivatives[0] = point[0];
                    derivatives[1] = 1.0;

                    // evaluated function f(x) = a * x + b
                    return a * point[0] + b;
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            for (var i = 0; i < LINE2_PARAMS; i++) {
                assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
            }
            assertNotNull(fitter.getCovar());
            assertEquals(LINE2_PARAMS, fitter.getCovar().getRows());
            assertEquals(LINE2_PARAMS, fitter.getCovar().getColumns());
            assertTrue(fitter.getChisq() > 0);

            final var chiSqrDegreesOfFreedom = npoints - LINE2_PARAMS;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var point = new double[evaluator.getNumberOfDimensions()];
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);


                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            final var standardDeviationA = Math.sqrt(cov.getElementAt(0, 0));
            final var standardDeviationB = Math.sqrt(cov.getElementAt(1, 1));

            LOGGER.log(Level.INFO, String.format("real parameter A: %f, estimated parameter A: %f",
                    params[0], fitter.getA()[0]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma A: %f, estimated parameter sigma A: %f",
                    sigmaA, standardDeviationA));


            LOGGER.log(Level.INFO, String.format("real parameter B: %f, estimated parameter B: %f",
                    params[1], fitter.getA()[1]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma B: %f, estimated parameter sigma B: %f",
                    sigmaB, standardDeviationB));

            assertEquals(standardDeviationA, sigmaA, SMALL_ABSOLUTE_ERROR);

            numValid++;

            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitUnidimensionalSineCovariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
            final var freq = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final var phase = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            final var sigmaAmplitude = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaFreq = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaPhase = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var varianceAmplitude = sigmaAmplitude * sigmaAmplitude;
            final var varianceFreq = sigmaFreq * sigmaFreq;
            final var variancePhase = sigmaPhase * sigmaPhase;

            final var params = new double[SINE_UNI_PARAMS];
            params[0] = amplitude;
            params[1] = freq;
            params[2] = phase;

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, 1.0);
            final var dist = new MultivariateNormalDist();
            final var covariance = Matrix.diagonal(new double[]{varianceAmplitude, varianceFreq, variancePhase});
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(final double[] params, final double[] y, final Matrix jacobian) {
                        // y[i] = amplitude * Math.sin(freq * xi + phase)
                        y[0] = params[0] * Math.sin(params[1] * xi + params[2]);

                        // derivatives respect parameters
                        jacobian.setElementAt(0, 0, Math.sin(params[1] * xi + params[2]));
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
                assertEquals(y[i], amplitude * Math.sin(freq * xi + phase), SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[SINE_UNI_PARAMS];
                    double error;
                    for (var i = 0; i < SINE_UNI_PARAMS; i++) {
                        error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives) {
                    final var amplitude = params[0];
                    final var freq = params[1];
                    final var phase = params[2];
                    final var y = amplitude * Math.sin(freq * point[0] + phase);

                    // derivative respect amplitude
                    derivatives[0] = Math.sin(freq * point[0] + phase);
                    // derivative respect frequency
                    derivatives[1] = amplitude * Math.cos(freq * point[0] + phase) * point[0];
                    // derivative respect phase
                    derivatives[2] = amplitude * Math.cos(freq * point[0] + phase);

                    return y;
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            var valid = true;
            for (var i = 0; i < SINE_UNI_PARAMS; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - SINE_UNI_PARAMS;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var point = new double[evaluator.getNumberOfDimensions()];
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);

                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            final var standardDeviationAmplitude = Math.sqrt(cov.getElementAt(0, 0));
            final var standardDeviationFreq = Math.sqrt(cov.getElementAt(1, 1));
            final var standardDeviationPhase = Math.sqrt(cov.getElementAt(2, 2));

            LOGGER.log(Level.INFO, String.format("real parameter Amplitude: %f, estimated parameter Amplitude: %f",
                    params[0], fitter.getA()[0]));
            LOGGER.log(Level.INFO, String.format(
                    "real parameter sigma Amplitude: %f, estimated parameter sigma Amplitude: %f",
                    sigmaAmplitude, standardDeviationAmplitude));

            LOGGER.log(Level.INFO, String.format("real parameter Freq: %f, estimated parameter Freq: %f",
                    params[1], fitter.getA()[1]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma Freq: %f, estimated parameter sigma Freq: %f",
                    sigmaFreq, standardDeviationFreq));

            LOGGER.log(Level.INFO, String.format("real parameter Phase: %f, estimated parameter Phase: %f",
                    params[2], fitter.getA()[2]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma Phase: %f, estimated parameter sigma Phase: %f",
                    sigmaPhase, standardDeviationPhase));

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitUnidimensionalGaussianCovariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final var numParams = numgaussians * GAUSS_UNI_PARAMS;

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[numParams];
            final var varianceParams = new double[numParams];
            for (var i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 1);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            final var dist = new MultivariateNormalDist();
            final var covariance = Matrix.diagonal(varianceParams);
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi);

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {
                    @Override
                    public void evaluate(final double[] params, final double[] y, final Matrix jacobian) {
                        y[0] = 0.0;
                        for (var k = 0; k < numgaussians; k++) {
                            final var b = params[k * GAUSS_UNI_PARAMS];
                            final var e = params[k * GAUSS_UNI_PARAMS + 1];
                            final var g = params[k * GAUSS_UNI_PARAMS + 2];
                            y[0] += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                        }

                        int i;
                        final var na = params.length;
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

                var yi = 0.0;
                for (var k = 0; k < numgaussians; k++) {
                    final var b = params[k * GAUSS_UNI_PARAMS];
                    final var e = params[k * GAUSS_UNI_PARAMS + 1];
                    final var g = params[k * GAUSS_UNI_PARAMS + 2];
                    yi += b * Math.exp(-Math.pow((xi - e) / g, 2.0));
                }

                y[i] = dist.getMean()[0];
                assertEquals(y[i], yi, SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                @Override
                public int getNumberOfDimensions() {
                    return 1;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[numParams];
                    for (var i = 0; i < numParams; i++) {
                        final var error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int pos, final double[] point, final double[] params, final double[] derivatives) {
                    int i;
                    final var na = params.length;
                    double fac;
                    double ex;
                    double arg;
                    var y = 0.0;
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

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            var valid = true;
            for (var i = 0; i < numParams; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - numParams;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var point = new double[evaluator.getNumberOfDimensions()];
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);

                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));


            for (var i = 0; i < numParams; i++) {
                final var standardDeviation = Math.sqrt(cov.getElementAt(i, i));

                LOGGER.log(Level.INFO, String.format("real parameter: %f, estimated parameter: %f",
                        params[i], fitter.getA()[i]));
                LOGGER.log(Level.INFO, String.format("real parameter sigma: %f, estimated parameter sigma: %f",
                        sigma, standardDeviation));
            }

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitMultidimensionalSineCovariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var amplitude = randomizer.nextDouble(MIN_SINE_AMPLITUDE, MAX_SINE_AMPLITUDE);
            final var freqx = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final var freqy = randomizer.nextDouble(MIN_SINE_FREQ, MAX_SINE_FREQ);
            final var phasex = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);
            final var phasey = randomizer.nextDouble(MIN_SINE_PHASE, MAX_SINE_PHASE);

            final var sigmaAmplitude = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaFreqx = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaFreqy = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaPhasex = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);
            final var sigmaPhasey = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var varianceAmplitude = sigmaAmplitude * sigmaAmplitude;
            final var varianceFreqx = sigmaFreqx * sigmaFreqx;
            final var varianceFreqy = sigmaFreqy * sigmaFreqy;
            final var variancePhasex = sigmaPhasex * sigmaPhasex;
            final var variancePhasey = sigmaPhasey * sigmaPhasey;

            final var params = new double[SINE_MULTI_PARAMS];
            params[0] = amplitude;
            params[1] = freqx;
            params[2] = freqy;
            params[3] = phasex;
            params[4] = phasey;

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 2);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, 1.0);
            final var point = new double[2];
            final var dist = new MultivariateNormalDist();
            final var covariance = Matrix.diagonal(new double[]{
                    varianceAmplitude, varianceFreqx, varianceFreqy, variancePhasex, variancePhasey});
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi0 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var xi1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi0);
                x.setElementAt(i, 1, xi1);
                point[0] = xi0;
                point[1] = xi1;

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {

                    final double[] derivatives = new double[params.length];

                    private final GradientEstimator gradientEstimator = new GradientEstimator(
                            params1 -> evaluateParams(point, params1));

                    @Override
                    public void evaluate(final double[] params, final double[] y, final Matrix jacobian) {
                        try {
                            gradientEstimator.gradient(params, derivatives);
                            y[0] = evaluateParams(point, params);
                            jacobian.fromArray(derivatives);
                        } catch (final Exception ignore) {
                            // never happens
                        }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }

                    double evaluateParams(final double[] point, final double[] params) {
                        final var amplitude = params[0];
                        final var freqx = params[1];
                        final var freqy = params[2];
                        final var phasex = params[3];
                        final var phasey = params[4];

                        return amplitude * Math.sin(freqx * point[0] + phasex) * Math.sin(freqy * point[1] + phasey);
                    }

                }, params, covariance, dist);

                // expression below is equal to:
                // y[i] = amplitude * Math.sin(freqx * x.getElementAt(i, 0) + phasex) *
                //      Math.sin(freqy * x.getElementAt(i, 1) + phasey)
                y[i] = dist.getMean()[0];
                assertEquals(y[i], amplitude * Math.sin(freqx * xi0 + phasex) * Math.sin(freqy * xi1 + phasey),
                        SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                private double[] point;

                private final GradientEstimator gradientEstimator =
                        new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                            @Override
                            public double evaluate(final double[] params) {
                                return evaluateParams(point, params);
                            }
                        });

                @Override
                public int getNumberOfDimensions() {
                    return 2;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[SINE_MULTI_PARAMS];
                    for (var i = 0; i < SINE_MULTI_PARAMS; i++) {
                        final var error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives)
                        throws EvaluationException {
                    this.point = point;
                    final var y = evaluateParams(point, params);
                    gradientEstimator.gradient(params, derivatives);

                    return y;
                }

                double evaluateParams(final double[] point, final double[] params) {
                    final var amplitude = params[0];
                    final var freqx = params[1];
                    final var freqy = params[2];
                    final var phasex = params[3];
                    final var phasey = params[4];

                    return amplitude * Math.sin(freqx * point[0] + phasex) * Math.sin(freqy * point[1] + phasey);
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            assertEquals(SINE_MULTI_PARAMS, fitter.getA().length);
            var valid = true;
            for (var i = 0; i < SINE_MULTI_PARAMS; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - SINE_MULTI_PARAMS;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);

                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            final var standardDeviationAmplitude = Math.sqrt(cov.getElementAt(0, 0));
            final var standardDeviationFreqx = Math.sqrt(cov.getElementAt(1, 1));
            final var standardDeviationFreqy = Math.sqrt(cov.getElementAt(2, 2));
            final var standardDeviationPhasex = Math.sqrt(cov.getElementAt(3, 3));
            final var standardDeviationPhasey = Math.sqrt(cov.getElementAt(4, 4));

            LOGGER.log(Level.INFO, String.format("real parameter Amplitude: %f, estimated parameter Amplitude: %f",
                    params[0], fitter.getA()[0]));
            LOGGER.log(Level.INFO, String.format(
                    "real parameter sigma Amplitude: %f, estimated parameter sigma Amplitude: %f",
                    sigmaAmplitude, standardDeviationAmplitude));

            LOGGER.log(Level.INFO, String.format("real parameter Freqx: %f, estimated parameter Freqx: %f",
                    params[1], fitter.getA()[1]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma Freqx: %f, estimated parameter sigma Freqx: %f",
                    sigmaFreqx, standardDeviationFreqx));

            LOGGER.log(Level.INFO, String.format("real parameter Freqy: %f, estimated parameter Freqy: %f",
                    params[2], fitter.getA()[2]));
            LOGGER.log(Level.INFO, String.format("real parameter sigma Freqy: %f, estimated parameter sigma Freqy: %f",
                    sigmaFreqy, standardDeviationFreqy));

            LOGGER.log(Level.INFO, String.format("real parameter Phasex: %f, estimated parameter Phasex: %f",
                    params[3], fitter.getA()[3]));
            LOGGER.log(Level.INFO, String.format(
                    "real parameter sigma Phasex: %f, estimated parameter sigma Phasex: %f",
                    sigmaPhasex, standardDeviationPhasex));

            LOGGER.log(Level.INFO, String.format("real parameter Phasey: %f, estimated parameter Phasey: %f",
                    params[4], fitter.getA()[4]));
            LOGGER.log(Level.INFO, String.format(
                    "real parameter sigma Phasey: %f, estimated parameter sigma Phasey: %f",
                    sigmaPhasey, standardDeviationPhasey));

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testFitMultidimensionalGaussianCovariance() throws Throwable {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = N_SAMPLES;
            final var numgaussians = randomizer.nextInt(MIN_GAUSSIANS, MAX_GAUSSIANS);
            final var numParams = numgaussians * GAUSS_MULTI_PARAMS;

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[numParams];
            final var varianceParams = new double[numParams];
            for (var i = 0; i < numParams; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            }
            Arrays.fill(varianceParams, sigma * sigma);

            final var y = new double[npoints];
            final var x = new Matrix(npoints, 2);
            final var sigmas = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, 1.0);
            final var point = new double[2];
            final var dist = new MultivariateNormalDist();
            final var covariance = Matrix.diagonal(varianceParams);
            double error;
            for (var i = 0; i < npoints; i++) {
                final var xi0 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                final var xi1 = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                x.setElementAt(i, 0, xi0);
                x.setElementAt(i, 1, xi1);
                point[0] = xi0;
                point[1] = xi1;

                // propagate standard deviation of a (sigma) for x[i] value:
                MultivariateNormalDist.propagate(new MultivariateNormalDist.JacobianEvaluator() {

                    final double[] derivatives = new double[params.length];

                    private final GradientEstimator gradientEstimator = new GradientEstimator(
                            params1 -> evaluateParams(point, params1));

                    @Override
                    public void evaluate(final double[] params, final double[] y, final Matrix jacobian) {
                        try {
                            gradientEstimator.gradient(params, derivatives);
                            y[0] = evaluateParams(point, params);
                            jacobian.fromArray(derivatives);
                        } catch (final Exception ignore) {
                            // never happens
                        }
                    }

                    @Override
                    public int getNumberOfVariables() {
                        return 1;
                    }

                    double evaluateParams(final double[] point, final double[] params) {
                        var y = 0.0;
                        for (var k = 0; k < numgaussians; k++) {
                            final var b = params[k * GAUSS_MULTI_PARAMS];
                            final var ex = params[k * GAUSS_MULTI_PARAMS + 1];
                            final var ey = params[k * GAUSS_MULTI_PARAMS + 2];
                            final var gx = params[k * GAUSS_MULTI_PARAMS + 3];
                            final var gy = params[k * GAUSS_MULTI_PARAMS + 4];
                            y += b * Math.exp(-(Math.pow(point[0] - ex, 2.0) + Math.pow(point[1] - ey, 2.0))
                                    / Math.pow(gx * gy, 2.0));
                        }

                        return y;
                    }

                }, params, covariance, dist);

                var yi = 0.0;
                for (var k = 0; k < numgaussians; k++) {
                    final var b = params[k * GAUSS_MULTI_PARAMS];
                    final var ex = params[k * GAUSS_MULTI_PARAMS + 1];
                    final var ey = params[k * GAUSS_MULTI_PARAMS + 2];
                    final var gx = params[k * GAUSS_MULTI_PARAMS + 3];
                    final var gy = params[k * GAUSS_MULTI_PARAMS + 4];
                    yi += b * Math.exp(-(Math.pow(x.getElementAt(i, 0)
                            - ex, 2.0) + Math.pow(x.getElementAt(i, 1) - ey, 2.0))
                            / Math.pow(gx * gy, 2.0));
                }

                y[i] = dist.getMean()[0];
                assertEquals(y[i], yi, SMALL_ABSOLUTE_ERROR);

                assertEquals(1, dist.getCovariance().getRows());
                assertEquals(1, dist.getCovariance().getColumns());

                sigmas[i] = Math.sqrt(dist.getCovariance().getElementAt(0, 0));

                errorRandomizer.setStandardDeviation(Math.max(sigmas[i], Double.MIN_VALUE));

                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LevenbergMarquardtMultiDimensionFunctionEvaluator() {

                private double[] point;

                private final GradientEstimator gradientEstimator =
                        new GradientEstimator(new MultiDimensionFunctionEvaluatorListener() {

                            @Override
                            public double evaluate(final double[] params) {
                                return evaluateParams(point, params);
                            }
                        });

                @Override
                public int getNumberOfDimensions() {
                    return NUM_DIMENSIONS;
                }

                @Override
                public double[] createInitialParametersArray() {
                    final var initParams = new double[numParams];
                    for (var i = 0; i < numParams; i++) {
                        final var error = errorRandomizer.nextDouble();
                        initParams[i] = params[i] + error;
                    }
                    return initParams;
                }

                @Override
                public double evaluate(
                        final int i, final double[] point, final double[] params, final double[] derivatives)
                        throws EvaluationException {
                    this.point = point;
                    final var y = evaluateParams(point, params);
                    gradientEstimator.gradient(params, derivatives);

                    return y;
                }

                double evaluateParams(final double[] point, final double[] params) {
                    var y = 0.0;
                    for (var k = 0; k < numgaussians; k++) {
                        final var b = params[k * GAUSS_MULTI_PARAMS];
                        final var ex = params[k * GAUSS_MULTI_PARAMS + 1];
                        final var ey = params[k * GAUSS_MULTI_PARAMS + 2];
                        final var gx = params[k * GAUSS_MULTI_PARAMS + 3];
                        final var gy = params[k * GAUSS_MULTI_PARAMS + 4];
                        y += b * Math.exp(-(Math.pow(point[0] - ex, 2.0) + Math.pow(point[1] - ey, 2.0))
                                / Math.pow(gx * gy, 2.0));
                    }

                    return y;
                }
            };

            final var fitter = new LevenbergMarquardtMultiDimensionFitter(evaluator, x, y, sigmas);
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
            var valid = true;
            for (var i = 0; i < numParams; i++) {
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

            final var chiSqrDegreesOfFreedom = npoints - numParams;
            final var chiSqr = fitter.getChisq();

            // probability that chi square can be smaller
            // (the smaller is p the better)
            final var p = ChiSqDist.cdf(chiSqr, chiSqrDegreesOfFreedom);
            assertEquals(fitter.getP(), p, SMALL_ABSOLUTE_ERROR);

            // measure of quality (1.0 indicates maximum quality)
            final var q = 1.0 - p;
            assertEquals(fitter.getQ(), q, SMALL_ABSOLUTE_ERROR);

            LOGGER.log(Level.INFO, String.format("chi sqr: %f, probability smaller chi sqr: %f %%, quality: %f %%",
                    chiSqr, p * 100.0, q * 100.0));

            final var mse = fitter.getMse();

            // Covariance of parameters is
            // Vp = (J’ * W * J)^-1

            // where J is the jacobian of measures, each row contains derivatives
            // for one sample, each column is the derivative for one parameter
            // W is the inverse of the measurement error covariance. Assuming
            // independent samples, W is diagonal with diagonal terms 1/sigma^2
            // where sigma is the standard deviation of each sample
            // More info: http://people.duke.edu/~hpgavin/ce281/lm.pdf

            final var invCov = new Matrix(params.length, params.length);
            final var tmp1 = new Matrix(params.length, 1);
            final var tmp2 = new Matrix(1, params.length);
            final var tmpInvCov = new Matrix(params.length, params.length);
            final var estimatedParams = fitter.getA();
            final var derivatives = new double[params.length];
            var mse2 = 0.0;
            for (var i = 0; i < npoints; i++) {
                x.getSubmatrixAsArray(i, 0, i, evaluator.getNumberOfDimensions() - 1,
                        point);
                final var yi = evaluator.evaluate(i, point, estimatedParams, derivatives);

                tmp1.fromArray(derivatives);
                tmp2.fromArray(derivatives);

                tmp1.multiply(tmp2, tmpInvCov);

                final var w = 1.0 / ((chiSqrDegreesOfFreedom + 1) * sigmas[i] * sigmas[i]);
                tmpInvCov.multiplyByScalar(w);
                invCov.add(tmpInvCov);

                mse2 += Math.pow(y[i] - yi, 2.0) / (npoints - params.length);
            }

            assertEquals(mse, mse2, SMALL_ABSOLUTE_ERROR);

            final var cov = Utils.inverse(invCov);
            assertTrue(cov.equals(fitter.getCovar(), SMALL_ABSOLUTE_ERROR));

            for (var i = 0; i < numParams; i++) {
                final var standardDeviation = Math.sqrt(cov.getElementAt(i, i));

                LOGGER.log(Level.INFO, String.format("real parameter: %f, estimated parameter: %f",
                        params[i], fitter.getA()[i]));
                LOGGER.log(Level.INFO, String.format("real parameter sigma: %f, estimated parameter sigma: %f",
                        sigma, standardDeviation));
            }

            if (valid) {
                numValid++;
                break;
            }
        }

        assertTrue(numValid > 0);
    }
}
