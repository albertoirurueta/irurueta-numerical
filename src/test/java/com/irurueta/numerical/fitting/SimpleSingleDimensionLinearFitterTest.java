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

import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class SimpleSingleDimensionLinearFitterTest {
    private static final int MIN_POLY_PARAMS = 1;
    private static final int MAX_POLY_PARAMS = 5;

    private static final int MIN_POINTS = 500;
    private static final int MAX_POINTS = 1000;

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;

    private static final double MIN_SIGMA_VALUE = 1e-4;
    private static final double MAX_SIGMA_VALUE = 1.0;

    private static final double ABSOLUTE_ERROR = 1e-1;

    private static final int TRIGO_PARAMS = 2;

    private static final int TIMES = 10;

    @Test
    public void testConstructor() throws FittingException {
        // test empty constructor
        SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());

        // test constructor with input data
        final double[] x = new double[2];
        final double[] y = new double[2];
        final double[] sig = new double[2];

        fitter = new SimpleSingleDimensionLinearFitter(x, y, sig);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertFalse(fitter.isReady());

        // Force IllegalArgumentException
        final double[] shortX = new double[1];
        final double[] shortY = new double[1];
        final double[] shortSig = new double[1];

        fitter = null;
        try {
            fitter = new SimpleSingleDimensionLinearFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SimpleSingleDimensionLinearFitter(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SimpleSingleDimensionLinearFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with input data (constant sigma)
        fitter = new SimpleSingleDimensionLinearFitter(x, y, 1.0);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertFalse(fitter.isReady());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new SimpleSingleDimensionLinearFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SimpleSingleDimensionLinearFitter(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator
        final LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createResultArray() {
                        return new double[2];
                    }

                    @Override
                    public void evaluate(final double point, final double[] result) {
                    }
                };

        fitter = new SimpleSingleDimensionLinearFitter(evaluator);

        // check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());

        // test constructor with evaluator and input data
        fitter = new SimpleSingleDimensionLinearFitter(evaluator, x, y, sig);

        // check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new SimpleSingleDimensionLinearFitter(evaluator, shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SimpleSingleDimensionLinearFitter(evaluator, x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SimpleSingleDimensionLinearFitter(evaluator, x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator and input data (constant sigma)
        fitter = new SimpleSingleDimensionLinearFitter(evaluator, x, y, 1.0);

        // check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(fitter.getSig()[i], 1.0, 0.0);
        }
        assertTrue(fitter.isReady());

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new SimpleSingleDimensionLinearFitter(evaluator, shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SimpleSingleDimensionLinearFitter(evaluator, x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);
    }

    @Test
    public void testGetSetFunctionEvaluator() throws FittingException {
        final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());

        // set new value
        final LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createResultArray() {
                        return new double[2];
                    }

                    @Override
                    public void evaluate(final double point, final double[] result) {
                    }
                };

        // set new value
        fitter.setFunctionEvaluator(evaluator);

        // check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
    }

    @Test
    public void testGetSetInputData() {
        final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter();

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
        final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter();

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
        final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final double[] x = new double[2];
        final double[] y = new double[2];
        final double[] sig = new double[2];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final LinearFitterSingleDimensionFunctionEvaluator evaluator =
                new LinearFitterSingleDimensionFunctionEvaluator() {

                    @Override
                    public double[] createResultArray() {
                        return new double[2];
                    }

                    @Override
                    public void evaluate(final double point, final double[] result) {
                    }
                };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());
    }

    @Test
    public void testFitPolynomial() throws FittingException, NotReadyException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int nparams = randomizer.nextInt(MIN_POLY_PARAMS,
                    MAX_POLY_PARAMS);
            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[nparams];
            for (int i = 0; i < nparams; i++) {
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
                for (int j = 0; j < nparams; j++) {
                    y[i] += params[j] * Math.pow(x[i], j);
                }
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LinearFitterSingleDimensionFunctionEvaluator evaluator =
                    new LinearFitterSingleDimensionFunctionEvaluator() {
                        @Override
                        public double[] createResultArray() {
                            return new double[nparams];
                        }

                        @Override
                        public void evaluate(final double point, final double[] result) {
                            for (int i = 0; i < result.length; i++) {
                                result[i] = Math.pow(point, i);
                            }
                        }
                    };

            final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter(evaluator, x, y,
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
            assertEquals(fitter.getA().length, nparams);
            boolean failed = false;
            for (int i = 0; i < nparams; i++) {
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

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testTrigo() throws FittingException, NotReadyException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[TRIGO_PARAMS];
            for (int i = 0; i < TRIGO_PARAMS; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            // function: y(x) = a * cos(x) + b * sin(x)
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                error = errorRandomizer.nextDouble();
                y[i] = params[0] * Math.sin(x[i]) + params[1] * Math.cos(x[i]) +
                        error;
            }

            final LinearFitterSingleDimensionFunctionEvaluator evaluator =
                    new LinearFitterSingleDimensionFunctionEvaluator() {

                        @Override
                        public double[] createResultArray() {
                            return new double[TRIGO_PARAMS];
                        }

                        @Override
                        public void evaluate(final double point, final double[] result) {
                            result[0] = Math.sin(point);
                            result[1] = Math.cos(point);
                        }
                    };

            final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter(
                    evaluator, x, y, 1.0);

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
            assertEquals(fitter.getA().length, TRIGO_PARAMS);
            boolean failed = false;
            for (int i = 0; i < TRIGO_PARAMS; i++) {
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

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testFitTrigoWithHoldAndFree() throws FittingException,
            NotReadyException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[TRIGO_PARAMS];
            for (int i = 0; i < TRIGO_PARAMS; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final double[] y = new double[npoints];
            final double[] x = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error;
            // function: y(x) = a * cos(x) + b * sin(x)
            for (int i = 0; i < npoints; i++) {
                x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
                error = errorRandomizer.nextDouble();
                y[i] = params[0] * Math.sin(x[i]) + params[1] * Math.cos(x[i]) +
                        error;
            }

            final LinearFitterSingleDimensionFunctionEvaluator evaluator =
                    new LinearFitterSingleDimensionFunctionEvaluator() {

                        @Override
                        public double[] createResultArray() {
                            return new double[TRIGO_PARAMS];
                        }

                        @Override
                        public void evaluate(final double point, final double[] result) {
                            result[0] = Math.sin(point);
                            result[1] = Math.cos(point);
                        }
                    };

            final SimpleSingleDimensionLinearFitter fitter = new SimpleSingleDimensionLinearFitter(evaluator,
                    x, y, 1.0);

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
            // first parameter is hold and matches exactly
            assertEquals(fitter.getA()[0], params[0], 0.0);
            if (Math.abs(fitter.getA()[1] - params[1]) > ABSOLUTE_ERROR) {
                continue;
            }
            assertEquals(fitter.getA()[1], params[1], ABSOLUTE_ERROR);
            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            // release first parameter
            fitter.free(0);

            // hold 2nd parameter
            fitter.hold(1, params[1]);

            // fit
            fitter.fit();

            // check correctness
            assertTrue(fitter.isResultAvailable());
            assertEquals(fitter.getA()[0], params[0], ABSOLUTE_ERROR);
            // second parameter is hold and matches exactly
            assertEquals(fitter.getA()[1], params[1], 0.0);
            assertNotNull(fitter.getCovar());
            assertTrue(fitter.getChisq() > 0);

            // Force FittingException (by holding all parameters)
            fitter.hold(0, params[0]);
            fitter.hold(1, params[1]);

            try {
                fitter.fit();
                fail("FittingException expected but not thrown");
            } catch (final FittingException ignore) {
            }

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }
}
