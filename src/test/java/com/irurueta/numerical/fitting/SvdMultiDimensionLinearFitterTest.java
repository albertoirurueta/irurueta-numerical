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
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.GaussianRandomizer;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.*;

public class SvdMultiDimensionLinearFitterTest {

    private static final int MIN_POINTS = 500;
    private static final int MAX_POINTS = 1000;

    private static final double MIN_RANDOM_VALUE = -100.0;
    private static final double MAX_RANDOM_VALUE = 100.0;

    private static final double MIN_SIGMA_VALUE = 1e-4;
    private static final double MAX_SIGMA_VALUE = 1.0;

    private static final double ABSOLUTE_ERROR = 1e-1;

    // For functions like: a*1 + b*x + c*y + d*x*y + e*x^2 + f*y^2
    private static final int NUM_QUADRATIC_PARAMS = 6;

    private static final int TIMES = 10;

    @Test
    public void testConstructor() throws FittingException, WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        // test empty constructor
        SvdMultiDimensionLinearFitter fitter =
                new SvdMultiDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // test constructor with input data
        final Matrix x = new Matrix(nPoints, 2);
        final double[] y = new double[nPoints];
        final double[] sig = new double[nPoints];

        fitter = new SvdMultiDimensionLinearFitter(x, y, sig);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertFalse(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        final Matrix shortX = new Matrix(nPoints - 1, 2);
        final double[] shortY = new double[nPoints - 1];
        final double[] shortSig = new double[nPoints - 1];

        fitter = null;
        try {
            fitter = new SvdMultiDimensionLinearFitter(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SvdMultiDimensionLinearFitter(x, shortY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SvdMultiDimensionLinearFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with input data (constant sigma)
        fitter = new SvdMultiDimensionLinearFitter(x, y, 1.0);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertFalse(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new SvdMultiDimensionLinearFitter(shortX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SvdMultiDimensionLinearFitter(x, shortY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator
        final LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

                    @Override
                    public int getNumberOfDimensions() {
                        return 2;
                    }

                    @Override
                    public double[] createResultArray() {
                        return new double[nPoints];
                    }

                    @Override
                    public void evaluate(final double[] point, final double[] result) {
                    }
                };

        fitter = new SvdMultiDimensionLinearFitter(evaluator);

        // check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // test constructor with evaluator and input data
        fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, sig);

        // check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new SvdMultiDimensionLinearFitter(evaluator, shortX, y,
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SvdMultiDimensionLinearFitter(evaluator, x, shortY,
                    sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y,
                    shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with evaluator and input data (constant sigma)
        fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, 1.0);

        // check default values
        assertSame(fitter.getFunctionEvaluator(), evaluator);
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNotNull(fitter.getSig());
        for (int i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertTrue(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        fitter = null;
        try {
            fitter = new SvdMultiDimensionLinearFitter(evaluator, shortX, y,
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter = new SvdMultiDimensionLinearFitter(evaluator, x, shortY,
                    1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);
    }

    @Test
    public void testGetSetFunctionEvaluator() throws FittingException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final SvdMultiDimensionLinearFitter fitter =
                new SvdMultiDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());

        // set new value
        final LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

                    @Override
                    public int getNumberOfDimensions() {
                        return 2;
                    }

                    @Override
                    public double[] createResultArray() {
                        return new double[nPoints];
                    }

                    @Override
                    public void evaluate(final double[] point, final double[] result) {
                    }
                };

        // set new value
        fitter.setFunctionEvaluator(evaluator);

        // check correctness
        assertSame(fitter.getFunctionEvaluator(), evaluator);
    }

    @Test
    public void testGetSetInputData() throws WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final SvdMultiDimensionLinearFitter fitter =
                new SvdMultiDimensionLinearFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final Matrix x = new Matrix(nPoints, 2);
        final double[] y = new double[nPoints];
        final double[] sig = new double[nPoints];

        // set input data
        fitter.setInputData(x, y, sig);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);

        // Force IllegalArgumentException
        final Matrix wrongX = new Matrix(nPoints - 1, 2);
        final double[] wrongY = new double[nPoints - 1];
        final double[] wrongSig = new double[nPoints - 1];

        try {
            fitter.setInputData(wrongX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, wrongY, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, y, wrongSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetInputDataWithConstantSigma()
            throws WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final SvdMultiDimensionLinearFitter fitter =
                new SvdMultiDimensionLinearFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final Matrix x = new Matrix(nPoints, 2);
        final double[] y = new double[nPoints];

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
        final Matrix wrongX = new Matrix(nPoints - 1, 2);
        final double[] wrongY = new double[nPoints - 1];

        try {
            fitter.setInputData(wrongX, y, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputData(x, wrongY, 1.0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testIsReady() throws FittingException, WrongSizeException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final SvdMultiDimensionLinearFitter fitter =
                new SvdMultiDimensionLinearFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final Matrix x = new Matrix(nPoints, 2);
        final double[] y = new double[nPoints];
        final double[] sig = new double[nPoints];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final LinearFitterMultiDimensionFunctionEvaluator evaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

                    @Override
                    public int getNumberOfDimensions() {
                        return 2;
                    }

                    @Override
                    public double[] createResultArray() {
                        return new double[nPoints];
                    }

                    @Override
                    public void evaluate(final double[] point, final double[] result) {
                    }
                };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());

        // test bad evaluator
        final LinearFitterMultiDimensionFunctionEvaluator badEvaluator =
                new LinearFitterMultiDimensionFunctionEvaluator() {

                    @Override
                    public int getNumberOfDimensions() {
                        return 3;
                    }

                    @Override
                    public double[] createResultArray() {
                        return new double[nPoints];
                    }

                    @Override
                    public void evaluate(final double[] point, final double[] result) {
                    }
                };

        fitter.setFunctionEvaluator(badEvaluator);

        assertFalse(fitter.isReady());
    }

    @Test
    public void testGetSetTol() {
        final SvdMultiDimensionLinearFitter fitter =
                new SvdMultiDimensionLinearFitter();

        // check default values
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // set new value
        fitter.setTol(1e-3);

        // check correctness
        assertEquals(1e-3, fitter.getTol(), 0.0);
    }

    @Test
    public void testFitQuadratic() throws FittingException, NotReadyException,
            WrongSizeException {

        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {
            final UniformRandomizer randomizer = new UniformRandomizer(new Random());

            final int npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

            final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final double[] params = new double[NUM_QUADRATIC_PARAMS];
            for (int i = 0; i < NUM_QUADRATIC_PARAMS; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE,
                        MAX_RANDOM_VALUE);
            }

            final Matrix x = Matrix.createWithUniformRandomValues(npoints, 2,
                    MIN_RANDOM_VALUE, MAX_RANDOM_VALUE, new Random());
            final double[] y = new double[npoints];
            final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                    new Random(), 0.0, sigma);
            double error, x1, x2;
            for (int i = 0; i < npoints; i++) {
                x1 = x.getElementAt(i, 0);
                x2 = x.getElementAt(i, 1);
                // function is: a*1 + b*x + c*y + d*x*y + e*x^2 + f*y^2
                y[i] = params[0] + params[1] * x1 + params[2] * x2 + params[3] * x1 * x2 +
                        params[4] * x1 * x1 + params[5] * x2 * x2;
                error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final LinearFitterMultiDimensionFunctionEvaluator evaluator =
                    new LinearFitterMultiDimensionFunctionEvaluator() {

                        @Override
                        public int getNumberOfDimensions() {
                            return 2; // x, y coordinates
                        }

                        @Override
                        public double[] createResultArray() {
                            // parameters a, b, c, d, e, f for function:
                            // a + b*x + c*y + d*x*y + e*x*x + f*y*y
                            return new double[NUM_QUADRATIC_PARAMS];
                        }

                        @Override
                        public void evaluate(final double[] point, final double[] result) {
                            final double x = point[0];
                            final double y = point[1];

                            result[0] = 1.0;
                            result[1] = x;
                            result[2] = y;
                            result[3] = x * y;
                            result[4] = x * x;
                            result[5] = y * y;
                        }
                    };

            final SvdMultiDimensionLinearFitter fitter =
                    new SvdMultiDimensionLinearFitter(evaluator, x, y, 1.0);

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
            assertEquals(NUM_QUADRATIC_PARAMS, fitter.getA().length);

            boolean failed = false;
            for (int i = 0; i < NUM_QUADRATIC_PARAMS; i++) {
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
}
