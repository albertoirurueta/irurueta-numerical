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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class SvdMultiDimensionLinearFitterTest {

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
    void testConstructor() throws FittingException, WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        // test empty constructor
        var fitter = new SvdMultiDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // test constructor with input data
        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];
        final var sig = new double[nPoints];

        fitter = new SvdMultiDimensionLinearFitter(x, y, sig);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        final var shortX = new Matrix(nPoints - 1, 2);
        final var shortY = new double[nPoints - 1];
        final var shortSig = new double[nPoints - 1];

        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(shortX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(x, shortY, sig));
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(x, y, shortSig));

        // test constructor with input data (constant sigma)
        fitter = new SvdMultiDimensionLinearFitter(x, y, 1.0);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertFalse(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(shortX, y, 1.0));
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(x, shortY, 1.0));

        // test constructor with evaluator
        final var evaluator = new LinearFitterMultiDimensionFunctionEvaluator() {

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
                // no action needed
            }
        };

        fitter = new SvdMultiDimensionLinearFitter(evaluator);

        // check default values
        assertSame(evaluator, fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // test constructor with evaluator and input data
        fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, sig);

        // check default values
        assertSame(evaluator, fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(evaluator, shortX, y,
                sig));
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(evaluator, x, shortY,
                sig));
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(evaluator, x, y,
                shortSig));

        // test constructor with evaluator and input data (constant sigma)
        fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, 1.0);

        // check default values
        assertSame(evaluator, fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertTrue(fitter.isReady());
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(evaluator, shortX, y,
                1.0));
        assertThrows(IllegalArgumentException.class, () -> new SvdMultiDimensionLinearFitter(evaluator, x, shortY,
                1.0));
    }

    @Test
    void testGetSetFunctionEvaluator() throws FittingException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new SvdMultiDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());

        // set new value
        final var evaluator = new LinearFitterMultiDimensionFunctionEvaluator() {

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
                // no action needed
            }
        };

        // set new value
        fitter.setFunctionEvaluator(evaluator);

        // check correctness
        assertSame(evaluator, fitter.getFunctionEvaluator());
    }

    @Test
    void testGetSetInputData() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new SvdMultiDimensionLinearFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];
        final var sig = new double[nPoints];

        // set input data
        fitter.setInputData(x, y, sig);

        // check correctness
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());

        // Force IllegalArgumentException
        final var wrongX = new Matrix(nPoints - 1, 2);
        final var wrongY = new double[nPoints - 1];
        final var wrongSig = new double[nPoints - 1];

        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(wrongX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, wrongY, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, y, wrongSig));
    }

    @Test
    void testGetSetInputDataWithConstantSigma() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new SvdMultiDimensionLinearFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];

        // set input data
        fitter.setInputData(x, y, 1.0);

        // check correctness
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }

        // Force IllegalArgumentException
        final var wrongX = new Matrix(nPoints - 1, 2);
        final var wrongY = new double[nPoints - 1];

        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(wrongX, y, 1.0));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, wrongY, 1.0));
    }

    @Test
    void testIsReady() throws FittingException, WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var fitter = new SvdMultiDimensionLinearFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final var x = new Matrix(nPoints, 2);
        final var y = new double[nPoints];
        final var sig = new double[nPoints];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final var evaluator = new LinearFitterMultiDimensionFunctionEvaluator() {

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
                // no action needed
            }
        };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());

        // test bad evaluator
        final var badEvaluator = new LinearFitterMultiDimensionFunctionEvaluator() {

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
                // no action needed
            }
        };

        fitter.setFunctionEvaluator(badEvaluator);

        assertFalse(fitter.isReady());
    }

    @Test
    void testGetSetTol() {
        final var fitter = new SvdMultiDimensionLinearFitter();

        // check default values
        assertEquals(SvdMultiDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // set new value
        fitter.setTol(1e-3);

        // check correctness
        assertEquals(1e-3, fitter.getTol(), 0.0);
    }

    @Test
    void testFitQuadratic() throws FittingException, NotReadyException, WrongSizeException {

        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {
            final var randomizer = new UniformRandomizer();

            final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

            final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

            final var params = new double[NUM_QUADRATIC_PARAMS];
            for (var i = 0; i < NUM_QUADRATIC_PARAMS; i++) {
                params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            }

            final var x = Matrix.createWithUniformRandomValues(npoints, 2, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var y = new double[npoints];
            final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
            for (var i = 0; i < npoints; i++) {
                final var x1 = x.getElementAt(i, 0);
                final var x2 = x.getElementAt(i, 1);
                // function is: a*1 + b*x + c*y + d*x*y + e*x^2 + f*y^2
                y[i] = params[0] + params[1] * x1 + params[2] * x2 + params[3] * x1 * x2 + params[4] * x1 * x1
                        + params[5] * x2 * x2;
                final var error = errorRandomizer.nextDouble();
                y[i] += error;
            }

            final var evaluator = new LinearFitterMultiDimensionFunctionEvaluator() {

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
                    final var x = point[0];
                    final var y = point[1];

                    result[0] = 1.0;
                    result[1] = x;
                    result[2] = y;
                    result[3] = x * y;
                    result[4] = x * x;
                    result[5] = y * y;
                }
            };

            final var fitter = new SvdMultiDimensionLinearFitter(evaluator, x, y, 1.0);

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

            var failed = false;
            for (var i = 0; i < NUM_QUADRATIC_PARAMS; i++) {
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
