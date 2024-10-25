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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class SvdSingleDimensionLinearFitterTest {
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

    @Test
    void testConstructor() throws FittingException {
        // test empty constructor
        var fitter = new SvdSingleDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // test constructor with input data
        final var x = new double[2];
        final var y = new double[2];
        final var sig = new double[2];

        fitter = new SvdSingleDimensionLinearFitter(x, y, sig);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        final var shortX = new double[1];
        final var shortY = new double[1];
        final var shortSig = new double[1];

        assertThrows(IllegalArgumentException.class, () -> new SvdSingleDimensionLinearFitter(shortX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> new SvdSingleDimensionLinearFitter(x, shortY, sig));
        assertThrows(IllegalArgumentException.class, () -> new SvdSingleDimensionLinearFitter(x, y, shortSig));

        // test constructor with input data (constant sigma)
        fitter = new SvdSingleDimensionLinearFitter(x, y, 1.0);

        // check default values
        assertNull(fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertFalse(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new SvdSingleDimensionLinearFitter(shortX, y, 1.0));
        assertThrows(IllegalArgumentException.class, () -> new SvdSingleDimensionLinearFitter(x, shortY, 1.0));

        // test constructor with evaluator
        final var evaluator = new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[2];
            }

            @Override
            public void evaluate(final double point, final double[] result) {
                // no action needed
            }
        };

        fitter = new SvdSingleDimensionLinearFitter(evaluator);

        // check default values
        assertSame(evaluator, fitter.getFunctionEvaluator());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // test constructor with evaluator and input data
        fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, sig);

        // check default values
        assertSame(evaluator, fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class,
                () -> new SvdSingleDimensionLinearFitter(evaluator, shortX, y, sig));
        assertThrows(IllegalArgumentException.class,
                () -> new SvdSingleDimensionLinearFitter(evaluator, x, shortY, sig));
        assertThrows(IllegalArgumentException.class,
                () -> new SvdSingleDimensionLinearFitter(evaluator, x, y, shortSig));

        // test constructor with evaluator and input data (constant sigma)
        fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, 1.0);

        // check default values
        assertSame(evaluator, fitter.getFunctionEvaluator());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNotNull(fitter.getSig());
        for (var i = 0; i < fitter.getSig().length; i++) {
            assertEquals(1.0, fitter.getSig()[i], 0.0);
        }
        assertTrue(fitter.isReady());
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class,
                () -> new SvdSingleDimensionLinearFitter(evaluator, shortX, y, 1.0));
        assertThrows(IllegalArgumentException.class,
                () -> new SvdSingleDimensionLinearFitter(evaluator, x, shortY, 1.0));
    }

    @Test
    void testGetSetFunctionEvaluator() throws FittingException {
        final var fitter = new SvdSingleDimensionLinearFitter();

        // check default values
        assertNull(fitter.getFunctionEvaluator());

        // set new value
        final var evaluator = new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[2];
            }

            @Override
            public void evaluate(final double point, final double[] result) {
                // no action needed
            }
        };

        // set new value
        fitter.setFunctionEvaluator(evaluator);

        // check correctness
        assertSame(evaluator, fitter.getFunctionEvaluator());
    }

    @Test
    void testGetSetInputData() {
        final var fitter = new SvdSingleDimensionLinearFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final var x = new double[2];
        final var y = new double[2];
        final var sig = new double[2];

        // set input data
        fitter.setInputData(x, y, sig);

        // check correctness
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());

        // Force IllegalArgumentException
        final var wrong = new double[1];
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(wrong, y, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, wrong, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, y, wrong));
    }

    @Test
    void testGetSetInputDataWithConstantSigma() {
        final var fitter = new SvdSingleDimensionLinearFitter();

        // check default value
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());

        final var x = new double[2];
        final var y = new double[2];

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
        final var wrong = new double[1];
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(wrong, y, 1.0));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(x, wrong, 1.0));
    }

    @Test
    void testIsReady() throws FittingException {
        final var fitter = new SvdSingleDimensionLinearFitter();

        // check default value
        assertFalse(fitter.isReady());

        // set new values
        final var x = new double[2];
        final var y = new double[2];
        final var sig = new double[2];

        fitter.setInputData(x, y, sig);

        assertFalse(fitter.isReady());

        final var evaluator = new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[2];
            }

            @Override
            public void evaluate(final double point, final double[] result) {
                // no action needed
            }
        };

        fitter.setFunctionEvaluator(evaluator);

        assertTrue(fitter.isReady());
    }

    @Test
    void testGetSetTol() {
        final var fitter = new SvdSingleDimensionLinearFitter();

        // check default values
        assertEquals(SvdSingleDimensionLinearFitter.DEFAULT_TOL, fitter.getTol(), 0.0);

        // set new value
        fitter.setTol(1e-3);

        // check correctness
        assertEquals(1e-3, fitter.getTol(), 0.0);
    }

    @Test
    void testFitPolynomial() throws FittingException, NotReadyException {
        final var randomizer = new UniformRandomizer();

        final var nparams = randomizer.nextInt(MIN_POLY_PARAMS, MAX_POLY_PARAMS);
        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[nparams];
        for (var i = 0; i < nparams; i++) {
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        }

        final var y = new double[npoints];
        final var x = new double[npoints];
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        for (var i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            y[i] = 0.0;
            for (var j = 0; j < nparams; j++) {
                y[i] += params[j] * Math.pow(x[i], j);
            }
            final var error = errorRandomizer.nextDouble();
            y[i] += error;
        }

        final var evaluator = new LinearFitterSingleDimensionFunctionEvaluator() {

            @Override
            public double[] createResultArray() {
                return new double[nparams];
            }

            @Override
            public void evaluate(final double point, final double[] result) {
                for (var i = 0; i < result.length; i++) {
                    result[i] = Math.pow(point, i);
                }
            }
        };

        final var fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, 1.0);

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
        assertEquals(nparams, fitter.getA().length);
        for (var i = 0; i < nparams; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);
    }

    @Test
    void testTrigo() throws FittingException, NotReadyException {
        final var randomizer = new UniformRandomizer();

        final var npoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);

        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var params = new double[TRIGO_PARAMS];
        for (var i = 0; i < TRIGO_PARAMS; i++) {
            params[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        }

        final var y = new double[npoints];
        final var x = new double[npoints];
        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);
        // function: y(x) = a * cos(x) + b * sin(x)
        for (var i = 0; i < npoints; i++) {
            x[i] = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
            final var error = errorRandomizer.nextDouble();
            y[i] = params[0] * Math.sin(x[i]) + params[1] * Math.cos(x[i]) + error;
        }

        final var evaluator = new LinearFitterSingleDimensionFunctionEvaluator() {

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

        final var fitter = new SvdSingleDimensionLinearFitter(evaluator, x, y, 1.0);

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
        assertEquals(TRIGO_PARAMS, fitter.getA().length);
        for (var i = 0; i < TRIGO_PARAMS; i++) {
            assertEquals(fitter.getA()[i], params[i], ABSOLUTE_ERROR);
        }
        assertNotNull(fitter.getCovar());
        assertTrue(fitter.getChisq() > 0);
    }
}
