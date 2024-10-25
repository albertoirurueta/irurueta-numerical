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

class StraightLineFitterTest {
    private static final double MIN_PARAM_VALUE = -100.0;
    private static final double MAX_PARAM_VALUE = 100.0;

    private static final double MIN_SIGMA_VALUE = 1e-3;
    private static final double MAX_SIGMA_VALUE = 3.0;

    private static final int MIN_POINTS = 1000;
    private static final int MAX_POINTS = 10000;

    private static final double MIN_DATA_VALUE = -100.0;
    private static final double MAX_DATA_VALUE = 100.0;

    private static final double ABSOLUTE_ERROR = 1e-1;

    @Test
    void testConstructor() {
        // test constructor without arguments
        var fitter = new StraightLineFitter();

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());
        assertEquals(0.0, fitter.getA(), 0.0);
        assertEquals(0.0, fitter.getB(), 0.0);
        assertEquals(0.0, fitter.getSigA(), 0.0);
        assertEquals(0.0, fitter.getSigB(), 0.0);
        assertEquals(0.0, fitter.getChi2(), 0.0);
        assertEquals(1.0, fitter.getQ(), 0.0);
        assertEquals(0.0, fitter.getSigdat(), 0.0);

        // test constructor with input data
        final var x = new double[2];
        final var y = new double[2];

        fitter = new StraightLineFitter(x, y);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(0.0, fitter.getA(), 0.0);
        assertEquals(0.0, fitter.getB(), 0.0);
        assertEquals(0.0, fitter.getSigA(), 0.0);
        assertEquals(0.0, fitter.getSigB(), 0.0);
        assertEquals(0.0, fitter.getChi2(), 0.0);
        assertEquals(1.0, fitter.getQ(), 0.0);
        assertEquals(0.0, fitter.getSigdat(), 0.0);

        // Force IllegalArgumentException
        final var shortX = new double[1];

        assertThrows(IllegalArgumentException.class, () -> new StraightLineFitter(shortX, y));

        // test constructor with input data and sigmas
        final var sig = new double[2];
        fitter = new StraightLineFitter(x, y, sig);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(0.0, fitter.getA(), 0.0);
        assertEquals(0.0, fitter.getB(), 0.0);
        assertEquals(0.0, fitter.getSigA(), 0.0);
        assertEquals(0.0, fitter.getSigB(), 0.0);
        assertEquals(0.0, fitter.getChi2(), 0.0);
        assertEquals(1.0, fitter.getQ(), 0.0);
        assertEquals(0.0, fitter.getSigdat(), 0.0);

        // Force IllegalArgumentException
        final var shortSig = new double[1];
        assertThrows(IllegalArgumentException.class, () -> new StraightLineFitter(x, y, shortSig));
    }

    @Test
    void testGetSetInputDataAndIsReady() {
        final var fitter = new StraightLineFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());

        // set new values
        final var x = new double[2];
        final var y = new double[2];

        fitter.setInputData(x, y);

        // check correctness
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());

        // Force IllegalArgumentException
        final var shortX = new double[1];
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputData(shortX, y));
    }

    @Test
    void testGetSetInputDataAndStandardDeviations() {
        final var fitter = new StraightLineFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());

        // set new values
        final var x = new double[2];
        final var y = new double[2];
        final var sig = new double[2];

        fitter.setInputDataAndStandardDeviations(x, y, sig);

        // check correctness
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertTrue(fitter.isReady());

        fitter.setInputDataAndStandardDeviations(x, y, null);

        // check correctness
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());

        // Force IllegalArgumentException
        final var shortX = new double[1];
        final var shortSig = new double[1];
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputDataAndStandardDeviations(shortX, y, sig));
        assertThrows(IllegalArgumentException.class, () -> fitter.setInputDataAndStandardDeviations(x, y, shortSig));
    }

    @Test
    void testFitNoSig() throws FittingException, NotReadyException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final var b = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);

        // generate data
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var x = new double[nPoints];
        final var y = new double[nPoints];
        for (var i = 0; i < nPoints; i++) {
            x[i] = randomizer.nextDouble(MIN_DATA_VALUE, MAX_DATA_VALUE);
            y[i] = a + b * x[i] + errorRandomizer.nextDouble();
        }

        var fitter = new StraightLineFitter(x, y);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(0.0, fitter.getA(), 0.0);
        assertEquals(0.0, fitter.getB(), 0.0);
        assertEquals(0.0, fitter.getSigA(), 0.0);
        assertEquals(0.0, fitter.getSigB(), 0.0);
        assertEquals(0.0, fitter.getChi2(), 0.0);
        assertEquals(1.0, fitter.getQ(), 0.0);
        assertEquals(0.0, fitter.getSigdat(), 0.0);

        // fit data
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertEquals(a, fitter.getA(), ABSOLUTE_ERROR);
        assertEquals(b, fitter.getB(), ABSOLUTE_ERROR);
        assertEquals(sigma, fitter.getSigdat(), ABSOLUTE_ERROR);
        assertTrue(fitter.getChi2() > 0);
        assertTrue(fitter.getQ() > 0);
        assertTrue(fitter.getSigA() > 0);
        assertEquals(fitter.getSigA(), fitter.getSigB() * Math.abs(fitter.getA()), ABSOLUTE_ERROR);
        assertTrue(fitter.getSigB() > 0);

        // Force NotReadyException
        final var fitter2 = new StraightLineFitter();

        assertFalse(fitter2.isReady());
        assertThrows(NotReadyException.class, fitter2::fit);
    }

    @Test
    void testFitSig() throws FittingException, NotReadyException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final var b = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final var sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final var errorRandomizer = new GaussianRandomizer(0.0, sigma);

        // generate data
        final var nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final var x = new double[nPoints];
        final var y = new double[nPoints];
        final var sig = new double[nPoints];
        for (var i = 0; i < nPoints; i++) {
            x[i] = randomizer.nextDouble(MIN_DATA_VALUE, MAX_DATA_VALUE);
            y[i] = a + b * x[i] + errorRandomizer.nextDouble();
            sig[i] = sigma;
        }

        var fitter = new StraightLineFitter(x, y, sig);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(x, fitter.getX());
        assertSame(y, fitter.getY());
        assertSame(sig, fitter.getSig());
        assertTrue(fitter.isReady());
        assertEquals(0.0, fitter.getA(), 0.0);
        assertEquals(0.0, fitter.getB(), 0.0);
        assertEquals(0.0, fitter.getSigA(), 0.0);
        assertEquals(0.0, fitter.getSigB(), 0.0);
        assertEquals(0.0, fitter.getChi2(), 0.0);
        assertEquals(1.0, fitter.getQ(), 0.0);
        assertEquals(0.0, fitter.getSigdat(), 0.0);

        // fit data
        fitter.fit();

        // check correctness
        assertTrue(fitter.isResultAvailable());
        assertEquals(a, fitter.getA(), 2.0 * ABSOLUTE_ERROR);
        assertEquals(b, fitter.getB(), 2.0 * ABSOLUTE_ERROR);
        assertEquals(0.0, fitter.getSigdat(), 0.0);
        assertTrue(fitter.getChi2() > 0);
        assertTrue(fitter.getQ() > 0);
        assertTrue(fitter.getSigA() > 0);
        assertEquals(fitter.getSigA(), fitter.getSigB() * Math.abs(fitter.getA()), ABSOLUTE_ERROR);
        assertTrue(fitter.getSigB() > 0);

        // Force NotReadyException
        final var fitter2 = new StraightLineFitter();

        assertFalse(fitter2.isReady());
        assertThrows(NotReadyException.class, fitter2::fit);
    }
}
