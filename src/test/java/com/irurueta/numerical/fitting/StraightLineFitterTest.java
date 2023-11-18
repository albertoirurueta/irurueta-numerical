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

public class StraightLineFitterTest {
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
    public void testConstructor() {
        // test constructor without arguments
        StraightLineFitter fitter = new StraightLineFitter();

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
        final double[] x = new double[2];
        final double[] y = new double[2];

        fitter = new StraightLineFitter(x, y);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
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
        final double[] shortX = new double[1];

        fitter = null;
        try {
            fitter = new StraightLineFitter(shortX, y);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);

        // test constructor with input data and sigmas
        final double[] sig = new double[2];
        fitter = new StraightLineFitter(x, y, sig);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());
        assertEquals(0.0, fitter.getA(), 0.0);
        assertEquals(0.0, fitter.getB(), 0.0);
        assertEquals(0.0, fitter.getSigA(), 0.0);
        assertEquals(0.0, fitter.getSigB(), 0.0);
        assertEquals(0.0, fitter.getChi2(), 0.0);
        assertEquals(1.0, fitter.getQ(), 0.0);
        assertEquals(0.0, fitter.getSigdat(), 0.0);

        // Force IllegalArgumentException
        final double[] shortSig = new double[1];

        fitter = null;
        try {
            fitter = new StraightLineFitter(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        assertNull(fitter);
    }

    @Test
    public void testGetSetInputDataAndIsReady() {
        final StraightLineFitter fitter = new StraightLineFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());

        // set new values
        final double[] x = new double[2];
        final double[] y = new double[2];

        fitter.setInputData(x, y);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());

        // Force IllegalArgumentException
        final double[] shortX = new double[1];

        try {
            fitter.setInputData(shortX, y);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetInputDataAndStandardDeviations() {
        final StraightLineFitter fitter = new StraightLineFitter();

        // check default values
        assertNull(fitter.getX());
        assertNull(fitter.getY());
        assertNull(fitter.getSig());
        assertFalse(fitter.isReady());

        // set new values
        final double[] x = new double[2];
        final double[] y = new double[2];
        final double[] sig = new double[2];

        fitter.setInputDataAndStandardDeviations(x, y, sig);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
        assertTrue(fitter.isReady());

        fitter.setInputDataAndStandardDeviations(x, y, null);

        // check correctness
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertNull(fitter.getSig());
        assertTrue(fitter.isReady());

        // Force IllegalArgumentException
        final double[] shortX = new double[1];
        final double[] shortSig = new double[1];
        try {
            fitter.setInputDataAndStandardDeviations(shortX, y, sig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
        try {
            fitter.setInputDataAndStandardDeviations(x, y, shortSig);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testFitNoSig() throws FittingException, NotReadyException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double a = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final double b = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);


        // generate data
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double[] x = new double[nPoints];
        final double[] y = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            x[i] = randomizer.nextDouble(MIN_DATA_VALUE, MAX_DATA_VALUE);
            y[i] = a + b * x[i] + errorRandomizer.nextDouble();
        }

        StraightLineFitter fitter = new StraightLineFitter(x, y);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
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
        assertEquals(fitter.getA(), a, ABSOLUTE_ERROR);
        assertEquals(fitter.getB(), b, ABSOLUTE_ERROR);
        assertEquals(fitter.getSigdat(), sigma, ABSOLUTE_ERROR);
        assertTrue(fitter.getChi2() > 0);
        assertTrue(fitter.getQ() > 0);
        assertTrue(fitter.getSigA() > 0);
        assertEquals(fitter.getSigA(),
                fitter.getSigB() * Math.abs(fitter.getA()), ABSOLUTE_ERROR);
        assertTrue(fitter.getSigB() > 0);

        // Force NotReadyException
        fitter = new StraightLineFitter();

        assertFalse(fitter.isReady());
        try {
            fitter.fit();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }
    }

    @Test
    public void testFitSig() throws FittingException, NotReadyException {
        final UniformRandomizer randomizer = new UniformRandomizer(new Random());
        final double a = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final double b = randomizer.nextDouble(MIN_PARAM_VALUE, MAX_PARAM_VALUE);
        final double sigma = randomizer.nextDouble(MIN_SIGMA_VALUE, MAX_SIGMA_VALUE);

        final GaussianRandomizer errorRandomizer = new GaussianRandomizer(
                new Random(), 0.0, sigma);


        // generate data
        final int nPoints = randomizer.nextInt(MIN_POINTS, MAX_POINTS);
        final double[] x = new double[nPoints];
        final double[] y = new double[nPoints];
        final double[] sig = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            x[i] = randomizer.nextDouble(MIN_DATA_VALUE, MAX_DATA_VALUE);
            y[i] = a + b * x[i] + errorRandomizer.nextDouble();
            sig[i] = sigma;
        }

        StraightLineFitter fitter = new StraightLineFitter(x, y, sig);

        // check default values
        assertFalse(fitter.isResultAvailable());
        assertSame(fitter.getX(), x);
        assertSame(fitter.getY(), y);
        assertSame(fitter.getSig(), sig);
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
        assertEquals(fitter.getA(), a, 2.0 * ABSOLUTE_ERROR);
        assertEquals(fitter.getB(), b, 2.0 * ABSOLUTE_ERROR);
        assertEquals(0.0, fitter.getSigdat(), 0.0);
        assertTrue(fitter.getChi2() > 0);
        assertTrue(fitter.getQ() > 0);
        assertTrue(fitter.getSigA() > 0);
        assertEquals(fitter.getSigA(),
                fitter.getSigB() * Math.abs(fitter.getA()), ABSOLUTE_ERROR);
        assertTrue(fitter.getSigB() > 0);

        // Force NotReadyException
        fitter = new StraightLineFitter();

        assertFalse(fitter.isReady());
        try {
            fitter.fit();
            fail("NotReadyException expected but not thrown");
        } catch (final NotReadyException ignore) {
        }
    }
}
