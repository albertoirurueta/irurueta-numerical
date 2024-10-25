/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.roots;

import com.irurueta.algebra.ArrayUtils;
import com.irurueta.algebra.Complex;
import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotAvailableException;
import com.irurueta.numerical.NotReadyException;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class LaguerrePolynomialRootsEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;

    private static final double TOLERANCE = 3e-8;

    private static final int TIMES = 100;

    @Test
    void testConstructor() throws NotAvailableException {

        final var randomizer = new UniformRandomizer();
        final var polishRoots = randomizer.nextBoolean();

        final var polyParams = new Complex[2];
        final var badPolyParams = new Complex[1];

        // test 1st constructor
        var estimator = new LaguerrePolynomialRootsEstimator();
        assertNotNull(estimator);

        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(LaguerrePolynomialRootsEstimator.DEFAULT_POLISH_ROOTS, estimator.areRootsPolished());
        assertThrows(NotReadyException.class, estimator::estimate);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());

        // test 2nd constructor
        estimator = new LaguerrePolynomialRootsEstimator(polishRoots);
        assertNotNull(estimator);
        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.areRootsPolished(), polishRoots);
        assertThrows(NotReadyException.class, estimator::estimate);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());

        // test 3rd constructor
        estimator = new LaguerrePolynomialRootsEstimator(polyParams);
        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(LaguerrePolynomialRootsEstimator.DEFAULT_POLISH_ROOTS, estimator.areRootsPolished());
        assertSame(polyParams, estimator.getPolynomialParameters());
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LaguerrePolynomialRootsEstimator(badPolyParams));

        // test 4th constructor
        estimator = new LaguerrePolynomialRootsEstimator(polyParams, polishRoots);
        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(polishRoots, estimator.areRootsPolished());
        assertSame(polyParams, estimator.getPolynomialParameters());
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new LaguerrePolynomialRootsEstimator(badPolyParams,
                polishRoots));
    }

    @Test
    void testGetSetPolishRoots() throws LockedException {

        final var randomizer = new UniformRandomizer();
        final var polishRoots = randomizer.nextBoolean();

        final var estimator = new LaguerrePolynomialRootsEstimator();

        // check default value
        assertEquals(LaguerrePolynomialRootsEstimator.DEFAULT_POLISH_ROOTS, estimator.areRootsPolished());

        // set new value
        estimator.setPolishRootsEnabled(polishRoots);
        // check correctness
        assertEquals(estimator.areRootsPolished(), polishRoots);
    }

    @Test
    void testGetSetPolynomialParameters() throws LockedException, NotAvailableException {

        final var polyParams = new Complex[2];
        final var badPolyParams = new Complex[1];

        final var estimator = new LaguerrePolynomialRootsEstimator();

        // check default values
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertFalse(estimator.arePolynomialParametersAvailable());

        // set polynomial parameters
        estimator.setPolynomialParameters(polyParams);
        // check correctness
        assertSame(polyParams, estimator.getPolynomialParameters());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams));
    }

    @Test
    void testEstimate() throws LockedException, NotReadyException, RootEstimationException, NotAvailableException {

        for (var t = 0; t < TIMES; t++) {

            final var randomizer = new UniformRandomizer();
            final var realRoot1 = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT / 3.0);
            final var realRoot2 = randomizer.nextDouble(MAX_EVAL_POINT / 3.0, 2.0 / 3.0 * MAX_EVAL_POINT);
            final var realRoot3 = randomizer.nextDouble(2.0 / 3.0 * MAX_EVAL_POINT, MAX_EVAL_POINT);

            final var root1 = new Complex();
            root1.setReal(randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT / 3.0));
            root1.setImaginary(randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT / 3.0));

            final var root2 = new Complex();
            root2.setReal(randomizer.nextDouble(MAX_EVAL_POINT / 3.0, 2.0 / 3.0 * MAX_EVAL_POINT));
            root2.setImaginary(randomizer.nextDouble(MAX_EVAL_POINT / 3.0, 2.0 / 3.0 * MAX_EVAL_POINT));

            final var root3 = new Complex();
            root3.setReal(randomizer.nextDouble(2.0 / 3.0 * MAX_EVAL_POINT, MAX_EVAL_POINT));
            root3.setImaginary(randomizer.nextDouble(2.0 / 3.0 * MAX_EVAL_POINT, MAX_EVAL_POINT));

            // also compute conjugates
            final var conjRoot1 = root1.conjugateAndReturnNew();
            final var conjRoot2 = root2.conjugateAndReturnNew();

            final var estimator = new LaguerrePolynomialRootsEstimator();

            // attempt set parameters for constant
            final var polyParams = generateConstantPolynomialParams(new Complex(realRoot1));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams));

            // set parameters for first degree polynomial with real root
            final var polyParams2 = generateFirstDegreePolynomialParams(new Complex(realRoot1));
            estimator.setPolynomialParameters(polyParams2);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            var roots = estimator.getRoots();

            assertEquals(1, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));

            // set parameters for first degree polynomial with complex root
            final var polyParams3 = generateFirstDegreePolynomialParams(root1);
            estimator.setPolynomialParameters(polyParams3);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(1, roots.length);
            assertTrue(roots[0].equals(root1, TOLERANCE));

            // set parameters for second degree polynomial with real roots
            final var polyParams4 = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2));
            estimator.setPolynomialParameters(polyParams4);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(2, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));

            // set parameters for second degree polynomial with complex conjugate
            // roots (and real coefficients)
            final var polyParams5 = generateSecondDegreePolynomialParams(root1, conjRoot1);
            estimator.setPolynomialParameters(polyParams5);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(2, roots.length);
            // because root[0] and root[1] might be exchanged, we check for their
            // real parts and absolute value of their imaginary parts (which are
            // the same but with opposite sign because they are complex
            // conjugates)
            assertEquals(roots[0].getReal(), root1.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[0].getImaginary()), Math.abs(root1.getImaginary()), TOLERANCE);
            assertEquals(roots[1].getReal(), conjRoot1.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[1].getImaginary()), Math.abs(conjRoot1.getImaginary()), TOLERANCE);

            // set parameters for second degree polynomial with double real roots
            final var polyParams6 = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1));
            estimator.setPolynomialParameters(polyParams6);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(2, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot1), TOLERANCE));

            // set parameters for third degree polynomial with real roots
            final var polyParams7 = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2), new Complex(realRoot3));
            estimator.setPolynomialParameters(polyParams7);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(3, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));
            assertTrue(roots[2].equals(new Complex(realRoot3), TOLERANCE));

            // set parameters for third degree polynomial with real root and two
            // complex conjugate roots
            final Complex[] polyParams8;
            if (realRoot1 < root2.getReal()) {
                polyParams8 = generateThirdDegreePolynomialParams(new Complex(realRoot1), root2, conjRoot2);
            } else {
                polyParams8 = generateThirdDegreePolynomialParams(root2, conjRoot2, new Complex(realRoot1));
            }
            estimator.setPolynomialParameters(polyParams8);
            estimator.estimate();

            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(3, roots.length);
            // because roots1 might be exchanged, we check for their
            // real parts and absolute value of their imaginary parts (which are
            // the same but with opposite sign because they are complex
            // conjugates)
            assertEquals(roots[0].getReal(), realRoot1, TOLERANCE);
            assertEquals(0.0, Math.abs(roots[0].getImaginary()), TOLERANCE);
            assertEquals(roots[1].getReal(), root2.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[1].getImaginary()), Math.abs(root2.getImaginary()), TOLERANCE);
            assertEquals(roots[2].getReal(), conjRoot2.getReal(), TOLERANCE);
            assertEquals(Math.abs(roots[2].getImaginary()), Math.abs(conjRoot2.getImaginary()), TOLERANCE);

            // set parameters for third degree polynomial with two double real
            // roots
            final var polyParams9 = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2), new Complex(realRoot2));
            estimator.setPolynomialParameters(polyParams9);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(3, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));
            assertTrue(roots[2].equals(new Complex(realRoot2), TOLERANCE));

            // set parameters for third degree polynomial with one triple real
            // roots
            final var polyParams10 = generateThirdDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1), new Complex(realRoot1));
            estimator.setPolynomialParameters(polyParams10);
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(3, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[2].equals(new Complex(realRoot1), TOLERANCE));
        }
    }

    private static double vectorNorm(final Complex[] v) {
        var normValue = 0.0;
        for (final Complex value : v) {
            final var real = value.getReal();
            final var imag = value.getImaginary();
            // square norm
            normValue += real * real + imag * imag;
        }

        return Math.sqrt(normValue);
    }

    private static Complex[] generateConstantPolynomialParams(final Complex param) {

        final var out = new Complex[1];
        out[0] = param;
        return out;
    }

    private static Complex[] generateFirstDegreePolynomialParams(final Complex root1) {

        final var out = new Complex[2];
        // p(x) = x - root1
        out[1] = new Complex(1.0, 0.0);
        out[0] = new Complex(-root1.getReal(), -root1.getImaginary());

        final var normValue = vectorNorm(out);
        // normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, normValue, out);

        return out;
    }

    private Complex[] generateSecondDegreePolynomialParams(final Complex root1, final Complex root2) {

        final var out = new Complex[3];
        // p(x) = (x - root1) * (x - root2) = x * x - (root1 * root2) * x +
        // root1 + root2
        out[2] = new Complex(1.0, 0.0);
        out[1] = root1.addAndReturnNew(root2).multiplyByScalarAndReturnNew(-1.0);
        out[0] = root1.multiplyAndReturnNew(root2);

        final var normValue = vectorNorm(out);
        // normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, normValue, out);

        return out;
    }

    private Complex[] generateThirdDegreePolynomialParams(
            final Complex root1, final Complex root2, final Complex root3) {

        final var out = new Complex[4];
        // p(x) = (x - root1) * (x - root2) * (x - root3) =
        // (x * x - (root1 + root2) * x + root1 * root2) * (x - root3) =
        // (x * x * x - (root1 + root2) * x * x + (root1 + root2) * x
        // - root3 * x * x + (root1 + root2) * root3 * x
        // - (root1 + root2) * root3 =

        // x * x * x - (root1 + root2 + root3) * x * x +
        // ((root1 * root2) + (root1 + root2) * root3) * x
        // - root1 * root2 * root3

        out[3] = new Complex(1.0, 0.0);
        out[2] = root1.addAndReturnNew(root2).addAndReturnNew(root3).multiplyByScalarAndReturnNew(-1.0);
        out[1] = root1.multiplyAndReturnNew(root2).addAndReturnNew(
                root1.addAndReturnNew(root2).multiplyAndReturnNew(root3));
        out[0] = root1.multiplyAndReturnNew(root2).multiplyAndReturnNew(root3).multiplyByScalarAndReturnNew(-1.0);

        final var normValue = vectorNorm(out);
        // normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, normValue, out);

        return out;
    }
}
