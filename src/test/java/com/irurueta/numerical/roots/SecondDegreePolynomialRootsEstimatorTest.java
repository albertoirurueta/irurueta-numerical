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

class SecondDegreePolynomialRootsEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;

    private static final double TOLERANCE = 3e-8;

    private static final int TIMES = 100;

    @Test
    void testConstructor() throws NotAvailableException, NotReadyException {

        final var polyParams = new double[3];
        // to ensure it's second degree
        polyParams[2] = 1.0;

        // test 1st constructor
        var estimator = new SecondDegreePolynomialRootsEstimator();
        assertNotNull(estimator);

        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertThrows(NotReadyException.class, estimator::estimate);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::isSecondDegree);
        assertThrows(NotReadyException.class, estimator::hasTwoDistinctRealRoots);
        assertThrows(NotReadyException.class, estimator::hasDoubleRoot);
        assertThrows(NotReadyException.class, estimator::hasTwoComplexConjugateRoots);

        // test 2nd constructor
        estimator = new SecondDegreePolynomialRootsEstimator(polyParams);
        assertNotNull(estimator);

        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertArrayEquals(estimator.getRealPolynomialParameters(), polyParams, TOLERANCE);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertTrue(estimator.isSecondDegree());
        estimator.hasTwoDistinctRealRoots();
        estimator.hasDoubleRoot();
        estimator.hasTwoComplexConjugateRoots();

        // Force IllegalArgumentException
        final var badPolyParams = new double[1];

        assertThrows(IllegalArgumentException.class, () -> new SecondDegreePolynomialRootsEstimator(badPolyParams));
    }

    @Test
    void testGetSetPolynomialParameters() throws LockedException, NotAvailableException {

        final var polyParams = new double[3];
        polyParams[2] = 1.0;
        final var polyParams2 = new double[4];
        polyParams2[2] = 1.0;
        polyParams2[3] = 0.0;
        final var badPolyParams = new double[2];
        final var badPolyParams2 = new double[3];
        badPolyParams2[2] = 0.0;

        final var estimator = new SecondDegreePolynomialRootsEstimator();

        // check default values
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRealPolynomialParameters);
        assertFalse(estimator.arePolynomialParametersAvailable());

        // set polynomial parameters
        estimator.setPolynomialParameters(polyParams);
        // check correctness
        assertArrayEquals(estimator.getRealPolynomialParameters(), polyParams, TOLERANCE);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);

        estimator.setPolynomialParameters(polyParams2);
        // check correctness
        assertArrayEquals(estimator.getRealPolynomialParameters(), polyParams2, TOLERANCE);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams));

        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams2));

        // attempting to use complex parameters will also raise an
        // IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(new Complex[3]));
    }

    @Test
    void testEstimate() throws LockedException, NotReadyException, RootEstimationException, NotAvailableException {

        for (var t = 0; t < TIMES; t++) {

            final var randomizer = new UniformRandomizer();
            final var realRoot1 = randomizer.nextDouble(MIN_EVAL_POINT, 0.5 * MAX_EVAL_POINT);
            final var realRoot2 = randomizer.nextDouble(0.5 * MAX_EVAL_POINT, MAX_EVAL_POINT);

            final var root1 = new Complex();
            root1.setReal(randomizer.nextDouble(MIN_EVAL_POINT, 0.5 * MAX_EVAL_POINT));
            root1.setImaginary(randomizer.nextDouble(MIN_EVAL_POINT, 0.5 * MAX_EVAL_POINT));

            final var root2 = new Complex();
            root2.setReal(randomizer.nextDouble(0.5 * MAX_EVAL_POINT, MAX_EVAL_POINT));
            root2.setImaginary(randomizer.nextDouble(0.5 * MAX_EVAL_POINT, MAX_EVAL_POINT));

            // also compute conjugates
            final var conjRoot1 = root1.conjugateAndReturnNew();

            final var estimator = new SecondDegreePolynomialRootsEstimator();

            // attempt set parameters for constant
            final var polyParams = generateConstantPolynomialParams(realRoot1);
            assertFalse(SecondDegreePolynomialRootsEstimator.isSecondDegree(polyParams));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams));
            assertThrows(NotReadyException.class, estimator::isSecondDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for first degree polynomial with real root
            final var polyParams2 = generateFirstDegreePolynomialParams(realRoot1);
            assertFalse(SecondDegreePolynomialRootsEstimator.isSecondDegree(polyParams2));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams));
            assertThrows(NotReadyException.class, estimator::isSecondDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for second degree polynomial with real roots
            final var polyParams3 = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot2));
            assertTrue(SecondDegreePolynomialRootsEstimator.isSecondDegree(polyParams3));
            estimator.setPolynomialParameters(polyParams3);
            assertTrue(estimator.isSecondDegree());
            assertTrue(estimator.hasTwoDistinctRealRoots());
            assertFalse(estimator.hasDoubleRoot());
            assertFalse(estimator.hasTwoComplexConjugateRoots());
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            var roots = estimator.getRoots();

            assertEquals(2, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot2), TOLERANCE));

            // set parameters for second degree polynomial with complex conjugate
            // roots (and real coefficients)
            final var polyParams4 = generateSecondDegreePolynomialParams(root1, conjRoot1);
            assertTrue(SecondDegreePolynomialRootsEstimator.isSecondDegree(polyParams4));
            estimator.setPolynomialParameters(polyParams4);
            assertTrue(estimator.isSecondDegree());
            assertFalse(estimator.hasTwoDistinctRealRoots());
            assertFalse(estimator.hasDoubleRoot());
            assertTrue(estimator.hasTwoComplexConjugateRoots());
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
            final var polyParams5 = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1));
            assertTrue(SecondDegreePolynomialRootsEstimator.isSecondDegree(polyParams5));
            estimator.setPolynomialParameters(polyParams5);
            assertTrue(estimator.isSecondDegree());
            assertFalse(estimator.hasTwoDistinctRealRoots());
            assertTrue(estimator.hasDoubleRoot());
            assertFalse(estimator.hasTwoComplexConjugateRoots());
            estimator.estimate();
            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(2, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
            assertTrue(roots[1].equals(new Complex(realRoot1), TOLERANCE));
        }
    }

    private static double vectorNorm(final double[] v) {
        var normValue = 0.0;
        for (var value : v) {
            // square norm
            normValue += value * value;
        }

        return Math.sqrt(normValue);
    }

    private static double[] generateConstantPolynomialParams(final double param) {
        final var out = new double[1];
        out[0] = param;
        return out;
    }

    private static double[] generateFirstDegreePolynomialParams(final double root1) {
        final var out = new double[2];
        // p(x) = x - root1
        out[1] = 1.0;
        out[0] = -root1;

        final var normValue = vectorNorm(out);
        // normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);
        return out;
    }

    private double[] generateSecondDegreePolynomialParams(final Complex root1, final Complex root2) {
        final var out = new double[3];
        // p(x) = (x - root1) * (x - root2) = x * x - (root1 + root2) * x +
        // root1 * root2
        out[2] = 1.0;
        out[1] = root1.addAndReturnNew(root2).multiplyByScalarAndReturnNew(-1.0).getReal();
        out[0] = root1.multiplyAndReturnNew(root2).getReal();

        final var normValue = vectorNorm(out);
        // normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);
        return out;
    }
}
