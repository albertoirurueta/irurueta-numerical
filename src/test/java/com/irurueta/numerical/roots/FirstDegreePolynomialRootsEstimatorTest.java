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

class FirstDegreePolynomialRootsEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;

    private static final double TOLERANCE = 3e-8;

    private static final int TIMES = 100;

    @Test
    void testConstructor() throws NotAvailableException, NotReadyException {

        final var polyParams = new double[2];
        // to ensure it's second degree
        polyParams[1] = 1.0;

        // test 1st constructor
        var estimator = new FirstDegreePolynomialRootsEstimator();
        assertNotNull(estimator);

        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertThrows(NotReadyException.class, estimator::estimate);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::isFirstDegree);
        assertTrue(estimator.isRealSolution());

        // test 2nd constructor
        estimator = new FirstDegreePolynomialRootsEstimator(polyParams);
        assertNotNull(estimator);

        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(estimator.getRealPolynomialParameters(), polyParams);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertTrue(estimator.isFirstDegree());
        assertTrue(estimator.isReady());

        // Force IllegalArgumentException
        final var badPolyParams = new double[1];
        assertThrows(IllegalArgumentException.class, () -> new FirstDegreePolynomialRootsEstimator(badPolyParams));
    }

    @Test
    void testGetSetPolynomialParameters() throws LockedException, NotAvailableException {

        final var polyParams = new double[2];
        polyParams[1] = 1.0;
        final var polyParams2 = new double[3];
        polyParams2[1] = 1.0;
        polyParams2[2] = 0.0;
        final var badPolyParams = new double[1];
        final var badPolyParams2 = new double[2];
        badPolyParams2[1] = 0.0;

        final var estimator = new FirstDegreePolynomialRootsEstimator();

        // check default values
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRealPolynomialParameters);
        assertFalse(estimator.arePolynomialParametersAvailable());

        // set polynomial parameters
        estimator.setPolynomialParameters(polyParams);
        // check correctness
        assertEquals(polyParams, estimator.getRealPolynomialParameters());
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);

        estimator.setPolynomialParameters(polyParams2);
        // check correctness
        assertEquals(estimator.getRealPolynomialParameters(), polyParams2);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams));
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams2));

        // attempting to use complex parameters will also raise an
        // IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(new Complex[3]));
    }

    @Test
    void testEstimate() throws LockedException, NotReadyException, NotAvailableException {

        for (var t = 0; t < TIMES; t++) {

            final var randomizer = new UniformRandomizer();
            final var realRoot1 = randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT);

            final var root1 = new Complex();
            root1.setReal(randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT));
            root1.setImaginary(randomizer.nextDouble(MIN_EVAL_POINT, MAX_EVAL_POINT));

            final var estimator = new FirstDegreePolynomialRootsEstimator();

            final Complex[] roots;

            // attempt set parameters for constant
            final var polyParams = generateConstantPolynomialParams(realRoot1);
            assertFalse(FirstDegreePolynomialRootsEstimator.isFirstDegree(
                    polyParams));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams));
            assertThrows(NotReadyException.class, estimator::isFirstDegree);
            assertThrows(NotReadyException.class, estimator::estimate);

            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for first degree polynomial with real root
            final var polyParams2 = generateFirstDegreePolynomialParams(realRoot1);
            assertTrue(FirstDegreePolynomialRootsEstimator.isFirstDegree(polyParams2));
            estimator.setPolynomialParameters(polyParams2);
            assertTrue(estimator.isFirstDegree());
            assertTrue(estimator.isRealSolution());
            estimator.estimate();

            // check correctness
            assertTrue(estimator.areRootsAvailable());
            roots = estimator.getRoots();

            assertEquals(1, roots.length);
            assertTrue(roots[0].equals(new Complex(realRoot1), TOLERANCE));
        }
    }

    private static double vectorNorm(double[] v) {
        var normValue = 0.0;
        for (final var value : v) {
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
}
