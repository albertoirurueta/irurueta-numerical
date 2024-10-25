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

class ThirdDegreePolynomialRootsEstimatorTest {

    private static final double MIN_EVAL_POINT = 0.0;
    private static final double MAX_EVAL_POINT = 1.0;

    private static final double TOLERANCE = 3e-7;

    private static final int TIMES = 100;

    @Test
    void testConstructor() throws NotAvailableException, NotReadyException {

        final var polyParams = new double[4];
        // to ensure it's third degree
        polyParams[3] = 1.0;

        // test 1st constructor
        var estimator = new ThirdDegreePolynomialRootsEstimator();
        assertNotNull(estimator);

        assertFalse(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertThrows(NotReadyException.class, estimator::estimate);
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertFalse(estimator.isReady());
        assertThrows(NotReadyException.class, estimator::isThirdDegree);
        assertThrows(NotReadyException.class, estimator::hasThreeDistinctRealRoots);
        assertThrows(NotReadyException.class, estimator::hasMultipleRealRoot);
        assertThrows(NotReadyException.class, estimator::hasOneRealRootAndTwoComplexConjugateRoots);

        // test 2nd constructor
        estimator = new ThirdDegreePolynomialRootsEstimator(polyParams);
        assertNotNull(estimator);

        assertTrue(estimator.arePolynomialParametersAvailable());
        assertFalse(estimator.areRootsAvailable());
        assertEquals(polyParams, estimator.getRealPolynomialParameters());
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);
        assertThrows(NotAvailableException.class, estimator::getRoots);
        assertFalse(estimator.isLocked());
        assertTrue(estimator.isReady());
        assertTrue(estimator.isThirdDegree());
        estimator.hasThreeDistinctRealRoots();
        estimator.hasMultipleRealRoot();
        estimator.hasOneRealRootAndTwoComplexConjugateRoots();

        // Force IllegalArgumentException
        final var badPolyParams = new double[1];

        assertThrows(IllegalArgumentException.class, () -> new ThirdDegreePolynomialRootsEstimator(badPolyParams));
    }

    @Test
    void testGetSetPolynomialParameters() throws LockedException, NotAvailableException {

        final var polyParams = new double[4];
        polyParams[3] = 1.0;
        final var polyParams2 = new double[5];
        polyParams2[3] = 1.0;
        polyParams2[4] = 0.0;
        final var badPolyParams = new double[3];
        final var badPolyParams2 = new double[4];
        badPolyParams2[3] = 0.0;

        final var estimator = new ThirdDegreePolynomialRootsEstimator();

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
        assertEquals(polyParams2, estimator.getRealPolynomialParameters());
        assertThrows(NotAvailableException.class, estimator::getPolynomialParameters);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams));

        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(badPolyParams2));

        // attempting to use complex parameters will also raise an
        // IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(new Complex[4]));
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

            final var estimator = new ThirdDegreePolynomialRootsEstimator();

            Complex[] roots;

            // attempt set parameters for constant
            final var polyParams = generateConstantPolynomialParams(realRoot1);
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams));
            assertThrows(NotReadyException.class, estimator::isThirdDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // attempt estimate root for first degree polynomial
            final var polyParams2 = generateFirstDegreePolynomialParams(realRoot1);
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams2));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams2));
            assertThrows(NotReadyException.class, estimator::isThirdDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for second degree polynomial with real roots
            final var polyParams3 = generateSecondDegreePolynomialParams(new Complex(realRoot1),
                    new Complex(realRoot2));
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams3));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams3));
            assertThrows(NotReadyException.class, estimator::isThirdDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for second degree polynomial with complex conjugate
            // roots (and real coefficients)
            final var polyParams4 = generateSecondDegreePolynomialParams(root1, conjRoot1);
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams4));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams4));
            assertThrows(NotReadyException.class, estimator::isThirdDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for second degree polynomial with double real roots
            final var polyParams5 = generateSecondDegreePolynomialParams(
                    new Complex(realRoot1), new Complex(realRoot1));
            assertFalse(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams5));
            assertThrows(IllegalArgumentException.class, () -> estimator.setPolynomialParameters(polyParams5));
            assertThrows(NotReadyException.class, estimator::isThirdDegree);
            assertThrows(NotReadyException.class, estimator::estimate);
            // check correctness
            assertFalse(estimator.areRootsAvailable());

            // set parameters for third degree polynomial with real roots
            final var polyParams6 = generateThirdDegreePolynomialParams(new Complex(realRoot1), new Complex(realRoot2),
                    new Complex(realRoot3));
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams6));
            estimator.setPolynomialParameters(polyParams6);
            assertTrue(estimator.isThirdDegree());
            assertTrue(estimator.hasThreeDistinctRealRoots());
            assertFalse(estimator.hasMultipleRealRoot());
            assertFalse(estimator.hasOneRealRootAndTwoComplexConjugateRoots());
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
            final var polyParams7 = generateThirdDegreePolynomialParams(new Complex(realRoot1), root2, conjRoot2);
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams7));
            estimator.setPolynomialParameters(polyParams7);
            assertTrue(estimator.isThirdDegree());
            assertFalse(estimator.hasThreeDistinctRealRoots());
            assertFalse(estimator.hasMultipleRealRoot());
            assertTrue(estimator.hasOneRealRootAndTwoComplexConjugateRoots());
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
            final var polyParams8 = generateThirdDegreePolynomialParams(new Complex(realRoot1), new Complex(realRoot2),
                    new Complex(realRoot2));
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams8));
            estimator.setPolynomialParameters(polyParams8);
            assertTrue(estimator.isThirdDegree());
            assertFalse(estimator.hasThreeDistinctRealRoots());
            assertTrue(estimator.hasMultipleRealRoot());
            assertFalse(estimator.hasOneRealRootAndTwoComplexConjugateRoots());
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
            final var polyParams9 = generateThirdDegreePolynomialParams(new Complex(realRoot1), new Complex(realRoot1),
                    new Complex(realRoot1));
            assertTrue(ThirdDegreePolynomialRootsEstimator.isThirdDegree(polyParams9));
            estimator.setPolynomialParameters(polyParams9);
            assertTrue(estimator.isThirdDegree());
            assertFalse(estimator.hasThreeDistinctRealRoots());
            assertTrue(estimator.hasMultipleRealRoot());
            assertFalse(estimator.hasOneRealRootAndTwoComplexConjugateRoots());
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

    private static double[] generateSecondDegreePolynomialParams(final Complex root1, final Complex root2) {
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

    private static double[] generateThirdDegreePolynomialParams(final Complex root1, final Complex root2,
                                                                final Complex root3) {

        final var out = new double[4];
        // p(x) = (x - root1) * (x - root2) * (x - root3) =
        // (x * x - (root1 + root2) * x + root1 * root2) * (x - root3) =
        // (x * x * x - (root1 + root2) * x * x + (root1 + root2) * x
        // - root3 * x * x + (root1 + root2) * root3 * x
        // - (root1 + root2) * root3 =

        // x * x * x - (root1 + root2 + root3) * x * x +
        // ((root1 * root2) + (root1 + root2) * root3) * x
        // - root1 * root2 * root3

        out[3] = 1.0;
        out[2] = root1.addAndReturnNew(root2).addAndReturnNew(root3).multiplyByScalarAndReturnNew(-1.0).getReal();
        out[1] = root1.multiplyAndReturnNew(root2).addAndReturnNew(
                root1.addAndReturnNew(root2).multiplyAndReturnNew(root3)).getReal();
        out[0] = root1.multiplyAndReturnNew(root2).multiplyAndReturnNew(root3).multiplyByScalarAndReturnNew(-1.0)
                .getReal();

        final var normValue = vectorNorm(out);
        // normalize vector of complex values
        ArrayUtils.multiplyByScalar(out, 1.0 / normValue, out);
        return out;
    }
}
