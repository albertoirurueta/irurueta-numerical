/*
 * Copyright (C) 2023 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.integration;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;

import org.junit.jupiter.api.Test;

class SimpsonLowerSquareRootMidPointQuadratureMatrixIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    /**
     * Fixed seed used to make randomized exponential integration tests deterministic.
     */
    private static final long SEED = 3_141_592_653L;

    /**
     * Number of randomized trials executed by exponential integration tests.
     */
    private static final int TIMES = 10;

    /**
     * Minimum length of the integration interval used in exponential integration tests.
     */
    private static final double MIN_INTERVAL = 0.5;

    /**
     * Minimum absolute value of sin(b) - sin(a) and cos(b) - cos(a) in exponential integration tests. The
     * integrator uses a relative convergence criterion based on the smallest absolute element of the result, which
     * cannot be met when any element of the expected integral is close to zero.
     */
    private static final double MIN_ABS_DIFF = 0.1;

    private static final double ABSOLUTE_ERROR_EXPONENTIAL = 1e-4;

    /**
     * Relative convergence accuracy used by the integrator in the exponential integration test. It is tighter than
     * {@link #ABSOLUTE_ERROR_EXPONENTIAL}, because convergence is only checked on the smallest absolute element of
     * the result and the remaining elements may still have a larger error.
     */
    private static final double EPS_EXPONENTIAL = 1e-7;

    private static final double ABSOLUTE_ERROR_IMPROPER_1 = 1e-5;

    @Test
    void integrate_whenExponential2_returnsExpectedResult() throws WrongSizeException, IntegrationException {
        // for single dimension functions, integral of f(t) = e^(a*t) is f(t) = 1/a*e^(a*t)
        // Consequently, for matrix function f(t) = e^(A*t) where A is a matrix, should be:
        // f(t) = A^-1 * e^(A*t)
        final var randomizer = new UniformRandomizer(new Random(SEED));

        // for matrix A = [0    1]
        //                [-1   0]

        // f(t) = e^(A*t) = [cos(t)     sin(t)]
        //                  [-sin(t)    cos(t)]

        // matrix A is antisymmetric and orthonormal, consequently, its inverse is its transpose:
        // A^-1 = A^T = [0  -1]
        //              [1   0]

        // A * A^T = [0   1][0  -1] = [0*0 + 1*1    -1*0 + 1*0 ] = [1  0]
        //           [-1  0][1   0]   [-1*0 + 0*1   -1*-1 + 0*0]   [0  1]

        // Consequently, the integral of f(t) is:
        // F(t) = A^-1*e^(A*t) = [0  -1][cos(t)     sin(t)] = [sin(t)  -cos(t)]
        //                       [1   0][-sin(t)    cos(t)]   [cos(t)   sin(t)]

        // Finally, integral of f(t) between a and b will be:
        // F(b) - F(a) = [sin(b) - sin(a)       -cos(b) + cos(a)]
        //               [cos(b) - cos(a)        sin(b) - sin(a)]

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double t, Matrix result) {
                // f(t) = e^(A*t) = [cost(t)  sin(t)]
                //                  [-sin(t)  cos(t)]
                final var c = Math.cos(t);
                final var s = Math.sin(t);
                result.setElementAtIndex(0, c);
                result.setElementAtIndex(1, -s);
                result.setElementAtIndex(2, s);
                result.setElementAtIndex(3, c);
            }

            @Override
            public int getRows() {
                return 2;
            }

            @Override
            public int getColumns() {
                return 2;
            }
        };

        // F(b) - F(a) = [sin(b) - sin(a)       -cos(b) + cos(a)]
        //               [cos(b) - cos(a)        sin(b) - sin(a)]
        for (var t = 0; t < TIMES; t++) {
            final var interval = drawInterval(randomizer);
            final var a = interval[0];
            final var b = interval[1];

            final var expectedResult = new Matrix(2, 2);
            final var sinDiff = Math.sin(b) - Math.sin(a);
            final var cosDiff = Math.cos(b) - Math.cos(a);
            expectedResult.setElementAtIndex(0, sinDiff);
            expectedResult.setElementAtIndex(1, cosDiff);
            expectedResult.setElementAtIndex(2, -cosDiff);
            expectedResult.setElementAtIndex(3, sinDiff);

            final var integrator = new SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator(a, b, listener,
                    EPS_EXPONENTIAL);

            final var integrationResult = new Matrix(2, 2);
            integrator.integrate(integrationResult);

            // check
            assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_EXPONENTIAL));
        }
    }

    @Test
    void integrate_whenImproperIntegrandWithSingularities_returnsExpectedResult() throws IntegrationException,
            WrongSizeException {

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, Math.log(point) * Math.log(1 - point));
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expected = 2.0 - Math.PI * Math.PI / 6.0;
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_IMPROPER_1));
    }

    @Test
    void getIntegratorType_returnsExpectedValue() throws WrongSizeException {
        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                // no action needed
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() throws WrongSizeException {
        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                // no action needed
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());
    }

    /**
     * Draws an integration interval [a, b], at least {@link #MIN_INTERVAL} long, where neither sin(b) - sin(a)
     * nor cos(b) - cos(a) is close to zero.
     *
     * @param randomizer randomizer.
     * @return array containing a and b.
     */
    private static double[] drawInterval(final UniformRandomizer randomizer) {
        double a;
        double b;
        do {
            a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE - MIN_INTERVAL);
            b = randomizer.nextDouble(a + MIN_INTERVAL, MAX_VALUE);
        } while (Math.abs(Math.sin(b) - Math.sin(a)) < MIN_ABS_DIFF
                || Math.abs(Math.cos(b) - Math.cos(a)) < MIN_ABS_DIFF);
        return new double[]{a, b};
    }

    /**
     * Maximum number of function evaluations expected when integrating a smooth function. Mid-point quadratures
     * triple the number of points on each refinement, so a 9/8 extrapolation with a reduced maximum number of
     * refinements must converge with fewer evaluations than the former 4/3 extrapolation did.
     * Note: 729 (3^6) evaluations is the minimum possible, since convergence is first checked
     * at the 7th stage (after JMIN); smooth integrands converge there. Before the fix this took
     * 6561 to 59049 evaluations.
     */
    private static final int MAX_EVALUATIONS = 1000;

    @Test
    void integrate_whenSmoothFunction_returnsExpectedResultWithFewEvaluations()
            throws IntegrationException, WrongSizeException {
        final var expected0 = Math.exp(2.0) - Math.exp(1.0);
        final var expected1 = 0.5 * (Math.exp(4.0) - Math.exp(2.0));
        final var counter = new int[1];

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(final double point, final Matrix result) {
                counter[0]++;
                result.setElementAtIndex(0, Math.exp(point));
                result.setElementAtIndex(1, Math.exp(2.0 * point));
            }

            @Override
            public int getRows() {
                return 2;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator(1.0, 2.0, listener, 1e-8);

        final var integrationResult = new Matrix(2, 1);
        integrator.integrate(integrationResult);

        assertEquals(expected0, integrationResult.getElementAtIndex(0),
                1e-6 * Math.max(1.0, Math.abs(expected0)));
        assertEquals(expected1, integrationResult.getElementAtIndex(1),
                1e-6 * Math.max(1.0, Math.abs(expected1)));
        assertTrue(counter[0] <= MAX_EVALUATIONS);
    }
}
