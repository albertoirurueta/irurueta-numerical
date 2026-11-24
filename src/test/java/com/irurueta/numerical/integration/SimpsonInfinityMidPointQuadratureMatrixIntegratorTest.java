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
import com.irurueta.statistics.NormalDist;

import org.junit.jupiter.api.Test;

class SimpsonInfinityMidPointQuadratureMatrixIntegratorTest {

    private static final double ABSOLUTE_ERROR_GAUSSIAN = 1e-3;

    @Test
    void integrate_whenGaussian_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        final var a = 10.0;
        final var b = 1.0;
        final var mu = 5.0;
        final var sigma = 0.1;

        final var expected = NormalDist.cdf(b, mu, sigma) - NormalDist.cdf(a, mu, sigma);

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, NormalDist.p(point, mu, sigma));
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

        final var integrator = new SimpsonInfinityMidPointQuadratureMatrixIntegrator(a, b, listener);

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_GAUSSIAN));
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

        final var integrator = new SimpsonInfinityMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
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

        final var integrator = new SimpsonInfinityMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());
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
        final var expected0 = Math.exp(-1.0) - Math.exp(-2.0);
        final var expected1 = 0.5 * (Math.exp(-2.0) - Math.exp(-4.0));
        final var counter = new int[1];

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(final double point, final Matrix result) {
                counter[0]++;
                result.setElementAtIndex(0, Math.exp(-point));
                result.setElementAtIndex(1, Math.exp(-2.0 * point));
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

        final var integrator = new SimpsonInfinityMidPointQuadratureMatrixIntegrator(1.0, 2.0, listener, 1e-8);

        final var integrationResult = new Matrix(2, 1);
        integrator.integrate(integrationResult);

        assertEquals(expected0, integrationResult.getElementAtIndex(0),
                1e-6 * Math.max(1.0, Math.abs(expected0)));
        assertEquals(expected1, integrationResult.getElementAtIndex(1),
                1e-6 * Math.max(1.0, Math.abs(expected1)));
        assertTrue(counter[0] <= MAX_EVALUATIONS);
    }
}
