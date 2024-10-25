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
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.Gamma;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class MidPointQuadratureMatrixIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double MIN_LAMBDA = -1.0;

    private static final double MAX_LAMBDA = 1.0;

    private static final double ABSOLUTE_ERROR_1 = 1e-10;

    private static final double ABSOLUTE_ERROR_EXPONENTIAL = 1e-6;

    private static final double ABSOLUTE_ERROR_IMPROPER_1 = 1e-8;

    private static final double ABSOLUTE_ERROR_IMPROPER_3 = 1e-5;

    private static final double ALMOST_INFINITY = 1e99;

    @Test
    void integrate_whenFirstDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration();
    }

    @Test
    void integrate_whenExponential1_returnsExpectedResult() throws WrongSizeException, IntegrationException {

        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var lambda = randomizer.nextDouble(MIN_LAMBDA, MAX_LAMBDA);

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double t, Matrix result) {
                result.setElementAtIndex(0, Math.exp(lambda * t));
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

        final var integrator = new MidPointQuadratureMatrixIntegrator(a, b, listener);

        final var expected = 1.0 / lambda * (Math.exp(lambda * b) - Math.exp(lambda * a));

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_EXPONENTIAL));
    }

    @Test
    void integrate_whenImproperIntegrandWithSingularities_returnsExpectedResult() throws IntegrationException,
            WrongSizeException {

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, Math.log(point) * Math.log(1.0 - point));
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

        final var integrator = new MidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);

        final var expected = 2.0 - Math.PI * Math.PI / 6.0;

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_IMPROPER_1));
    }

    @Test
    void integrate_whenImproperIntegralFromZeroToInfinity3_returnsWrongResult() throws IntegrationException,
            WrongSizeException {

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, Math.pow(point, -2.0 / 7.0) * Math.exp(-point * point));
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

        final var integrator = new MidPointQuadratureMatrixIntegrator(0.0, ALMOST_INFINITY, listener);

        final var expected = 0.5 * Math.exp(Gamma.gammln(5.0 / 14.0));
        assertEquals(1.24663, expected, ABSOLUTE_ERROR_IMPROPER_3);

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertFalse(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_IMPROPER_3));
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

        final var integrator = new MidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
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

        final var integrator = new MidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());
    }

    private static void assertPolynomialIntegration() throws IntegrationException, WrongSizeException {
        final var polynomial = buildPolynomial();
        final var integrationPolynomial = polynomial.integrationAndReturnNew();

        // set integration interval
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final var expected = integrationPolynomial.evaluate(b) - integrationPolynomial.evaluate(a);

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, polynomial.evaluate(point));
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

        final var integrator = new MidPointQuadratureMatrixIntegrator(a, b, listener);

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, MidPointQuadratureMatrixIntegratorTest.ABSOLUTE_ERROR_1));
    }

    private static Polynomial buildPolynomial() {
        final var randomizer = new UniformRandomizer();
        final var root = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        return new Polynomial(-root, 1.0);
    }
}