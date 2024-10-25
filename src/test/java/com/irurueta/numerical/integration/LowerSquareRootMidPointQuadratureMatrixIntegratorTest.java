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

import org.junit.jupiter.api.Test;

class LowerSquareRootMidPointQuadratureMatrixIntegratorTest {

    private static final double ABSOLUTE_ERROR_IMPROPER_1 = 1e-5;

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

        final var integrator = new LowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);

        final var expected = 2.0 - Math.PI * Math.PI / 6.0;

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
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

        final var integrator = new LowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
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

        final var integrator = new LowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());
    }
}