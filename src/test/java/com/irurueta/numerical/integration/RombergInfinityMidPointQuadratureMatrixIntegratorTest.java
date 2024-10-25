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
import static org.junit.jupiter.api.Assertions.assertThrows;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;

import org.junit.jupiter.api.Test;

class RombergInfinityMidPointQuadratureMatrixIntegratorTest {

    private static final double ALMOST_INFINITY = 1e99;

    @Test
    void integrate_whenImproperIntegrandWithSingularities_returnsExpectedResult() throws WrongSizeException {

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(final double point, final Matrix result) {
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

        final var integrator = new RombergInfinityMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);

        final var integrationResult = new Matrix(1, 1);
        assertThrows(IntegrationException.class, () -> integrator.integrate(integrationResult));
    }

    @Test
    void integrate_whenImproperIntegralFromZeroToInfinity3_returnsWrongResult() throws WrongSizeException {

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(final double point, final Matrix result) {
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

        final var integrator = new RombergInfinityMidPointQuadratureMatrixIntegrator(0.0, ALMOST_INFINITY, listener);

        final var integrationResult = new Matrix(1, 1);
        assertThrows(IntegrationException.class, () -> integrator.integrate(integrationResult));
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

        final var integrator = new RombergInfinityMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
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

        final var integrator = new RombergInfinityMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());
    }
}