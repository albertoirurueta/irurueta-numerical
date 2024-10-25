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

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.jupiter.api.Assertions.*;

class MatrixIntegratorTest {
    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double MIN_LAMBDA = -1.0;

    private static final double MAX_LAMBDA = 1.0;

    private static final double EPS = 1e-6;

    @Test
    void create_whenAccuracyIntegratorAndQuadratureTypes_returnsExpectedIntegrator() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

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

        var integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(TrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE, QuadratureType.MID_POINT);
        assertInstanceOf(MidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(InfinityMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(LowerSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(UpperSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(DoubleExponentialRuleQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> MatrixIntegrator.create(a, b, listener, EPS,
                IntegratorType.QUADRATURE, QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(SimpsonTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON, QuadratureType.MID_POINT);
        assertInstanceOf(SimpsonMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(SimpsonInfinityMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonUpperSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(SimpsonDoubleExponentialRuleQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> MatrixIntegrator.create(a, b, listener, EPS,
                IntegratorType.SIMPSON, QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(RombergTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG, QuadratureType.MID_POINT);
        assertInstanceOf(RombergMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(RombergInfinityMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergLowerSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergUpperSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.EXPONENTIAL_MID_POINT);
        assertInstanceOf(RombergExponentialMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.EXPONENTIAL_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(RombergDoubleExponentialRuleQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());
    }

    @Test
    void create_whenDefaultAccuracyIntegratorAndQuadratureTypes_returnsExpectedIntegrator() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

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

        var integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(TrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE, QuadratureType.MID_POINT);
        assertInstanceOf(MidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(InfinityMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(LowerSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(UpperSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(DoubleExponentialRuleQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> MatrixIntegrator.create(a, b, listener,
                IntegratorType.QUADRATURE, QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(SimpsonTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.MID_POINT);
        assertInstanceOf(SimpsonMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(SimpsonInfinityMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonLowerSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonUpperSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(SimpsonDoubleExponentialRuleQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> MatrixIntegrator.create(a, b, listener,
                IntegratorType.SIMPSON, QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(RombergTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.MID_POINT);
        assertInstanceOf(RombergMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(RombergInfinityMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergLowerSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergUpperSquareRootMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG,
                QuadratureType.EXPONENTIAL_MID_POINT);
        assertInstanceOf(RombergExponentialMidPointQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.EXPONENTIAL_MID_POINT, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(RombergDoubleExponentialRuleQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());
    }

    @Test
    void create_whenAccuracyAndIntegratorType_returnsExpectedIntegrator() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

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

        var integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE);
        assertInstanceOf(TrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.SIMPSON);
        assertInstanceOf(SimpsonTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, EPS, IntegratorType.ROMBERG);
        assertInstanceOf(RombergTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void create_whenDefaultAccuracyAndIntegratorType_returnsExpectedIntegrator() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

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

        var integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.QUADRATURE);
        assertInstanceOf(TrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.SIMPSON);
        assertInstanceOf(SimpsonTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = MatrixIntegrator.create(a, b, listener, IntegratorType.ROMBERG);
        assertInstanceOf(RombergTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void create_whenAccuracyAndDefaultIntegratorType_returnsExpectedIntegrator() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

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

        final var integrator = MatrixIntegrator.create(a, b, listener, EPS);
        assertInstanceOf(RombergTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void create_whenDefaultAccuracyAndDefaultIntegratorType_returnsExpectedIntegrator() throws WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

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

        final var integrator = MatrixIntegrator.create(a, b, listener);
        assertInstanceOf(RombergTrapezoidalQuadratureMatrixIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void comparePerformance() throws IntegrationException, WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var lambda = randomizer.nextDouble(MIN_LAMBDA, MAX_LAMBDA);

        final var evaluations = new int[1];
        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                evaluations[0] += 1;
                result.setElementAtIndex(0, Math.exp(lambda * point));
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

        final var expected = 1.0 / lambda * (Math.exp(lambda * b) - Math.exp(lambda * a));

        var bestEvaluations = Integer.MAX_VALUE;
        IntegratorType bestEvaluationsIntegratorType = null;
        QuadratureType bestEvaluationsQuadratureType = null;

        var bestDuration = Long.MAX_VALUE;
        IntegratorType bestDurationIntegratorType = null;
        QuadratureType bestDurationQuadratureType = null;

        var bestError = Double.MAX_VALUE;
        IntegratorType bestErrorIntegratorType = null;
        QuadratureType bestErrorQuadratureType = null;

        for (final var integratorType : IntegratorType.values()) {
            for (final var quadratureType : QuadratureType.values()) {
                final var integrator = MatrixIntegrator.create(a, b, listener);
                assertNotNull(integrator);

                evaluations[0] = 0;
                final var start = System.nanoTime();
                final var result = new Matrix(1, 1);
                integrator.integrate(result);
                final var end = System.nanoTime();
                final var duration = end - start;
                final var error = Math.abs(expected - result.getElementAtIndex(0));

                Logger.getGlobal().log(Level.INFO, String.format(
                        "Integrator type: %s, Quadrature type: %s, evaluations: %d, duration: %d ns, error: %f",
                        integratorType, quadratureType, evaluations[0], duration, error));

                if (evaluations[0] < bestEvaluations) {
                    bestEvaluations = evaluations[0];
                    bestEvaluationsIntegratorType = integratorType;
                    bestEvaluationsQuadratureType = quadratureType;
                }

                if (duration < bestDuration) {
                    bestDuration = duration;
                    bestDurationIntegratorType = integratorType;
                    bestDurationQuadratureType = quadratureType;
                }

                if (error < bestError) {
                    bestError = error;
                    bestErrorIntegratorType = integratorType;
                    bestErrorQuadratureType = quadratureType;
                }
            }
        }

        Logger.getGlobal().log(Level.INFO, String.format(
                "Best evaluations - Integrator type: %s, Quadrature type: %s, evaluations: %d",
                bestEvaluationsIntegratorType, bestEvaluationsQuadratureType, bestEvaluations));

        Logger.getGlobal().log(Level.INFO, String.format(
                "Best duration - Integrator type: %s, Quadrature type: %s, duration: %d ns",
                bestDurationIntegratorType, bestDurationQuadratureType, bestDuration));

        Logger.getGlobal().log(Level.INFO, String.format(
                "Best error - Integrator type: %s, Quadrature type: %s, error: %f",
                bestErrorIntegratorType, bestErrorQuadratureType, bestError));
    }
}