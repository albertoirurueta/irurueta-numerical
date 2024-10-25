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

import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.jupiter.api.Assertions.*;

class IntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double MIN_LAMBDA = -1.0;

    private static final double MAX_LAMBDA = 1.0;

    private static final double EPS = 1e-6;

    @Test
    void create_whenAccuracyIntegratorAndQuadratureTypes_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE, QuadratureType.MID_POINT);
        assertInstanceOf(MidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(InfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(LowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(UpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(DoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> Integrator.create(a, b, listener, EPS,
                IntegratorType.QUADRATURE, QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(SimpsonTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON, QuadratureType.MID_POINT);
        assertInstanceOf(SimpsonMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(SimpsonInfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonLowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonUpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(SimpsonDoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> Integrator.create(a, b, listener, EPS,
                IntegratorType.SIMPSON, QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(RombergTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG, QuadratureType.MID_POINT);
        assertInstanceOf(RombergMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(RombergInfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergLowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergUpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.EXPONENTIAL_MID_POINT);
        assertInstanceOf(RombergExponentialMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.EXPONENTIAL_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(RombergDoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());
    }

    @Test
    void create_whenDefaultAccuracyIntegratorAndQuadratureTypes_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE, QuadratureType.MID_POINT);
        assertInstanceOf(MidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(InfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(LowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(UpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(DoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> Integrator.create(a, b, listener, IntegratorType.QUADRATURE,
                QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(SimpsonTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.MID_POINT);
        assertInstanceOf(SimpsonMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(SimpsonInfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonLowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(SimpsonUpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON, QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(SimpsonDoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> Integrator.create(a, b, listener, IntegratorType.SIMPSON,
                QuadratureType.EXPONENTIAL_MID_POINT));

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(RombergTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.MID_POINT);
        assertInstanceOf(RombergMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(RombergInfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG,
                QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergLowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG,
                QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(RombergUpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.EXPONENTIAL_MID_POINT);
        assertInstanceOf(RombergExponentialMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.EXPONENTIAL_MID_POINT, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG, QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(RombergDoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());
    }

    @Test
    void create_whenAccuracyAndIntegratorType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = Integrator.create(a, b, listener, EPS, IntegratorType.QUADRATURE);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.SIMPSON);
        assertInstanceOf(SimpsonTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, EPS, IntegratorType.ROMBERG);
        assertInstanceOf(RombergTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void create_whenDefaultAccuracyAndIntegratorType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = Integrator.create(a, b, listener, IntegratorType.QUADRATURE);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.SIMPSON);
        assertInstanceOf(SimpsonTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = Integrator.create(a, b, listener, IntegratorType.ROMBERG);
        assertInstanceOf(RombergTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void create_whenAccuracyAndDefaultIntegratorType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        final var integrator = Integrator.create(a, b, listener, EPS);
        assertInstanceOf(RombergTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void create_whenDefaultAccuracyAndDefaultIntegratorType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        final var integrator = Integrator.create(a, b, listener);
        assertInstanceOf(RombergTrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    @Test
    void comparePerformance() throws IntegrationException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var lambda = randomizer.nextDouble(MIN_LAMBDA, MAX_LAMBDA);

        final var evaluations = new int[1];
        final SingleDimensionFunctionEvaluatorListener listener =
                point -> {
                    evaluations[0] += 1;
                    return Math.exp(lambda * point);
                };

        final var expected = 1.0 / lambda * (Math.exp(lambda * b) - Math.exp(lambda * a));

        int bestEvaluations = Integer.MAX_VALUE;
        IntegratorType bestEvaluationsIntegratorType = null;
        QuadratureType bestEvaluationsQuadratureType = null;

        long bestDuration = Long.MAX_VALUE;
        IntegratorType bestDurationIntegratorType = null;
        QuadratureType bestDurationQuadratureType = null;

        double bestError = Double.MAX_VALUE;
        IntegratorType bestErrorIntegratorType = null;
        QuadratureType bestErrorQuadratureType = null;

        for (final var integratorType : IntegratorType.values()) {
            for (final var quadratureType : QuadratureType.values()) {
                final var integrator = Integrator.create(a, b, listener);
                assertNotNull(integrator);

                evaluations[0] = 0;
                final var start = System.nanoTime();
                final var result = integrator.integrate();
                final var end = System.nanoTime();
                final var duration = end - start;
                final var error = Math.abs(expected - result);

                Logger.getGlobal().log(Level.INFO, String.format(
                        "Integrator type: %s, Quadrature type: %s, evaluations: %s, duration: %d" + "ns, error: %f",
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
                "Best duration - Integrator type: %s, Quadrature type: %s, duration: %d" + "ns",
                bestDurationIntegratorType, bestDurationQuadratureType, bestDuration));

        Logger.getGlobal().log(Level.INFO, String.format(
                "Best error - Integrator type: %s, Quadrature type: %s, error: %f",
                bestErrorIntegratorType, bestErrorQuadratureType, bestError));
    }
}