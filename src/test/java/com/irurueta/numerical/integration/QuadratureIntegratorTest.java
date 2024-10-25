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

import static org.junit.jupiter.api.Assertions.*;

class QuadratureIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double EPS = 1e-6;

    @Test
    void create_whenAccuracyAndQuadratureType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = QuadratureIntegrator.create(a, b, listener, EPS,
                QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, EPS, QuadratureType.MID_POINT);
        assertInstanceOf(MidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, EPS, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(InfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, EPS, QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(LowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, EPS, QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(UpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, EPS, QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(DoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class,
                () -> QuadratureIntegrator.create(a, b, listener, EPS, QuadratureType.EXPONENTIAL_MID_POINT));
    }

    @Test
    void create_whenDefaultAccuracyAndQuadratureType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = QuadratureIntegrator.create(a, b, listener, QuadratureType.TRAPEZOIDAL);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, QuadratureType.MID_POINT);
        assertInstanceOf(MidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, QuadratureType.INFINITY_MID_POINT);
        assertInstanceOf(InfinityMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, QuadratureType.LOWER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(LowerSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, QuadratureType.UPPER_SQUARE_ROOT_MID_POINT);
        assertInstanceOf(UpperSquareRootMidPointQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.UPPER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener, QuadratureType.DOUBLE_EXPONENTIAL_RULE);
        assertInstanceOf(DoubleExponentialRuleQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.DOUBLE_EXPONENTIAL_RULE, integrator.getQuadratureType());

        assertThrows(IllegalArgumentException.class, () -> QuadratureIntegrator.create(a, b, listener,
                QuadratureType.EXPONENTIAL_MID_POINT));
    }

    @Test
    void create_whenNoQuadratureType_returnsExpectedIntegrator() {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final SingleDimensionFunctionEvaluatorListener listener = point -> 0.0;

        var integrator = QuadratureIntegrator.create(a, b, listener, EPS);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());

        integrator = QuadratureIntegrator.create(a, b, listener);
        assertInstanceOf(TrapezoidalQuadratureIntegrator.class, integrator);
        assertEquals(IntegratorType.QUADRATURE, integrator.getIntegratorType());
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }
}