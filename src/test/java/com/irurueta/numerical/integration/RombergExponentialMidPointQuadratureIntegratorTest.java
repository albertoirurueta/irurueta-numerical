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

import com.irurueta.statistics.Gamma;

import org.junit.jupiter.api.Test;

class RombergExponentialMidPointQuadratureIntegratorTest {

    private static final double ABSOLUTE_ERROR_IMPROPER_3 = 1e-4;

    @Test
    void integrate_whenImproperIntegralFromZeroToInfinity3_returnsExpectedResult() throws IntegrationException {
        final var expected = 0.5 * Math.exp(Gamma.gammln(5.0 / 14.0));
        assertEquals(1.24663, expected, ABSOLUTE_ERROR_IMPROPER_3);

        final var integrator = new RombergExponentialMidPointQuadratureIntegrator(0.0,
                point -> Math.pow(point, -2.0 / 7.0) * Math.exp(-point * point));
        final var result = integrator.integrate();

        assertEquals(expected, result, ABSOLUTE_ERROR_IMPROPER_3);
    }

    @Test
    void getIntegratorType_returnsExpectedValue() {
        final var integrator = new RombergExponentialMidPointQuadratureIntegrator(0.0, null);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() {
        final var integrator = new RombergExponentialMidPointQuadratureIntegrator(0.0, null);
        assertEquals(QuadratureType.EXPONENTIAL_MID_POINT, integrator.getQuadratureType());
    }
}