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

import org.junit.jupiter.api.Test;

class SimpsonLowerSquareRootMidPointQuadratureIntegratorTest {

    private static final double ABSOLUTE_ERROR_IMPROPER_1 = 1e-5;

    @Test
    void integrate_whenImproperIntegrandWithSingularities_returnsExpectedResult() throws IntegrationException {
        final var expected = 2.0 - Math.PI * Math.PI / 6.0;

        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureIntegrator(0.0, 1.0,
                point -> Math.log(point) * Math.log(1 - point));
        final var result = integrator.integrate();

        assertEquals(expected, result, ABSOLUTE_ERROR_IMPROPER_1);
    }

    @Test
    void getIntegratorType_returnsExpectedValue() {
        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() {
        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());
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
    void integrate_whenSmoothFunction_returnsExpectedResultWithFewEvaluations() throws IntegrationException {
        final var expected = Math.exp(2.0) - Math.exp(1.0);
        final var counter = new int[1];

        final var integrator = new SimpsonLowerSquareRootMidPointQuadratureIntegrator(1.0, 2.0, point -> {
            counter[0]++;
            return Math.exp(point);
        }, 1e-8);
        final var result = integrator.integrate();

        assertEquals(expected, result, 1e-6 * Math.max(1.0, Math.abs(expected)));
        assertTrue(counter[0] <= MAX_EVALUATIONS);
    }
}
