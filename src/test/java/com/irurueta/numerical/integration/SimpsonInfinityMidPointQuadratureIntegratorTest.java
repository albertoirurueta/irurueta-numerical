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

import com.irurueta.statistics.NormalDist;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;

import org.junit.jupiter.api.Test;

class SimpsonInfinityMidPointQuadratureIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double ABSOLUTE_ERROR_GAUSSIAN = 1e-10;

    /**
     * Fixed seed used to make the randomized Gaussian integration test deterministic.
     */
    private static final long SEED = 2_718_281_828L;

    /**
     * Number of randomized trials executed by the Gaussian integration test.
     */
    private static final int TIMES = 20;

    /**
     * Minimum length of the integration interval, and minimum value of its start point, used in the Gaussian
     * integration test.
     */
    private static final double MIN_INTERVAL = 0.5;

    /**
     * Minimum expected integral in the Gaussian integration test. The integrator uses a relative convergence
     * criterion, which cannot be met when the expected integral is close to zero.
     */
    private static final double MIN_EXPECTED = 1e-3;

    @Test
    void integrate_whenGaussian_returnsExpectedResult() throws IntegrationException {
        final var randomizer = new UniformRandomizer(new Random(SEED));
        for (var t = 0; t < TIMES; t++) {
            assertGaussianIntegration(randomizer);
        }
    }

    @Test
    void getIntegratorType_returnsExpectedValue() {
        final var integrator = new SimpsonInfinityMidPointQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() {
        final var integrator = new SimpsonInfinityMidPointQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(QuadratureType.INFINITY_MID_POINT, integrator.getQuadratureType());
    }

    private static void assertGaussianIntegration(final UniformRandomizer randomizer) throws IntegrationException {
        // Infinity Mid-Point requires that a * b > 0.0
        // (either a and b are positive, or both a and b are negative)
        double a;
        double b;
        double mu;
        double sigma;
        double expected;
        do {
            a = randomizer.nextDouble(MIN_INTERVAL, MAX_VALUE - MIN_INTERVAL);
            b = randomizer.nextDouble(a + MIN_INTERVAL, MAX_VALUE);
            mu = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            sigma = ABSOLUTE_ERROR_GAUSSIAN + Math.abs(randomizer.nextDouble(a, MAX_VALUE));
            expected = NormalDist.cdf(b, mu, sigma) - NormalDist.cdf(a, mu, sigma);
        } while (expected < MIN_EXPECTED);

        final var finalMu = mu;
        final var finalSigma = sigma;
        final var integrator = new SimpsonInfinityMidPointQuadratureIntegrator(a, b,
                point -> NormalDist.p(point, finalMu, finalSigma));
        final var result = integrator.integrate();

        assertEquals(expected, result, ABSOLUTE_ERROR_GAUSSIAN);
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
        final var expected = Math.exp(-1.0) - Math.exp(-2.0);
        final var counter = new int[1];

        final var integrator = new SimpsonInfinityMidPointQuadratureIntegrator(1.0, 2.0, point -> {
            counter[0]++;
            return Math.exp(-point);
        }, 1e-8);
        final var result = integrator.integrate();

        assertEquals(expected, result, 1e-6 * Math.max(1.0, Math.abs(expected)));
        assertTrue(counter[0] <= MAX_EVALUATIONS);
    }
}
