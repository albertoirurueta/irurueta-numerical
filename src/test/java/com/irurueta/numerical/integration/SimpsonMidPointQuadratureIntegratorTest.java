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
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.Gamma;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;

import org.junit.jupiter.api.Test;

class SimpsonMidPointQuadratureIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double MAX_LAMBDA = 1.0;

    /**
     * Fixed seed used to make randomized exponential integration tests deterministic.
     */
    private static final long SEED = 7_736_412_987L;

    /**
     * Number of randomized trials executed by exponential integration tests.
     */
    private static final int TIMES = 20;

    /**
     * Minimum length of the integration interval used in exponential integration tests.
     */
    private static final double MIN_INTERVAL = 0.5;

    /**
     * Minimum absolute value of lambda used in exponential integration tests, to avoid ill-conditioned values
     * close to zero.
     */
    private static final double MIN_ABS_LAMBDA = 0.1;

    /**
     * Relative convergence accuracy used by the integrator in exponential integration tests.
     */
    private static final double EPS_EXPONENTIAL = 1e-8;

    /**
     * Minimum absolute value of lambda used in the small lambda (near-constant integrand) exponential test.
     */
    private static final double MIN_SMALL_LAMBDA = 1e-6;

    /**
     * Minimum distance between the root of the first degree polynomial and the middle of the integration interval,
     * so that the expected integral is not close to zero (the relative convergence criterion of the integrator
     * cannot be met for integrals close to zero).
     */
    private static final double MIN_ROOT_DISTANCE = 1.0;

    private static final double ABSOLUTE_ERROR_1 = 1e-10;

    private static final double ABSOLUTE_ERROR_EXPONENTIAL = 1e-6;

    private static final double ABSOLUTE_ERROR_IMPROPER_1 = 1e-8;

    private static final double ABSOLUTE_ERROR_IMPROPER_3 = 1e-5;

    private static final double ALMOST_INFINITY = 1e99;

    @Test
    void integrate_whenFirstDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration();
    }

    @Test
    void integrate_whenExponential_returnsExpectedResult() throws IntegrationException {
        final var randomizer = new UniformRandomizer(new Random(SEED));
        for (var t = 0; t < TIMES; t++) {
            final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE - MIN_INTERVAL);
            final var b = randomizer.nextDouble(a + MIN_INTERVAL, MAX_VALUE);
            final var lambda = randomizer.nextDouble(MIN_ABS_LAMBDA, MAX_LAMBDA)
                    * (randomizer.nextBoolean() ? 1.0 : -1.0);

            final var expected = 1.0 / lambda * (Math.exp(lambda * b) - Math.exp(lambda * a));

            final var integrator = new SimpsonMidPointQuadratureIntegrator(a, b, point -> Math.exp(lambda * point),
                    EPS_EXPONENTIAL);
            final var result = integrator.integrate();

            // integrator convergence is relative, so scale tolerance with result magnitude
            final var tolerance = ABSOLUTE_ERROR_EXPONENTIAL * Math.max(1.0, Math.abs(expected));
            assertEquals(expected, result, tolerance);
        }
    }

    @Test
    void integrate_whenExponentialWithSmallLambda_returnsExpectedResult() throws IntegrationException {
        final var randomizer = new UniformRandomizer(new Random(SEED));
        for (var t = 0; t < TIMES; t++) {
            assertSmallLambdaExponentialIntegration(randomizer);
        }
    }

    @Test
    void integrate_whenImproperIntegrandWithSingularities_returnsExpectedResult() throws IntegrationException {
        final var expected = 2.0 - Math.PI * Math.PI / 6.0;

        final var integrator = new SimpsonMidPointQuadratureIntegrator(0.0, 1.0,
                point -> Math.log(point) * Math.log(1 - point));
        final var result = integrator.integrate();

        assertEquals(expected, result, ABSOLUTE_ERROR_IMPROPER_1);
    }

    @Test
    void integrate_whenImproperIntegralFromZeroToInfinity3_returnsWrongResult() throws IntegrationException {
        final var expected = 0.5 * Math.exp(Gamma.gammln(5.0 / 14.0));
        assertEquals(1.24663, expected, ABSOLUTE_ERROR_IMPROPER_3);

        final var integrator = new SimpsonMidPointQuadratureIntegrator(0.0, ALMOST_INFINITY,
                point -> Math.pow(point, -2.0 / 7.0) * Math.exp(-point * point));
        final var result = integrator.integrate();

        assertNotEquals(expected, result, ABSOLUTE_ERROR_IMPROPER_3);
    }

    @Test
    void getIntegratorType_returnsExpectedValue() {
        final var integrator = new SimpsonMidPointQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(IntegratorType.SIMPSON, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() {
        final var integrator = new SimpsonMidPointQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(QuadratureType.MID_POINT, integrator.getQuadratureType());
    }

    private static void assertPolynomialIntegration() throws IntegrationException {
        final var randomizer = new UniformRandomizer(new Random(SEED));
        for (var t = 0; t < TIMES; t++) {
            assertPolynomialIntegration(randomizer);
        }
    }

    private static void assertPolynomialIntegration(final UniformRandomizer randomizer) throws IntegrationException {
        // set integration interval
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE - MIN_INTERVAL);
        final var b = randomizer.nextDouble(a + MIN_INTERVAL, MAX_VALUE);

        final var polynomial = buildPolynomial(randomizer, a, b);
        final var integrationPolynomial = polynomial.integrationAndReturnNew();

        final var expected = integrationPolynomial.evaluate(b) - integrationPolynomial.evaluate(a);

        final var integrator = new SimpsonMidPointQuadratureIntegrator(a, b, polynomial::evaluate);
        final var result = integrator.integrate();

        assertEquals(expected, result, SimpsonMidPointQuadratureIntegratorTest.ABSOLUTE_ERROR_1);
    }

    private static void assertSmallLambdaExponentialIntegration(final UniformRandomizer randomizer)
            throws IntegrationException {
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE - MIN_INTERVAL);
        final var b = randomizer.nextDouble(a + MIN_INTERVAL, MAX_VALUE);
        final var lambda = randomizer.nextDouble(MIN_SMALL_LAMBDA, MIN_ABS_LAMBDA)
                * (randomizer.nextBoolean() ? 1.0 : -1.0);

        // expm1 avoids the cancellation of exp(lambda * b) - exp(lambda * a) when lambda is small
        final var expected = Math.exp(lambda * a) * Math.expm1(lambda * (b - a)) / lambda;

        final var integrator = new SimpsonMidPointQuadratureIntegrator(a, b, point -> Math.exp(lambda * point),
                EPS_EXPONENTIAL);
        final var result = integrator.integrate();

        // integrator convergence is relative, so scale tolerance with result magnitude
        final var tolerance = ABSOLUTE_ERROR_EXPONENTIAL * Math.max(1.0, Math.abs(expected));
        assertEquals(expected, result, tolerance);
    }

    /**
     * Builds a first degree polynomial whose root is at least {@link #MIN_ROOT_DISTANCE} away from the middle of
     * the integration interval [a, b], so that its integral over that interval is not close to zero.
     *
     * @param randomizer randomizer.
     * @param a          start of integration interval.
     * @param b          end of integration interval.
     * @return first degree polynomial.
     */
    private static Polynomial buildPolynomial(final UniformRandomizer randomizer, final double a, final double b) {
        final var middle = 0.5 * (a + b);
        final var distance = randomizer.nextDouble(MIN_ROOT_DISTANCE, MAX_VALUE);
        final var root = randomizer.nextBoolean() ? middle + distance : middle - distance;
        return new Polynomial(-root, 1.0);
    }

    /**
     * Maximum number of function evaluations expected when integrating a smooth function. Mid-point quadratures
     * triple the number of points on each refinement, so a 9/8 extrapolation with a reduced maximum number of
     * refinements must converge with fewer evaluations than the former 4/3 extrapolation did.
     */
    private static final int MAX_EVALUATIONS = 1000;

    @Test
    void integrate_whenSmoothFunction_returnsExpectedResultWithFewEvaluations() throws IntegrationException {
        final var expected = Math.exp(2.0) - Math.exp(1.0);
        final var counter = new int[1];

        final var integrator = new SimpsonMidPointQuadratureIntegrator(1.0, 2.0, point -> {
            counter[0]++;
            return Math.exp(point);
        }, 1e-8);
        final var result = integrator.integrate();

        assertEquals(expected, result, 1e-6 * Math.max(1.0, Math.abs(expected)));
        assertTrue(counter[0] <= MAX_EVALUATIONS);
    }
}
