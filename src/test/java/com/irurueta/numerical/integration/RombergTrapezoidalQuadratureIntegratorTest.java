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

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.NormalDist;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class RombergTrapezoidalQuadratureIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double MIN_LAMBDA = -1.0;

    private static final double MAX_LAMBDA = 1.0;

    private static final double ABSOLUTE_ERROR_1 = 1e-9;

    private static final double ABSOLUTE_ERROR_6 = 1e-6;

    private static final double ABSOLUTE_ERROR_GAUSSIAN = 1e-8;

    private static final double ABSOLUTE_ERROR_EXPONENTIAL = 1e-7;

    private static final double ALMOST_INFINITY = 1e99;

    @Test
    void integrate_whenFirstDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration(1, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenSecondDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration(2, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenThirdDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration(3, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenFourthDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration(4, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenFifthDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration(5, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenSixthDegreePolynomial_returnsExpectedResult() throws IntegrationException {
        assertPolynomialIntegration(6, ABSOLUTE_ERROR_6);
    }

    @Test
    void integrate_whenGaussian_returnsExpectedResult() throws IntegrationException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var mu = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var sigma = ABSOLUTE_ERROR_GAUSSIAN + Math.abs(randomizer.nextDouble(a, MAX_VALUE));

        final var expected = NormalDist.cdf(b, mu, sigma) - NormalDist.cdf(a, mu, sigma);

        final var integrator = new RombergTrapezoidalQuadratureIntegrator(a, b,
                point -> NormalDist.p(point, mu, sigma));
        final var result = integrator.integrate();

        assertEquals(expected, result, ABSOLUTE_ERROR_GAUSSIAN);
    }

    @Test
    void integrate_whenExponential_returnsExpectedResult() throws IntegrationException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var lambda = randomizer.nextDouble(MIN_LAMBDA, MAX_LAMBDA);

        final var expected = 1.0 / lambda * (Math.exp(lambda * b) - Math.exp(lambda * a));

        final var integrator = new RombergTrapezoidalQuadratureIntegrator(a, b,
                point -> Math.exp(lambda * point));
        final var result = integrator.integrate();

        assertEquals(expected, result, ABSOLUTE_ERROR_EXPONENTIAL);
    }

    @Test
    void integrate_whenImproperIntegrandWithSingularities_throwsIntegrationException() {

        final var integrator = new RombergTrapezoidalQuadratureIntegrator(0.0, 1.0,
                point -> Math.log(point) * Math.log(1 - point));
        assertThrows(IntegrationException.class, integrator::integrate);
    }

    @Test
    void integrate_whenImproperIntegralFromZeroToInfinity3_throwsIntegrationException() {

        final var integrator = new RombergTrapezoidalQuadratureIntegrator(0.0, ALMOST_INFINITY,
                point -> Math.pow(point, -2.0 / 7.0) * Math.exp(-point * point));
        assertThrows(IntegrationException.class, integrator::integrate);
    }

    @Test
    void getIntegratorType_returnsExpectedValue() {
        final var integrator = new RombergTrapezoidalQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() {
        final var integrator = new RombergTrapezoidalQuadratureIntegrator(0.0, 1.0, null);
        assertEquals(QuadratureType.TRAPEZOIDAL, integrator.getQuadratureType());
    }

    private void assertPolynomialIntegration(final int degree, final double error) throws IntegrationException {
        final var polynomial = buildPolynomial(degree);
        final var integrationPolynomial = polynomial.integrationAndReturnNew();

        // set integration interval
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final var expected = integrationPolynomial.evaluate(b) - integrationPolynomial.evaluate(a);

        final var integrator = new RombergTrapezoidalQuadratureIntegrator(a, b, polynomial::evaluate);
        final var result = integrator.integrate();

        assertEquals(expected, result, error);
    }

    private static Polynomial buildPolynomial(final int degree) {
        final var randomizer = new UniformRandomizer();
        final var result = new Polynomial(1.0);
        for (var i = 0; i < degree; i++) {
            final var root = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            final var poly = new Polynomial(-root, 1.0);
            result.multiply(poly);
        }

        return result;
    }
}