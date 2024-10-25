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
package com.irurueta.numerical.interpolation;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class AccurateInterpolatingPolynomialEstimatorTest {

    private static final double MIN_VALUE = -1.0;

    private static final double MAX_VALUE = 1.0;

    private static final double ABSOLUTE_ERROR = 1e-8;

    @Test
    void estimate_whenFirstDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimation(1);
    }

    @Test
    void estimate_whenSecondDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimation(2);
    }

    @Test
    void estimate_whenThirdDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimation(3);
    }

    @Test
    void estimate_whenFourthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimation(4);
    }

    @Test
    void estimate_whenFifthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimation(5);
    }

    @Test
    void estimate_whenSixthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimation(6);
    }

    @Test
    void estimateCoefficients_whenFirstDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(1);
    }

    @Test
    void estimateCoefficients_whenSecondDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(2);
    }

    @Test
    void estimateCoefficients_whenThirdDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(3);
    }

    @Test
    void estimateCoefficients_whenFourthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(4);
    }

    @Test
    void estimateCoefficients_whenFifthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(5);
    }

    @Test
    void estimateCoefficients_whenSixthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(6);
    }

    private static void assertEstimation(final int degree) throws InterpolationException {
        final var polynomial = buildPolynomial(degree);

        final var samples = degree + 1;
        final var x = new double[samples];
        final var y = new double[samples];

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < samples; i++) {
            x[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            y[i] = polynomial.evaluate(x[i]);
        }

        final var estimator = new AccurateInterpolatingPolynomialEstimator();

        final var result = new double[samples];
        estimator.estimate(x, y, result);

        assertArrayEquals(polynomial.getPolyParams(), result, ABSOLUTE_ERROR);
    }

    private static void assertEstimationCoefficients(final int degree) throws InterpolationException {
        final var polynomial = buildPolynomial(degree);

        final var samples = degree + 1;
        final var x = new double[samples];
        final var y = new double[samples];

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < samples; i++) {
            x[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            y[i] = polynomial.evaluate(x[i]);
        }

        final var estimator = new AccurateInterpolatingPolynomialEstimator();

        final var result = estimator.estimateCoefficients(x, y);

        assertArrayEquals(polynomial.getPolyParams(), result, ABSOLUTE_ERROR);
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