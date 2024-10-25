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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.sorting.Sorter;
import com.irurueta.sorting.SortingException;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

 class CubicSplineInterpolatorTest {

    private static final double MIN_VALUE = -100.0;

    private static final double MAX_VALUE = 100.0;

    private static final int SAMPLES = 1000;

    private static final int INTERPOLATIONS = 1000;

    private static final double ABSOLUTE_ERROR_1 = 1e-7;

    private static final double ABSOLUTE_ERROR_2 = 20.0;

    private static final double ABSOLUTE_ERROR_3 = 1e-1;

    @Test
    void interpolate_whenFirstDegreePolynomial_returnsExpectedResult() throws InterpolationException, SortingException {
        assertInterpolation(1, SAMPLES, ABSOLUTE_ERROR_1);
    }

    @Test
    void interpolate_whenFirstDegreePolynomialMinimumSamples_returnsExpectedResult() throws InterpolationException,
            SortingException {
        assertInterpolation(1, 2, ABSOLUTE_ERROR_1);
    }

    @Test
    void interpolate_whenSecondDegreePolynomial_returnsExpectedResult() throws InterpolationException,
            SortingException {
        assertInterpolation(2, SAMPLES, ABSOLUTE_ERROR_2);
    }

    @Test
    void interpolate_whenSecondDegreePolynomialMinimumSamples_returnsExpectedResult() throws InterpolationException,
            SortingException {
        assertInterpolation(2, 3, ABSOLUTE_ERROR_2);
    }

    @Test
    void interpolate_whenThirdDegreePolynomial_returnsExpectedResult() throws InterpolationException, SortingException {
        assertInterpolation(3, SAMPLES, ABSOLUTE_ERROR_3);
    }

    @Test
    void interpolate_whenThirdDegreePolynomialMinimumSamples_returnsExpectedResult() throws InterpolationException,
            SortingException {
        assertInterpolation(3, 4, ABSOLUTE_ERROR_3);
    }

    @Test
    void interpolate_whenMismatchedLength_throwsIllegalArgumentException() {
        final var x = new double[SAMPLES];
        final var y = new double[SAMPLES + 1];
        assertThrows(IllegalArgumentException.class, () -> new CubicSplineInterpolator(x, y));
    }

    @Test
    void interpolate_whenNotEnoughSamples_throwsIllegalArgumentException() {
        final var x = new double[1];
        final var y = new double[1];
        assertThrows(IllegalArgumentException.class, () -> new CubicSplineInterpolator(x, y));
    }

    private static void assertInterpolation(final int degree, final int samples, final double error)
            throws InterpolationException, SortingException {
        final var roots = new double[degree];
        final var polynomial = buildPolynomial(degree, roots);

        for (var i = 0; i < degree; i++) {
            assertEquals(0.0, polynomial.evaluate(roots[i]), error);
        }

        // create multiple samples and evaluations
        final var randomizer = new UniformRandomizer();
        final var unorderedX = new double[samples];
        final var unorderedY = new double[samples];
        for (var i = 0; i < samples; i++) {
            unorderedX[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            unorderedY[i] = polynomial.evaluate(unorderedX[i]);
        }

        // order values by ascending x
        final var x = Arrays.copyOf(unorderedX, samples);
        final var y = new double[samples];
        final var indices = Sorter.create().sortWithIndices(x);
        for (var i = 0; i < samples; i++) {
            y[i] = unorderedY[indices[i]];
        }

        // check data
        for (var i = 0; i < samples; i++) {
            assertEquals(y[i], polynomial.evaluate(x[i]), 0.0);
        }

        // get derivatives at lowest and highest endpoints
        final var derivativePolynomial = polynomial.derivativeAndReturnNew();
        final var yp1 = derivativePolynomial.evaluate(x[0]);
        final var ypn = derivativePolynomial.evaluate(x[samples - 1]);

        final var interpolator = new CubicSplineInterpolator(x, y, yp1, ypn);

        // check that interpolator passes through provided points
        for (var i = 0; i < samples; i++) {
            assertEquals(y[i], interpolator.interpolate(x[i]), 0.0);
        }

        // check random values
        for (var i = 0; i < INTERPOLATIONS; i++) {
            final var xi = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            assertEquals(polynomial.evaluate(xi), interpolator.interpolate(xi), error);
        }
    }

    private static Polynomial buildPolynomial(final int degree, final double[] roots) {
        final var randomizer = new UniformRandomizer();
        final var result = new Polynomial(1.0);
        for (var i = 0; i < degree; i++) {
            final var root = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            final var poly = new Polynomial(-root, 1.0);
            result.multiply(poly);

            roots[i] = root;
        }

        return result;
    }
}