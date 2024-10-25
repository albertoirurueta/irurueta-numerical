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

import static org.junit.jupiter.api.Assertions.*;

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.sorting.Sorter;
import com.irurueta.sorting.SortingException;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

class BarycentricRationalInterpolatorTest {

    private static final double MIN_VALUE = -100.0;

    private static final double MAX_VALUE = 100.0;

    private static final int SAMPLES = 1000;

    @Test
    void interpolate_whenFirstDegreePolynomial_returnsExpectedResult() throws SortingException {
        assertInterpolation(1);
    }

    @Test
    void interpolate_whenSecondDegreePolynomial_returnsExpectedResult() throws SortingException {
        assertInterpolation(2);
    }

    @Test
    void interpolate_whenThirdDegreePolynomial_returnsExpectedResult() throws SortingException {
        assertInterpolation(3);
    }

    @Test
    void interpolate_whenFourthDegreePolynomial_returnsExpectedResult() throws SortingException {
        assertInterpolation(4);
    }

    @Test
    void interpolate_whenFifthDegreePolynomial_returnsExpectedResult() throws SortingException {
        assertInterpolation(5);
    }

    @Test
    void interpolate_whenMismatchedLength_throwsIllegalArgumentException() {
        final var x = new double[SAMPLES];
        final var y = new double[SAMPLES + 1];
        assertThrows(IllegalArgumentException.class, () -> new BarycentricRationalInterpolator(x, y, 1));
    }

    @Test
    void interpolate_whenNotEnoughSamples_throwsIllegalArgumentException() {
        final var x = new double[1];
        final var y = new double[1];
        assertThrows(IllegalArgumentException.class, () -> new BarycentricRationalInterpolator(x, y, 1));
    }

    @Test
    void interpolate_whenTooManyPointsToTakeIntoAccount_throwsIllegalArgumentException() {
        final var x = new double[SAMPLES];
        final var y = new double[SAMPLES];
        assertThrows(IllegalArgumentException.class, () -> new BarycentricRationalInterpolator(x, y, SAMPLES + 1));
    }

    private static void assertInterpolation(final int degree) throws SortingException {
        final var roots = new double[degree];
        final var poles = new double[degree];
        final var numPolynomial = buildPolynomial(degree, roots);
        final var denomPolynomial = buildPolynomial(degree, poles);

        // create multiple samples and evaluations
        final var randomizer = new UniformRandomizer();
        final var unorderedX = new double[SAMPLES];
        final var unorderedY = new double[SAMPLES];
        for (var i = 0; i < SAMPLES; i++) {
            unorderedX[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            unorderedY[i] = numPolynomial.evaluate(unorderedX[i]) / denomPolynomial.evaluate(unorderedX[i]);
        }

        // order values by ascending x
        final var x = Arrays.copyOf(unorderedX, SAMPLES);
        final var y = new double[SAMPLES];
        final var indices = Sorter.create().sortWithIndices(x);
        for (var i = 0; i < SAMPLES; i++) {
            y[i] = unorderedY[indices[i]];
        }

        // check data
        for (var i = 0; i < SAMPLES; i++) {
            assertEquals(y[i], numPolynomial.evaluate(x[i]) / denomPolynomial.evaluate(x[i]), 0.0);
        }

        final var interpolator = new BarycentricRationalInterpolator(x, y, degree + 1);

        // check that interpolator passes through provided points
        for (var i = 0; i < SAMPLES; i++) {
            assertEquals(y[i], interpolator.interpolate(x[i]), 0.0);
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