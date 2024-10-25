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
import java.util.logging.Level;
import java.util.logging.Logger;

class RationalInterpolatorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final int SAMPLES = 1000;

    private static final int INTERPOLATIONS = 1000;

    private static final double ABSOLUTE_ERROR_1 = 1e-7;

    private static final double ABSOLUTE_ERROR_2 = 1e-1;

    private static final double ABSOLUTE_ERROR_3 = 1e-2;

    private static final double ABSOLUTE_ERROR_4 = 1.0;

    private static final double ABSOLUTE_ERROR_5 = 2.0;

    @Test
    void interpolate_whenFirstDegreePolynomial_returnsExpectedResult() {
        assertThrows(InterpolationException.class, () -> assertInterpolation(1, ABSOLUTE_ERROR_1));
    }

    @Test
    void interpolate_whenSecondDegreePolynomial_returnsExpectedResult() throws InterpolationException,
            SortingException {
        assertInterpolation(2, ABSOLUTE_ERROR_2);
    }

    @Test
    void interpolate_whenThirdDegreePolynomial_returnsExpectedResult() throws InterpolationException, SortingException {
        assertInterpolation(3, ABSOLUTE_ERROR_3);
    }

    @Test
    void interpolate_whenFourthDegreePolynomial_returnsExpectedResult() throws InterpolationException,
            SortingException {
        assertInterpolation(4, ABSOLUTE_ERROR_4);
    }

    @Test
    void interpolate_whenFifthDegreePolynomial_returnsExpectedResult() throws InterpolationException, SortingException {
        assertInterpolation(5, ABSOLUTE_ERROR_5);
    }

    @Test
    void interpolate_whenMismatchedLength_throwsIllegalArgumentException() {
        final var x = new double[SAMPLES];
        final var y = new double[SAMPLES + 1];
        assertThrows(IllegalArgumentException.class, () -> new RationalInterpolator(x, y));
    }

    @Test
    void interpolate_whenNotEnoughSamples_throwsIllegalArgumentException() {
        final var x = new double[1];
        final var y = new double[1];
        assertThrows(IllegalArgumentException.class, () -> new RationalInterpolator(x, y));
    }

    @Test
    void interpolate_whenNotEnoughSamplesToTakeIntoAccount_throwsIllegalArgumentException() {
        final var x = new double[SAMPLES];
        final var y = new double[SAMPLES];
        assertThrows(IllegalArgumentException.class, () -> new RationalInterpolator(x, y, 1));
    }

    @Test
    void interpolate_whenTooManyPointsToTakeIntoAccount_throwsIllegalArgumentException() {
        final var x = new double[SAMPLES];
        final var y = new double[SAMPLES];
        assertThrows(IllegalArgumentException.class, () -> new RationalInterpolator(x, y, SAMPLES + 1));
    }

    private static void assertInterpolation(final int degree, final double error) throws InterpolationException,
            SortingException {
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
        for (int i = 0; i < SAMPLES; i++) {
            assertEquals(y[i], numPolynomial.evaluate(x[i]) / denomPolynomial.evaluate(x[i]), 0.0);
        }

        final var interpolator = new RationalInterpolator(x, y);

        // check that interpolator passes through provided points
        for (var i = 0; i < SAMPLES; i++) {
            assertEquals(y[i], interpolator.interpolate(x[i]), 0.0);
        }

        // check random values
        var accuracy = 0.0;
        for (var i = 0; i < INTERPOLATIONS; i++) {
            final var xi = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            assertEquals(numPolynomial.evaluate(xi) / denomPolynomial.evaluate(xi),
                    interpolator.interpolate(xi), error);
            final var dy = interpolator.getDy();
            accuracy += dy / INTERPOLATIONS;
        }

        Logger.getLogger(PolynomialInterpolatorTest.class.getName())
                .log(Level.INFO, String.format("accuracy: %f degree: %d", accuracy, degree));
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