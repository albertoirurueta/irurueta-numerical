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

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.sorting.Sorter;
import com.irurueta.sorting.SortingException;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

class BilinearInterpolatorTest {

    private static final double MIN_VALUE = -100.0;

    private static final double MAX_VALUE = 100.0;

    private static final int SAMPLES = 1000;

    private static final int DEGREE = 1;

    private static final int MIN_SAMPLES = 2;

    private static final double ABSOLUTE_ERROR = 1e-8;

    @Test
    void interpolate_whenFirstDegree2DPolynomial_returnsExpectedResult() throws SortingException, WrongSizeException {
        assertInterpolation(SAMPLES);
    }

    @Test
    void interpolate_whenFirstDegree2DPolynomialMinimumSamples_returnsExpectedResult() throws SortingException,
            WrongSizeException {
        assertInterpolation(MIN_SAMPLES);
    }

    @Test
    void interpolate_whenNotEnoughSamples_throwsIllegalArgumentException() throws WrongSizeException {
        final var x1 = new double[1];
        final var x2 = new double[1];
        final var y = new Matrix(SAMPLES, SAMPLES);
        assertThrows(IllegalArgumentException.class, () -> new BilinearInterpolator(x1, x2, y));
    }

    @Test
    void getM_returnsExpectedValue() throws WrongSizeException {
        final var x1 = new double[SAMPLES];
        final var x2 = new double[SAMPLES + 1];
        final var y = new Matrix(SAMPLES, SAMPLES + 1);
        final var interpolator = new BilinearInterpolator(x1, x2, y);

        assertEquals(SAMPLES, interpolator.getM());
    }

    @Test
    void getN_returnsExpectedValue() throws WrongSizeException {
        final var x1 = new double[SAMPLES];
        final var x2 = new double[SAMPLES + 1];
        final var y = new Matrix(SAMPLES, SAMPLES + 1);
        final var interpolator = new BilinearInterpolator(x1, x2, y);

        assertEquals(SAMPLES + 1, interpolator.getN());
    }

    private static void assertInterpolation(final int samples) throws SortingException, WrongSizeException {
        final var roots1 = new double[DEGREE];
        final var polynomial1 = buildPolynomial(roots1);

        final var roots2 = new double[DEGREE];
        final var polynomial2 = buildPolynomial(roots2);

        for (var i = 0; i < DEGREE; i++) {
            assertEquals(0.0, polynomial1.evaluate(roots1[i]), 0.0);
            assertEquals(0.0, polynomial2.evaluate(roots2[i]), 0.0);
        }

        // create multiple samples and evaluations
        final var randomizer = new UniformRandomizer();
        final var unorderedX1 = new double[samples];
        final var unorderedX2 = new double[samples];
        for (var i = 0; i < samples; i++) {
            unorderedX1[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            unorderedX2[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        }

        // order values by ascending x and create evaluations
        final var x1 = Arrays.copyOf(unorderedX1, samples);
        final var x2 = Arrays.copyOf(unorderedX2, samples);
        final var y = new Matrix(samples, samples);
        final var sorter = Sorter.<Double>create();
        sorter.sort(x1);
        sorter.sort(x2);
        for (var i = 0; i < samples; i++) {
            for (var j = 0; j < samples; j++) {
                y.setElementAt(i, j, polynomial1.evaluate(x1[i]) * polynomial2.evaluate(x2[j]));
            }
        }

        // check data
        for (var i = 0; i < samples; i++) {
            for (var j = 0; j < samples; j++) {
                assertEquals(y.getElementAt(i, j), polynomial1.evaluate(x1[i]) * polynomial2.evaluate(x2[j]),
                        0.0);
            }
        }

        final var interpolator = new BilinearInterpolator(x1, x2, y);

        assertEquals(samples, interpolator.getM());
        assertEquals(samples, interpolator.getN());

        // check that interpolator passes through provided points
        for (var i = 0; i < samples; i++) {
            for (var j = 0; j < samples; j++) {
                assertEquals(y.getElementAt(i, j), interpolator.interpolate(x1[i], x2[j]), 0.0);
            }
        }

        // check random values
        for (var i = 0; i < SAMPLES; i++) {
            final var xi = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            for (var j = 0; j < SAMPLES; j++) {
                final var xj = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
                assertEquals(polynomial1.evaluate(xi) * polynomial2.evaluate(xj),
                        interpolator.interpolate(xi, xj), ABSOLUTE_ERROR);
            }
        }
    }

    private static Polynomial buildPolynomial(final double[] roots) {
        final var randomizer = new UniformRandomizer();
        final var result = new Polynomial(1.0);
        for (var i = 0; i < DEGREE; i++) {
            final var root = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            final var poly = new Polynomial(-root, 1.0);
            result.multiply(poly);

            roots[i] = root;
        }

        return result;
    }
}