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
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class ShepardInterpolatorTest {

    private static final double MIN_VALUE = -1.0;

    private static final double MAX_VALUE = 1.0;

    private static final int SAMPLES = 10000;

    private static final double ABSOLUTE_ERROR_1 = 1e-2;

    private static final double ABSOLUTE_ERROR_2 = 2.0;

    @Test
    void interpolate_dim1_returnsExpectedResult() throws WrongSizeException {
        assertInterpolation(1, ABSOLUTE_ERROR_1);
    }

    @Test
    void interpolate_dim2_returnsExpectedResult() throws WrongSizeException {
        assertInterpolation(2, ABSOLUTE_ERROR_2);
    }

    private static void assertInterpolation(final int dim, final double error) throws WrongSizeException {
        final var roots = new double[dim];
        final var polynomials = buildPolynomials(dim, roots);

        for (var i = 0; i < dim; i++) {
            assertEquals(0.0, polynomials[i].evaluate(roots[i]), 0.0);
        }
        assertEquals(0.0, evaluate(polynomials, roots), 0.0);

        // create multiple samples and evaluations
        final var randomizer = new UniformRandomizer();
        final var point = new double[dim];
        final var pts = new Matrix(SAMPLES, dim);
        final var values = new double[SAMPLES];
        for (var i = 0; i < SAMPLES; i++) {
            randomizer.fill(point, MIN_VALUE, MAX_VALUE);
            pts.setSubmatrix(i, 0, i, dim -1, point);
            values[i] = evaluate(polynomials, point);
        }

        // check data
        for (var i = 0; i < SAMPLES; i++) {
            pts.getSubmatrixAsArray(i, 0, i, dim - 1, point);
            assertEquals(values[i], evaluate(polynomials, point), 0.0);
        }

        final var interpolator = new ShepardInterpolator(pts, values);

        // check that interpolator at provided points
        for (var i = 0; i < SAMPLES; i++) {
            pts.getSubmatrixAsArray(i, 0, i, dim - 1, point);
            assertEquals(values[i], interpolator.interpolate(point), error);
        }

        // check random values
        for (var i = 0; i < SAMPLES; i++) {
            randomizer.fill(point, MIN_VALUE, MAX_VALUE);
            assertEquals(evaluate(polynomials, point), interpolator.interpolate(point), error);
        }
    }

    private static double evaluate(final Polynomial[] polynomials, double[] point) {
        final var dim = polynomials.length;
        var result = 1.0;
        for(var i = 0; i < dim; i++) {
            result *= polynomials[i].evaluate(point[i]);
        }
        return result;
    }

    private static Polynomial[] buildPolynomials(int dim, final double[] roots) {
        final var r = new double[1];
        final var result = new Polynomial[dim];
        for (var i = 0; i < dim; i++) {
            result[i] = buildPolynomial(r);
            roots[i] = r[0];
        }
        return result;
    }

    private static Polynomial buildPolynomial(final double[] roots) {
        final var randomizer = new UniformRandomizer();
        final var root = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var result = new Polynomial(-root, 1.0);
        roots[0] = root;

        return result;
    }
}