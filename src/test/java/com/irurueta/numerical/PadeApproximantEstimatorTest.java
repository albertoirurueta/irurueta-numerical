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
package com.irurueta.numerical;

import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class PadeApproximantEstimatorTest {

    private static final double[] TAYLOR_COEFFICIENTS = {
            2.0, 1.0 / 9.0, 1.0 / 81.0, -49.0 / 8748.0, 175.0 / 78732
    };

    private static final double MIN_VALUE = 1.0;
    private static final double MAX_VALUE = 10.0;

    private static final double MIN_EXP_VALUE = -2.0;

    private static final double MAX_EXP_VALUE = 2.0;

    private static final int EXPONENTIAL_TAYLOR_ORDER = 14;

    @Test
    void constructor_whenInvalidNumberOfTimes_throwsIllegalArgumentException() {
        assertThrows(IllegalArgumentException.class, () -> new PadeApproximantEstimator(-1));
    }

    @Test
    void estimatePadeCoefficients_whenInvalidTaylorSeries_throwsIllegalArgumentException() {
        final var estimator = new PadeApproximantEstimator();
        assertThrows(IllegalArgumentException.class, () -> estimator.estimatePadeCoefficients(new double[1]));
    }

    @Test
    void estimatePadeCoefficients_whenInvalidNumeratorCoefficients_throwsException() {
        final var estimator = new PadeApproximantEstimator();
        assertThrows(IllegalArgumentException.class,
                () -> estimator.estimatePadeCoefficients(TAYLOR_COEFFICIENTS, new double[1], new double[3]));
    }

    @Test
    void estimatePadeCoefficients_whenInvalidDenominatorCoefficients_throwsException() {
        final var estimator = new PadeApproximantEstimator();
        assertThrows(IllegalArgumentException.class,
                () -> estimator.estimatePadeCoefficients(TAYLOR_COEFFICIENTS, new double[3], new double[1]));
    }

    @Test
    void estimatePadeCoefficients_whenTestFunction_returnsExpectedValue() throws NumericalException {
        final var estimator = new PadeApproximantEstimator();
        final var result = estimator.estimatePadeCoefficients(TAYLOR_COEFFICIENTS);
        final var numerators = new double[3];
        final var denominators = new double[3];
        estimator.estimatePadeCoefficients(TAYLOR_COEFFICIENTS, numerators, denominators);

        assertArrayEquals(numerators, result.getNumerators(), 0.0);
        assertArrayEquals(denominators, result.getDenominators(), 0.0);

        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var exactValue = testFunction(x);
        final var taylorValue = testTaylorApproximation(x, TAYLOR_COEFFICIENTS);
        final var padeValue = padeApproximation(x, numerators, denominators);

        final var taylorError = Math.abs(taylorValue - exactValue);
        final var padeError = Math.abs(padeValue - exactValue);

        assertTrue(padeError < taylorError);
    }

    @Test
    void estimatePadeCoefficients_whenExponentialFunction_returnsExpectedValue() throws NumericalException {
        final var taylorCoefficients = exponentialTaylorCoefficients();

        final var estimator = new PadeApproximantEstimator();
        final var result = estimator.estimatePadeCoefficients(taylorCoefficients);

        final var randomizer = new UniformRandomizer();
        final var x = randomizer.nextDouble(MIN_EXP_VALUE, MAX_EXP_VALUE);
        final var exactValue = Math.exp(x);
        final var taylorValue = testTaylorApproximation(x, taylorCoefficients);
        final var padeValue = padeApproximation(x, result.getNumerators(), result.getDenominators());

        final var taylorError = Math.abs(taylorValue - exactValue);
        final var padeError = Math.abs(padeValue - exactValue);

        assertTrue(taylorError >= 0.0);
        assertTrue(padeError >= 0.0);
    }

    private static double[] exponentialTaylorCoefficients() {
        final var factorialEstimator = new DoubleFactorialEstimator();
        final var result = new double[EXPONENTIAL_TAYLOR_ORDER];
        for (var i = 0; i < EXPONENTIAL_TAYLOR_ORDER; i++) {
            result[i] = 1.0 / factorialEstimator.factorial(i);
        }

        return result;
    }

    private static double testFunction(final double x) {
        return Math.pow(7.0 + Math.pow(1 + x, 4.0 / 3.0), 1.0 / 3.0);
    }

    private static double testTaylorApproximation(final double x, final double[] taylorCoefficients) {
        var result = 0.0;
        var prevX = 1.0;
        for (final var taylorCoefficient : taylorCoefficients) {
            result += taylorCoefficient * prevX;
            prevX *= x;
        }

        return result;
    }

    private static double padeApproximation(
            final double x, double[] numeratorCoefficients, double[] denominatorCoefficients) {
        var numerator = 0.0;
        var prevXNumerator = 1.0;
        for (final var numeratorCoefficient : numeratorCoefficients) {
            numerator += numeratorCoefficient * prevXNumerator;
            prevXNumerator *= x;
        }

        var denominator = 0.0;
        var prevXDenominator = 1.0;
        for (final var denominatorCoefficient : denominatorCoefficients) {
            denominator += denominatorCoefficient * prevXDenominator;
            prevXDenominator *= x;
        }

        return numerator / denominator;
    }
}