package com.irurueta.numerical.interpolation;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class SimpleInterpolatingPolynomialEstimatorTest {

    private static final double MIN_VALUE = -1.0;

    private static final double MAX_VALUE = 1.0;

    private static final double ABSOLUTE_ERROR = 1e-5;

    @Test
    void estimate_whenFirstDegree_returnsExpectedResult() {
        assertEstimation(1);
    }

    @Test
    void estimate_whenSecondDegree_returnsExpectedResult() {
        assertEstimation(2);
    }

    @Test
    void estimate_whenThirdDegree_returnsExpectedResult() {
        assertEstimation(3);
    }

    @Test
    void estimate_whenFourthDegree_returnsExpectedResult() {
        assertEstimation(4);
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

    private static void assertEstimation(final int degree) {
        final var polynomial = buildPolynomial(degree);

        final var samples = degree + 1;
        final var x = new double[samples];
        final var y = new double[samples];

        final var randomizer = new UniformRandomizer();
        for (var i = 0; i < samples; i++) {
            x[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            y[i] = polynomial.evaluate(x[i]);
        }

        final var estimator = new SimpleInterpolatingPolynomialEstimator();

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

        final var estimator = new SimpleInterpolatingPolynomialEstimator();

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