package com.irurueta.numerical.interpolation;

import static org.junit.Assert.assertArrayEquals;

import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.Test;

public class SimpleInterpolatingPolynomialEstimatorTest {

    private static final double MIN_VALUE = -1.0;

    private static final double MAX_VALUE = 1.0;

    private static final double ABSOLUTE_ERROR = 1e-5;

    @Test
    public void estimate_whenFirstDegree_returnsExpectedResult() {
        assertEstimation(1);
    }

    @Test
    public void estimate_whenSecondDegree_returnsExpectedResult() {
        assertEstimation(2);
    }

    @Test
    public void estimate_whenThirdDegree_returnsExpectedResult() {
        assertEstimation(3);
    }

    @Test
    public void estimate_whenFourthDegree_returnsExpectedResult() {
        assertEstimation(4);
    }

    @Test
    public void estimateCoefficients_whenFirstDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(1);
    }

    @Test
    public void estimateCoefficients_whenSecondDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(2);
    }

    @Test
    public void estimateCoefficients_whenThirdDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(3);
    }

    @Test
    public void estimateCoefficients_whenFourthDegree_returnsExpectedResult() throws InterpolationException {
        assertEstimationCoefficients(4);
    }

    private static void assertEstimation(final int degree) {
        final Polynomial polynomial = buildPolynomial(degree);

        final int samples = degree + 1;
        final double[] x = new double[samples];
        final double[] y = new double[samples];

        final UniformRandomizer randomizer = new UniformRandomizer();
        for (int i = 0; i < samples; i++) {
            x[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            y[i] = polynomial.evaluate(x[i]);
        }

        final SimpleInterpolatingPolynomialEstimator estimator =
                new SimpleInterpolatingPolynomialEstimator();

        final double[] result = new double[samples];
        estimator.estimate(x, y, result);

        assertArrayEquals(polynomial.getPolyParams(), result, ABSOLUTE_ERROR);
    }

    private static void assertEstimationCoefficients(final int degree) throws InterpolationException {
        final Polynomial polynomial = buildPolynomial(degree);

        final int samples = degree + 1;
        final double[] x = new double[samples];
        final double[] y = new double[samples];

        final UniformRandomizer randomizer = new UniformRandomizer();
        for (int i = 0; i < samples; i++) {
            x[i] = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            y[i] = polynomial.evaluate(x[i]);
        }

        final SimpleInterpolatingPolynomialEstimator estimator =
                new SimpleInterpolatingPolynomialEstimator();

        final double[] result = estimator.estimateCoefficients(x, y);

        assertArrayEquals(polynomial.getPolyParams(), result, ABSOLUTE_ERROR);
    }

    private static Polynomial buildPolynomial(final int degree) {
        final UniformRandomizer randomizer = new UniformRandomizer();
        final Polynomial result = new Polynomial(1.0);
        for (int i = 0; i < degree; i++) {
            final double root = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
            final Polynomial poly = new Polynomial(-root, 1.0);
            result.multiply(poly);
        }

        return result;
    }
}