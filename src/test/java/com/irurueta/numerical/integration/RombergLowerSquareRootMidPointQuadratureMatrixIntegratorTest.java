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
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.EvaluationException;
import com.irurueta.numerical.ExponentialMatrixEstimator;
import com.irurueta.numerical.polynomials.Polynomial;
import com.irurueta.statistics.NormalDist;
import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class RombergLowerSquareRootMidPointQuadratureMatrixIntegratorTest {

    private static final double MIN_VALUE = -10.0;

    private static final double MAX_VALUE = 10.0;

    private static final double MIN_LAMBDA = -1.0;

    private static final double MAX_LAMBDA = 1.0;

    private static final double ABSOLUTE_ERROR_1 = 1e-9;

    private static final double ABSOLUTE_ERROR_5 = 1e-1;

    private static final double ABSOLUTE_ERROR_GAUSSIAN = 1e-7;

    private static final double ABSOLUTE_ERROR_EXPONENTIAL = 1e-4;

    private static final double ABSOLUTE_ERROR_IMPROPER_1 = 1e-5;

    @Test
    void integrate_whenFirstDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration(1, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenSecondDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration(2, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenThirdDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration(3, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenFourthDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration(4, ABSOLUTE_ERROR_1);
    }

    @Test
    void integrate_whenFifthDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration(5, ABSOLUTE_ERROR_5);
    }

    @Test
    void integrate_whenSixthDegreePolynomial_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        assertPolynomialIntegration(6, ABSOLUTE_ERROR_5);
    }

    @Test
    void integrate_whenGaussian_returnsExpectedResult() throws IntegrationException, WrongSizeException {
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var mu = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var sigma = ABSOLUTE_ERROR_GAUSSIAN + Math.abs(randomizer.nextDouble(a, MAX_VALUE));

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, NormalDist.p(point, mu, sigma));
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(a, b, listener);

        final var expected = NormalDist.cdf(b, mu, sigma) - NormalDist.cdf(a, mu, sigma);

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_GAUSSIAN));
    }

    @Test
    void integrate_whenExponential1_returnsExpectedResult() throws WrongSizeException, IntegrationException {

        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);
        final var lambda = randomizer.nextDouble(MIN_LAMBDA, MAX_LAMBDA);

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double t, Matrix result) {
                result.setElementAtIndex(0, Math.exp(lambda * t));
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(a, b, listener);

        final var expected = 1.0 / lambda * (Math.exp(lambda * b) - Math.exp(lambda * a));

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_EXPONENTIAL));
    }

    @Test
    void integrate_whenExponential2_returnsExpectedResult() throws WrongSizeException, IntegrationException {
        // for single dimension functions, integral of f(t) = e^(a*t) is f(t) = 1/a*e^(a*t)
        // Consequently, for matrix function f(t) = e^(A*t) where A is a matrix, should be:
        // f(t) = A^-1 * e^(A*t)
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        // for matrix A = [0    1]
        //                [-1   0]

        // f(t) = e^(A*t) = [cos(t)     sin(t)]
        //                  [-sin(t)    cos(t)]

        // matrix A is antisymmetric and orthonormal, consequently, its inverse is its transpose:
        // A^-1 = A^T = [0  -1]
        //              [1   0]

        // A * A^T = [0   1][0  -1] = [0*0 + 1*1    -1*0 + 1*0 ] = [1  0]
        //           [-1  0][1   0]   [-1*0 + 0*1   -1*-1 + 0*0]   [0  1]

        // Consequently, the integral of f(t) is:
        // F(t) = A^-1*e^(A*t) = [0  -1][cos(t)     sin(t)] = [sin(t)  -cos(t)]
        //                       [1   0][-sin(t)    cos(t)]   [cos(t)   sin(t)]

        // Finally, integral of f(t) between a and b will be:
        // F(b) - F(a) = [sin(b) - sin(a)       -cos(b) + cos(a)]
        //               [cos(b) - cos(a)        sin(b) - sin(a)]

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double t, Matrix result) {
                // f(t) = e^(A*t) = [cost(t)  sin(t)]
                //                  [-sin(t)  cos(t)]
                final var c = Math.cos(t);
                final var s = Math.sin(t);
                result.setElementAtIndex(0, c);
                result.setElementAtIndex(1, -s);
                result.setElementAtIndex(2, s);
                result.setElementAtIndex(3, c);
            }

            @Override
            public int getRows() {
                return 2;
            }

            @Override
            public int getColumns() {
                return 2;
            }
        };

        // F(b) - F(a) = [sin(b) - sin(a)       -cos(b) + cos(a)]
        //               [cos(b) - cos(a)        sin(b) - sin(a)]
        final var expectedResult = new Matrix(2, 2);
        final var sinDiff = Math.sin(b) - Math.sin(a);
        final var cosDiff = Math.cos(b) - Math.cos(a);
        expectedResult.setElementAtIndex(0, sinDiff);
        expectedResult.setElementAtIndex(1, cosDiff);
        expectedResult.setElementAtIndex(2, -cosDiff);
        expectedResult.setElementAtIndex(3, sinDiff);

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(a, b, listener);

        final var integrationResult = new Matrix(2, 2);
        integrator.integrate(integrationResult);

        // check
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_EXPONENTIAL));
    }

    @Test
    void integrate_whenExponential3_returnsExpectedResult() throws WrongSizeException, IntegrationException {
        // for single dimension functions, integral of f(t) = e^(a*t) is f(t) = 1/a*e^(a*t)
        // Consequently, for matrix function f(t) = e^(A*t) where A is a matrix, should be:
        // f(t) = A^-1 * e^(A*t)
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        // for matrix A = [0    1]
        //                [-1   0]

        // f(t) = e^(A*t) = [cos(t)     sin(t)]
        //                  [-sin(t)    cos(t)]

        // matrix A is antisymmetric and orthonormal, consequently, its inverse is its transpose:
        // A^-1 = A^T = [0  -1]
        //              [1   0]

        // A * A^T = [0   1][0  -1] = [0*0 + 1*1    -1*0 + 1*0 ] = [1  0]
        //           [-1  0][1   0]   [-1*0 + 0*1   -1*-1 + 0*0]   [0  1]

        // Consequently, the integral of f(t) is:
        // F(t) = A^-1*e^(A*t) = [0  -1][cos(t)     sin(t)] = [sin(t)  -cos(t)]
        //                       [1   0][-sin(t)    cos(t)]   [cos(t)   sin(t)]

        // Finally, integral of f(t) between a and b will be:
        // F(b) - F(a) = [sin(b) - sin(a)       -cos(b) + cos(a)]
        //               [cos(b) - cos(a)        sin(b) - sin(a)]

        final var matrixA = new Matrix(2, 2);
        matrixA.setElementAtIndex(0, 0.0);
        matrixA.setElementAtIndex(1, -1.0);
        matrixA.setElementAtIndex(2, 1.0);
        matrixA.setElementAtIndex(3, 0.0);

        final var tmp = new Matrix(2, 2);

        final var exponentialMatrixEstimator = new ExponentialMatrixEstimator();

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double t, Matrix result) throws EvaluationException {
                try {
                    // f(t) = e^(A*t) = [cost(t)  sin(t)]
                    //                  [-sin(t)  cos(t)]
                    tmp.copyFrom(matrixA);
                    tmp.multiplyByScalar(t);
                    exponentialMatrixEstimator.exponential(tmp, result);
                } catch (final AlgebraException ex) {
                    throw new EvaluationException(ex);
                }
            }

            @Override
            public int getRows() {
                return 2;
            }

            @Override
            public int getColumns() {
                return 2;
            }
        };

        // F(b) - F(a) = [sin(b) - sin(a)       -cos(b) + cos(a)]
        //               [cos(b) - cos(a)        sin(b) - sin(a)]
        final var expectedResult = new Matrix(2, 2);
        final var sinDiff = Math.sin(b) - Math.sin(a);
        final var cosDiff = Math.cos(b) - Math.cos(a);
        expectedResult.setElementAtIndex(0, sinDiff);
        expectedResult.setElementAtIndex(1, cosDiff);
        expectedResult.setElementAtIndex(2, -cosDiff);
        expectedResult.setElementAtIndex(3, sinDiff);

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(a, b, listener);

        final var integrationResult = new Matrix(2, 2);
        integrator.integrate(integrationResult);

        // check
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_EXPONENTIAL));
    }

    @Test
    void integrate_whenImproperIntegrandWithSingularities_returnsExpectedResult() throws IntegrationException,
            WrongSizeException {

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, Math.log(point) * Math.log(1.0 - point));
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);

        final var expected = 2.0 - Math.PI * Math.PI / 6.0;

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, ABSOLUTE_ERROR_IMPROPER_1));
    }

    @Test
    void getIntegratorType_returnsExpectedValue() throws WrongSizeException {
        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                // no action needed
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(IntegratorType.ROMBERG, integrator.getIntegratorType());
    }

    @Test
    void getQuadratureType_returnsExpectedValue() throws WrongSizeException {
        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                // no action needed
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(0.0, 1.0, listener);
        assertEquals(QuadratureType.LOWER_SQUARE_ROOT_MID_POINT, integrator.getQuadratureType());
    }

    private static void assertPolynomialIntegration(final int degree, final double error) throws IntegrationException,
            WrongSizeException {
        final var polynomial = buildPolynomial(degree);
        final var integrationPolynomial = polynomial.integrationAndReturnNew();

        // set integration interval
        final var randomizer = new UniformRandomizer();
        final var a = randomizer.nextDouble(MIN_VALUE, MAX_VALUE);
        final var b = randomizer.nextDouble(a, MAX_VALUE);

        final var expected = integrationPolynomial.evaluate(b) - integrationPolynomial.evaluate(a);

        final var listener = new MatrixSingleDimensionFunctionEvaluatorListener() {

            @Override
            public void evaluate(double point, Matrix result) {
                result.setElementAtIndex(0, polynomial.evaluate(point));
            }

            @Override
            public int getRows() {
                return 1;
            }

            @Override
            public int getColumns() {
                return 1;
            }
        };

        final var integrator = new RombergLowerSquareRootMidPointQuadratureMatrixIntegrator(a, b, listener);

        final var integrationResult = new Matrix(1, 1);
        integrator.integrate(integrationResult);

        // check
        final var expectedResult = new Matrix(1, 1);
        expectedResult.setElementAtIndex(0, expected);
        assertTrue(expectedResult.equals(integrationResult, error));
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