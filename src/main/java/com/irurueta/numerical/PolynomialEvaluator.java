/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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

import com.irurueta.algebra.Complex;

/**
 * Utility class to evaluate polynomials having either real or complex
 * coefficients.
 */
public class PolynomialEvaluator {

    /**
     * Empty constructor.
     */
    protected PolynomialEvaluator() {
    }

    /**
     * Evaluates polynomial formed by provided polynomial parameters at provided
     * point x. The polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * hence, the array of parameters is [a0, a1, ... a(n-1), an].
     *
     * @param polyParams array of polynomial parameters.
     * @param x          point where polynomial is evaluated.
     * @return result of evaluation.
     * @throws IllegalArgumentException if provided array is null or has length
     *                                  0.
     */
    public static double evaluate(final double[] polyParams, final double x) {
        if (polyParams == null || polyParams.length == 0) {
            throw new IllegalArgumentException();
        }

        final int length = polyParams.length;

        double result = 0.0;
        double powX = 1.0;
        for (int i = length - 1; i >= 0; i--) {
            result += polyParams[i] * powX;
            powX *= x;
        }

        return result;
    }

    /**
     * Evaluates polynomial formed by provided polynomial parameters at provided
     * point x. The polynomial of degree n is defined as:
     * p(x) = a0 * x^n + a1 * x^(n - 1) + ... a(n-1) * x + an
     * hence, the array of parameters is [a0, a1, ... a(n-1), an].
     *
     * @param polyParams array of polynomial parameters.
     * @param x          point where polynomial is evaluated.
     * @return result of evaluation.
     * @throws IllegalArgumentException if provided array is null or has length
     *                                  0.
     */
    public static Complex evaluate(final Complex[] polyParams, final Complex x) {
        if (polyParams == null || polyParams.length == 0) {
            throw new IllegalArgumentException();
        }

        final int length = polyParams.length;

        final Complex result = new Complex();
        final Complex powX = new Complex(1.0, 0.0);
        final Complex tmp = new Complex();
        for (int i = length - 1; i >= 0; i--) {

            polyParams[i].multiply(powX, tmp);
            result.add(tmp);

            powX.multiply(x);
        }

        return result;
    }
}
