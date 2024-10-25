/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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

/**
 * Class to estimate the gradient of a multidimensional function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the gradient at such point
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtained.
 * This implementation is more accurate although slower than
 * GradientEstimator.
 */
public class SymmetricGradientEstimator extends GradientEstimator {

    /**
     * Internal array containing one point to sample close to the original one.
     */
    private double[] xh1;

    /**
     * Internal array containing one point to sample close to the original one.
     */
    private double[] xh2;

    /**
     * Constructor.
     *
     * @param listener Listener to evaluate a multidimensional function.
     */
    public SymmetricGradientEstimator(final MultiDimensionFunctionEvaluatorListener listener) {
        super(listener);
    }

    /**
     * Sets estimated gradient in provided result array of a multidimensional
     * function at provided point.
     * This method is preferred respect to gradient(double[]) because result
     * array can be reused and hence is more memory efficient.
     *
     * @param point  Input point.
     * @param result Output parameter containing estimated array. This parameter
     *               must be an array of length equal to point.
     * @throws EvaluationException      Raised if function cannot be evaluated.
     * @throws IllegalArgumentException Raised if length of result and point are
     *                                  not equal.
     */
    @Override
    public void gradient(final double[] point, final double[] result) throws EvaluationException {
        final var n = point.length;
        if (result.length != n) {
            throw new IllegalArgumentException();
        }

        if (xh1 == null || xh1.length != n) {
            xh1 = new double[n];
            System.arraycopy(point, 0, xh1, 0, n);
        }
        if (xh2 == null || xh2.length != n) {
            xh2 = new double[n];
            System.arraycopy(point, 0, xh2, 0, n);
        }

        double temp;
        double h;
        double h1;
        double h2;
        double hh;
        double fh1;
        double fh2;

        for (var j = 0; j < n; j++) {
            temp = point[j];
            h = EPS * Math.abs(temp);
            if (h == 0.0) {
                // Trick to reduce finite-precision error
                h = EPS;
            }
            xh1[j] = temp + h;
            xh2[j] = temp - h;
            // because of machine precision h could be different in both cases

            h1 = xh1[j] - temp;
            h2 = temp - xh2[j];

            // this is more or less equal to 2.0 * h
            hh = h1 + h2;

            fh1 = listener.evaluate(xh1);
            fh2 = listener.evaluate(xh2);

            xh1[j] = temp;
            xh2[j] = temp;

            result[j] = (fh1 - fh2) / hh;
        }
    }
}
