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

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.ArrayUtils;
import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.SingularValueDecomposer;
import com.irurueta.algebra.Utils;

/**
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtained.
 * This implementation is more robust against small discontinuities than
 * SymmetricDerivativeEstimator, but it is also slower to compute.
 * This method interpolates the sampled function values into a polynomial of
 * 2nd degree (parabolic), whose derivative is known.
 * Because a linear system of equations has to be solved to determine such
 * polynomial, this method might be less accurate when large values are involved
 * due to limited machine precision.
 */
@SuppressWarnings("WeakerAccess")
public class SavitzkyGolayDerivativeEstimator extends DerivativeEstimator {

    /**
     * Number of required point to evaluate to compute derivative.
     */
    public static final int N_POINTS = 3;

    /**
     * Constructor.
     *
     * @param listener listener to evaluate a single dimension function.
     */
    public SavitzkyGolayDerivativeEstimator(
            final SingleDimensionFunctionEvaluatorListener listener) {
        super(listener);
    }

    /**
     * Computes the function derivative at provided point x.
     *
     * @param x Point where derivative is estimated.
     * @return Derivative of function at provided point.
     * @throws EvaluationException Raised if function cannot be properly
     *                             evaluated.
     */
    @Override
    @SuppressWarnings("Duplicates")
    public double derivative(final double x) throws EvaluationException {
        // fit a polynomial of degree 2 by evaluating function at x-h, x and x+h
        double h = EPS * Math.abs(x);
        if (h == 0.0) {
            // Trick to reduce finite-precision error
            h = EPS;
        }

        final double xh1 = x + h;
        final double xh2 = x - h;

        final double f = listener.evaluate(x);
        final double fh1 = listener.evaluate(xh1);
        final double fh2 = listener.evaluate(xh2);

        // express the problem as:
        // a * x^2 + b * x + c = f(x)
        // b * xh1^2 + b * xh1 + c = f(xh1)
        // c * xh2^2 + b * xh2 + c = f(xh2)

        final Matrix a;
        final double aParam;
        final double bParam;
        try {
            a = new Matrix(N_POINTS, N_POINTS);

            a.setElementAt(0, 0, x * x);
            a.setElementAt(1, 0, xh1 * xh1);
            a.setElementAt(2, 0, xh2 * xh2);

            a.setElementAt(0, 1, x);
            a.setElementAt(1, 1, xh1);
            a.setElementAt(2, 1, xh2);

            a.setElementAt(0, 2, 1.0);
            a.setElementAt(1, 2, 1.0);
            a.setElementAt(2, 2, 1.0);

            final double[] b = new double[N_POINTS];


            // normalize to increase accuracy
            final double normA = Utils.normF(a);
            a.multiplyByScalar(1.0 / normA);

            b[0] = f;
            b[1] = fh1;
            b[2] = fh2;

            // normalize to increase accuracy
            ArrayUtils.multiplyByScalar(b, 1.0 / normA, b);

            final SingularValueDecomposer decomposer = new SingularValueDecomposer(a);

            decomposer.decompose();

            // now solve the system of equations in Least Mean Squared Error
            // because SVD allows the system of equations to be solved using the
            // pseudo-inverse
            final double[] params = decomposer.solve(b);
            aParam = params[0];
            bParam = params[1];

        } catch (final AlgebraException e) {
            return Double.NaN;
        }

        // and c = params[2], but we don't need it

        // because we have fitted the function into a polynomial that has
        // expression: a * x^2 + b * x + c, then its derivative is:
        // 2.0 * a * x + b, therefore:
        return 2.0 * aParam * x + bParam;
    }
}
