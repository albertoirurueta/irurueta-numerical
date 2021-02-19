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
 * Class to estimate the derivative of a single dimension function at a given
 * point.
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is faster although less accurate than
 * SymmetricDerivativeEstimator.
 */
@SuppressWarnings("WeakerAccess")
public class DerivativeEstimator {

    /**
     * Constant defining machine precision for this algorithm.
     */
    public static final double EPS = 1e-8;

    /**
     * Listener to evaluate a single dimension function.
     */
    protected SingleDimensionFunctionEvaluatorListener listener;

    /**
     * Constructor
     *
     * @param listener listener to evaluate a single dimension function
     */
    public DerivativeEstimator(
            final SingleDimensionFunctionEvaluatorListener listener) {
        this.listener = listener;
    }

    /**
     * Computes the function derivative at provided point x.
     *
     * @param x Point where derivative is estimated
     * @return Derivative of function at provided point
     * @throws EvaluationException Raised if function cannot be properly
     *                             evaluated
     */
    public double derivative(final double x) throws EvaluationException {
        final double fold = listener.evaluate(x);

        double h = EPS * Math.abs(x);
        if (h == 0.0) {
            // Trick to reduce finite-precision error
            h = EPS;
        }
        final double xh = x + h;
        h = xh - x;

        final double fh = listener.evaluate(xh);
        return (fh - fold) / h;
    }
}
