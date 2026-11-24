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

import com.irurueta.numerical.SingleDimensionFunctionEvaluatorListener;

/**
 * Computes function integration by using Simpson's rule an Upper Square Root mid-point quadrature
 * when upper integration bound lies at a function singularity.
 * <p>
 * Because mid-point quadratures refine the step by a factor of 3 on each iteration, successive
 * estimates are combined using a 9/8 extrapolation (instead of the 4/3 used with trapezoidal quadratures),
 * and the maximum number of refinements is reduced accordingly.
 * See Numerical Recipes, 3rd ed., section 4.4, p. 169.
 */
public class SimpsonUpperSquareRootMidPointQuadratureIntegrator
        extends SimpsonIntegrator<UpperSquareRootMidPointQuadrature> {

    /**
     * Constructor.
     *
     * @param a        Lower limit of integration.
     * @param b        Upper limit of integration.
     * @param listener listener to evaluate a single dimension function at required points.
     * @param eps      required accuracy.
     */
    public SimpsonUpperSquareRootMidPointQuadratureIntegrator(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener, final double eps) {
        super(new UpperSquareRootMidPointQuadrature(a, b, listener), eps, TRIPLING_RATIO, TRIPLING_JMAX);
    }

    /**
     * Constructor.
     *
     * @param a        Lower limit of integration.
     * @param b        Upper limit of integration.
     * @param listener listener to evaluate a single dimension function at required points.
     */
    public SimpsonUpperSquareRootMidPointQuadratureIntegrator(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener) {
        this(a, b, listener, EPS);
    }

    /**
     * Gets type of quadrature.
     *
     * @return type of quadrature.
     */
    @Override
    public QuadratureType getQuadratureType() {
        return QuadratureType.UPPER_SQUARE_ROOT_MID_POINT;
    }
}
