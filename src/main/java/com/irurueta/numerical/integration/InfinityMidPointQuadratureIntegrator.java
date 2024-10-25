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
 * Computes function integration by using an infinity mid-point quadrature.
 * If you need to integrate from a negative lower limit to positive infinity, you do this by
 * breaking the integral into two pieces at some positive value:
 * InfinityMidPointQuadratureIntegrator integrator1 =
 * new InfinityMidPointQuadratureIntegrator(-5.0, 2.0, listener);
 * InfinityMidPointQuadratureIntegrator integrator2 =
 * new InfinityMidPointQuadratureIntegrator(2.0, 1e99, listener);
 * double integral = integrator1.integrate() + integrator2.integrate();
 */
public class InfinityMidPointQuadratureIntegrator extends QuadratureIntegrator<InfinityMidPointQuadrature> {

    /**
     * Constructor.
     *
     * @param a        Lower limit of integration.
     * @param b        Upper limit of integration.
     * @param listener listener to evaluate a single dimension function at required points.
     * @param eps      required accuracy.
     */
    public InfinityMidPointQuadratureIntegrator(
            final double a, final double b, final SingleDimensionFunctionEvaluatorListener listener, final double eps) {
        super(new InfinityMidPointQuadrature(a, b, listener), eps);
    }

    /**
     * Constructor.
     *
     * @param a        Lower limit of integration.
     * @param b        Upper limit of integration.
     * @param listener listener to evaluate a single dimension function at required points.
     */
    public InfinityMidPointQuadratureIntegrator(
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
        return QuadratureType.INFINITY_MID_POINT;
    }
}
