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
package com.irurueta.numerical.interpolation;

/**
 * Gaussian Radial Basis Function implementation.
 */
public class GaussianRadialBasisFunction implements RadialBasisFunction {

    /**
     * Scale factor.
     */
    private final double r0;

    /**
     * Constructor.
     *
     * @param scale scale factor.
     */
    public GaussianRadialBasisFunction(final double scale) {
        r0 = scale;
    }

    /**
     * Constructor.
     * Uses default scale factor, which is 1.0.
     */
    public GaussianRadialBasisFunction() {
        this(1.0);
    }

    /**
     * Evaluates RBF at provided distance between two points.
     *
     * @param r distance between two points.
     * @return result of evaluating RBF.
     */
    @Override
    public double evaluate(double r) {
        final var value = r / r0;
        return Math.exp(-0.5 * value * value);
    }
}
