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
 * Multi-quadric Radial Function Basis implementation.
 */
public class MultiQuadricRadialBasisFunction implements RadialBasisFunction {

    /**
     * Squared r0 value, which is a scale factor.
     */
    private final double r02;

    /**
     * Constructor.
     *
     * @param scale scale factor.
     */
    public MultiQuadricRadialBasisFunction(final double scale) {
        r02 = scale * scale;
    }

    /**
     * Constructor.
     * Uses default scale factor, which is 1.0.
     */
    public MultiQuadricRadialBasisFunction() {
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
        return Math.sqrt((r * r) + r02);
    }
}
