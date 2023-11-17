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

import com.irurueta.numerical.EvaluationException;

/**
 * Interface to define how single dimension functions can be evaluated in Double Exponential Rule
 * Quadrature function integrators.
 */
public interface DoubleExponentialSingleDimensionFunctionEvaluatorListener {

    /**
     * Evaluates a single dimension function such as f(x) at provided point and returns the result.
     *
     * @param x     point where function will be evaluated.
     * @param delta value to handle singularities. If the function has no singularities, or the
     *              singularities are “mild” (e.g., no worse than logarithmic), this can be ignored.
     * @return value returned by the function.
     * @throws EvaluationException raised if something failed during the evaluation.
     */
    double evaluate(final double x, final double delta) throws EvaluationException;
}
