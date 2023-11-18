/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.fitting;

import com.irurueta.numerical.EvaluationException;

/**
 * Interface to evaluate linear multidimensional functions
 * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
 * Where the linear function is composed of a linear combination of a basis of
 * functions f0, f1, ... fM
 * For each evaluation at a given point (x1, x2, ...), this interface will
 * return an array containing the evaluations of the basis functions at such
 * point f0(x1, x2, ...), f1(x1, x2, ...), ..., fM(x1, x2, ...)
 */
public interface LinearFitterMultiDimensionFunctionEvaluator {

    /**
     * Number of dimensions of points (i.e. length of arrays) evaluated by
     * this function evaluator
     *
     * @return number of dimensions of points
     */
    int getNumberOfDimensions();

    /**
     * Creates array where basis function results will be stored
     *
     * @return array where basis function results will be stored
     */
    double[] createResultArray();

    /**
     * Evaluates a linear multi dimension function at provided point and
     * returns the evaluations of the basis functions at such point
     *
     * @param point  point where function will be evaluated
     * @param result array where result of evaluation of basis functions is
     *               stored
     * @throws EvaluationException raised if something failed during the evaluation
     */
    void evaluate(final double[] point, final double[] result) throws EvaluationException;
}
