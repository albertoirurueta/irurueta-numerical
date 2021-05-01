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

import com.irurueta.algebra.Matrix;
import com.irurueta.numerical.EvaluationException;

/**
 * Interface to evaluate non-linear multi variate and multi dimensional
 * functions.
 * Evaluation of functions requires both function value at provided point x and
 * function jacobian respect to its parameters (i.e. derivatives respect to its
 * parameters for each function output or variable)
 */
public interface LevenbergMarquardtMultiVariateFunctionEvaluator {

    /**
     * Number of dimensions of points (i.e. length of input points arrays)
     * evaluated by this function evaluator
     *
     * @return number of dimensions of points
     */
    int getNumberOfDimensions();

    /**
     * Number of variables of function f. This is equal to the length of the
     * array obtained as function evaluations. Hence, a function f can
     * be expressed as f = [f1, f2, ... fN], and the number of variables would
     * be N
     *
     * @return number of variables of function f
     */
    int getNumberOfVariables();

    /**
     * Creates array where estimated parameters will be stored.
     * This array MUST contain the initial guessed solution for the Levenberg-
     * Marquardt algorithm
     *
     * @return array where estimated parameters will be stored
     */
    double[] createInitialParametersArray();

    /**
     * Evaluates a non-linear multi variate function at provided point using
     * provided parameters and returns its evaluation and jacobian of the
     * function respect the function parameters
     *
     * @param i        number of sample being evaluated
     * @param point    point where function will be evaluated
     * @param result   result of function evaluation. Its length is equal to the
     *                 number of function variables
     * @param params   initial parameters estimation to be tried. These will
     *                 change as the Levenberg-Marquard algorithm iterates to the best solution.
     *                 These are used as input parameters along with point to evaluate function
     * @param jacobian jacobian containing partial derivatives of the function
     *                 respect to each provided parameter for each function output or variable
     * @throws EvaluationException raised if something failed during the evaluation
     */
    void evaluate(final int i, final double[] point, final double[] result,
                  final double[] params, final Matrix jacobian) throws EvaluationException;

}
