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

/**
 * Interface to evaluate non-linear multi dimensional functions.
 * Evaluation of functions requires both function value at provided point x and
 * function gradient respect to its parameters (i.e. derivatives respect to its
 * parameters).
 */
public interface LevenbergMarquardtMultiDimensionFunctionEvaluator {
    
    /**
     * Number of dimensions of points (i.e. length of arrays) evaluated by
     * this function evaluator.
     * @return number of dimensions of points.
     */
    int getNumberOfDimensions();
    
    /**
     * Creates array where estimated parameters will be stored.
     * This array MUST contain the initial guessed solution for the LEvenberg-
     * Marquardt algorithm.
     * @return array where estimated parameters will be stored.
     */
    double[] createInitialParametersArray();
    
    /**
     * Evaluates a non-linear multi dimension function at provided point using 
     * provided parameters and returns its evaluation and derivatives of the 
     * function respect the function parameters.
     * @param i number of sample being evaluated.
     * @param point point where function will be evaluated.
     * @param params initial parameters estimation to be tried. These will 
     * change as the Levenberg-Marquard algorithm iterates to the best solution.
     * These are used as input parameters along with point to evaluate function.
     * @param derivatives partial derivatives of the function respect to each
     * provided parameter.
     * @return function evaluation at provided point.
     * @throws Throwable raised if something failed during the evaluation.
     */
    double evaluate(int i, double[] point, double[] params,
            double[] derivatives) throws Throwable;    
}
