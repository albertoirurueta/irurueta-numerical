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
 * Interface to evaluate linear single dimensional functions 
 * f(x) = a * f0(x) + b * f1(x) + ...
 * Where the linear function is composed of a linear combination of a basis of
 * functions f0, f1, ... fM
 * For each evaluation at a given point x, this interface will return an array
 * containing the evaluations of the basis functions at such point f0(x), f1(x),
 * ..., fM(x)
 */
public interface LinearFitterSingleDimensionFunctionEvaluator {
    
    /**
     * Creates array where basis function results will be stored
     * @return array where basis function results will be stored
     */
    double[] createResultArray();
    
    /**
     * Evaluates a linear single dimension function at provided point and
     * returns the evaluations of the basis functions at such point
     * @param point point where function will be evaluated
     * @param result array where result of evaluation of basis functions is 
     * stored
     * @throws Throwable raised if something failed during the evaluation
     */
    void evaluate(double point, double[] result) throws Throwable;
}
