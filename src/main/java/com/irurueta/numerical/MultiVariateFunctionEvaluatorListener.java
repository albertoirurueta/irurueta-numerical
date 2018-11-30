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
package com.irurueta.numerical;

/**
 * Interface to define how multivariate functions can be evaluated.
 */
public interface MultiVariateFunctionEvaluatorListener {
    /**
     * Evaluates a multi variate function such as f1(x1, x2, ...), 
     * f2(x1, x2, ...) at providd multidimensional point and returns the result
     * as a vectorial value
     * @param point multidimensional point where function will be evaluated
     * @param result vector where function evaluation will be stored
     * @throws Throwable if something failed during the evaluation
     */
    void evaluate(double[] point, double[] result) throws Exception;
    
    /**
     * Number of variables of function f. This is equal to the length of the
     * array obtained as function evaluations. Hence, a function f can
     * be expressed as f = [f1, f2, ... fN], and the number of variables would
     * be N
     * @return number of variables of function f
     */
    int getNumberOfVariables();
}
