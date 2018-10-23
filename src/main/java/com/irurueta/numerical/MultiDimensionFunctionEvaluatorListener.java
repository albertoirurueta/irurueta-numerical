/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
 * Interface to define how multi dimension functions can be evaluated.
 * This interface is used in several algorithms to provide methods to evaluate
 * functions and retrieve their minima/maxima, etc.
 */
public interface MultiDimensionFunctionEvaluatorListener {
    
    /**
     * Evaluates a multi dimension function such as f([x1, x2, ..., xn]) at 
     * provided multidimensional point and returns the result as a scalar value.
     * @param point Multidimensional point where function will be evaluated. 
     * This must be an array of length equal to the dimensionality of the 
     * function.
     * @return Value returned by the function.
     * @throws Throwable Raised if something failed during the evaluation.
     */    
    double evaluate(double[] point) throws Throwable;
}
