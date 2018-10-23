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
 * Listener to evaluate/retrieve a multidimensional function's gradient.
 */
public interface GradientFunctionEvaluatorListener {
    
    /**
     * Computes/retrieves a multidimensional function's gradient.
     * @param params Array of input parameters of a multidimensional function
     * @param result Array containing estimated gradient. This parameter must 
     * contain an already instantiated array having the same length as params.
     * The values of this gradient will be rewritten after executing this method
     * and this array will contain the estimated or evaluated gradient
     * @throws Throwable Raised if something fails.
     */
    void evaluateGradient(double[] params, double[] result)
            throws Throwable;
}
