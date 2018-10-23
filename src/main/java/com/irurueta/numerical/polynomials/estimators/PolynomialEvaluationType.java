/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.polynomials.estimators;

/**
 * Determines different types of polynomial evaluations that can be used
 * to estimate a polynomial.
 */
public enum PolynomialEvaluationType {
    /**
     * A direct evaluation of a polynomial.
     */
    DIRECT_EVALUATION,
    
    /**
     * Evaluation of the nth-derivative of a polynomial.
     */
    DERIVATIVE_EVALUATION,
    
    /**
     * Evaluation of the nth-integral of a polynomial (assuming a given 
     * constant value)
     */
    INTEGRAL_EVALUATION,
    
    /**
     * Interval nth-integration of a polynomial.
     */
    INTEGRAL_INTERVAL
}
