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
 * Listener to be notified when estimation starts, finishes or any progress
 * changes.
 */
public interface PolynomialEstimatorListener {
    
    /**
     * Called when an estimator starts the polynomial estimation process.
     * @param estimator reference to a polynomial estimator.
     */
    void onEstimateStart(PolynomialEstimator estimator);
    
    /**
     * Called when an estimator ends the polynomial estimation process.
     * @param estimator reference to a polynomial estimator.
     */
    void onEstimateEnd(PolynomialEstimator estimator);
}
