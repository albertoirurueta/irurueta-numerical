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
 * Listener to be notified of events such as when estimation starts, ends or
 * when progress changes.
 */
public interface PolynomialRobustEstimatorListener {

    /**
     * Called when an estimator starts the polynomial estimation process.
     *
     * @param estimator reference to a polynomial robust estimator.
     */
    void onEstimateStart(final PolynomialRobustEstimator estimator);

    /**
     * Called when an estimator ends the polynomial estimation process.
     *
     * @param estimator reference to a polynomial robust estimator.
     */
    void onEstimateEnd(final PolynomialRobustEstimator estimator);

    /**
     * Called when estimator iterates to refine a possible solution.
     *
     * @param estimator reference to a polynomial robust estimator.
     * @param iteration current iteration.
     */
    void onEstimateNextIteration(final PolynomialRobustEstimator estimator, final int iteration);

    /**
     * Called when estimation progress changes significantly.
     *
     * @param estimator reference to a polynomial robust estimator.
     * @param progress  progress of estimation expressed as a value between 0.0
     *                  and 1.0.
     */
    void onEstimateProgressChange(final PolynomialRobustEstimator estimator, final float progress);
}
