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
package com.irurueta.numerical.robust;

import java.util.List;

/**
 * Listener to get data samples and residuals for LMedS method.
 *
 * @param <T> type of object to be estimated.
 */
public interface LMedSRobustEstimatorListener<T>
        extends RobustEstimatorListener<T> {
    /**
     * Returns total number of samples to be randomly processed.
     *
     * @return total number of samples.
     */
    int getTotalSamples();

    /**
     * Returns size of subsets of samples to be selected.
     *
     * @return size of subsets of samples to be selected.
     */
    int getSubsetSize();

    /**
     * Estimates a list of possible preliminar solutions to be used during an
     * iteration of LMedS robust estimator.
     *
     * @param samplesIndices indices of random subset of samples that have been
     *                       picked.
     * @param solutions      list where possible preliminar solutions to be used
     *                       during an iteration of LMedS robust estimator will be stored. Provided
     *                       list will always be empty, and it is up to the implementor to fill it
     *                       with preliminar solutions based on provided sample indices.
     */
    void estimatePreliminarSolutions(
            final int[] samplesIndices, final List<T> solutions);

    /**
     * Computes residual for sample located at i-th position using estimation
     * on current iteration.
     *
     * @param currentEstimation a preliminar estimation that has been found for
     *                          current iteration.
     * @param i                 position of sample to be checked.
     * @return residual for i-th sample.
     */
    double computeResidual(final T currentEstimation, final int i);
}
