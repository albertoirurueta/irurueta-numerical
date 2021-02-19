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

/**
 * Listener to get data samples and residuals for PROSAC method
 *
 * @param <T> type of object to be estimated.
 */
public interface PROSACRobustEstimatorListener<T>
        extends RANSACRobustEstimatorListener<T> {

    /**
     * Returns quality scores corresponding to each sample.
     * The larger the score the better the quality.
     *
     * @return quality scores for all samples.
     */
    double[] getQualityScores();
}
