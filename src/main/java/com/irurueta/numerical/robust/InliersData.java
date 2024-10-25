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
package com.irurueta.numerical.robust;

import java.util.BitSet;

/**
 * Base class defining inlier data for a robust estimator.
 */
public abstract class InliersData {

    /**
     * Residuals obtained for each sample of data.
     */
    protected double[] residuals;

    /**
     * Number of inliers found on current iteration.
     */
    protected int numInliers;

    /**
     * Returns efficient array indicating which samples are considered inliers
     * and which ones aren't or null if inliers are not kept.
     *
     * @return array indicating which samples are considered inliers and which
     * ones aren't, or null.
     */
    public abstract BitSet getInliers();

    /**
     * Returns residuals obtained for each sample of data or null if residuals
     * are not kept.
     *
     * @return residuals obtained for each sample of data.
     */
    public double[] getResiduals() {
        return residuals;
    }

    /**
     * Returns number of inliers found.
     *
     * @return number of inliers found.
     */
    public int getNumInliers() {
        return numInliers;
    }
}
