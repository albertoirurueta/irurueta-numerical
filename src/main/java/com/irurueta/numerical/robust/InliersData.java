/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.InliersData
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date December 26, 2016.
 */
package com.irurueta.numerical.robust;

import java.util.BitSet;

/**
 * Base clase defining inlier data for a robust estimator.
 */
public abstract class InliersData {
    
    /**
     * Residuals obtained for each sample of data.
     */
    protected double[] mResiduals;
    
    /**
     * Number of inliers found on current iteration.
     */
    protected int mNumInliers;
    
    /**
     * Returns efficient array indicating which samples are considered inliers
     * and which ones aren't or null if inliers are not kept.
     * @return array indicating which samples are considered inliers and which 
     * ones aren't, or null.
     */
    public abstract BitSet getInliers();

    /**
     * Returns residuals obtained for each sample of data or null if residuals
     * are not kept.
     * @return residuals obtained for each sample of data.
     */
    public double[] getResiduals() {
        return mResiduals;
    }
    
    /**
     * Returns number of inliers found.
     * @return number of inliers found.
     */
    public int getNumInliers() {
        return mNumInliers;
    }
}
