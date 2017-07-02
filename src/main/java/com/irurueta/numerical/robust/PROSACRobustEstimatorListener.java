/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.PROSACRobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 6, 2015
 */
package com.irurueta.numerical.robust;

/**
 * Listener to get data samples and residuals for PROSAC method
 * @param <T> type of object to be estimated
 */
public interface PROSACRobustEstimatorListener<T> 
    extends RANSACRobustEstimatorListener<T> {
    
    /**
     * Returns quality scores corresponding to each sample.
     * The larger the score the better the quality
     * @return quality scores for all samples
     */
    public double[] getQualityScores();
}
