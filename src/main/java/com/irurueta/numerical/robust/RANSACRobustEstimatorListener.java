/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.RANSACRobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date April 9, 2013
 */
package com.irurueta.numerical.robust;

/**
 *
 * Listener to get data samples and residuals for RANSAC method
 * @param <T> type of object to be estimated
 */
public interface RANSACRobustEstimatorListener<T> 
    extends LMedSRobustEstimatorListener<T>{
        
    /**
     * Threshold to determine whether samples are inliers or not.
     * Threshold must be a positive value, otherwise RANSAC algorithm will fail
     * @return threshold to determine whether samples are inliers or not.
     */
    public double getThreshold();
}
