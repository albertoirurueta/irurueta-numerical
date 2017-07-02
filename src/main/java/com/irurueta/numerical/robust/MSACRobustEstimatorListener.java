/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.MSACRobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 6, 2015
 */
package com.irurueta.numerical.robust;

/**
 * Listener to get data samples and residuals for MSAC method
 * @param <T> type of object to be estimated
 */
public interface MSACRobustEstimatorListener<T> extends 
        RANSACRobustEstimatorListener<T> {    }
