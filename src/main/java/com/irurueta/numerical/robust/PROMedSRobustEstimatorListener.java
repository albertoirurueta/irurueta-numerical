/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.PROMedSRobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 7, 2015
 */
package com.irurueta.numerical.robust;

/**
 * Listener to get data samples and residuals for PROMedS method
 * @param <T> type of object to be estimated
 */
public interface PROMedSRobustEstimatorListener<T> 
    extends PROSACRobustEstimatorListener<T> {   
}
