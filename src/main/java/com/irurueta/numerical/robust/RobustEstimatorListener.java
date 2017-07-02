/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.RobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date March 4, 2013
 */
package com.irurueta.numerical.robust;

/**
 *
 * Listener to be notified of events on a robust estimator such as when 
 * estimation starts, ends or when progress changes
 * @param <T> type of instance being estimated
 */
public interface RobustEstimatorListener<T> {
    
    /**
     * Called to determine when a robust estimator is ready to start estimation.
     * This is true when all required data to start estimation is available
     * @return true if robust estimator is ready to start estimation, false
     * otherwise
     */
    public boolean isReady();
   
    /**
     * Called when estimation starts
     * @param estimator reference to robust estimator
     */
    public void onEstimateStart(RobustEstimator<T> estimator);
    
    /**
     * Called when estimation ends
     * @param estimator reference to robust estimator
     */
    public void onEstimateEnd(RobustEstimator<T> estimator);
    
    /**
     * Called when estimator iterates to refine a possible solution
     * @param estimator reference to robust estimator
     * @param iteration current iteration
     */
    public void onEstimateNextIteration(RobustEstimator<T> estimator, 
            int iteration);
    
    /**
     * Called when estimation progress changes significantly.
     * @param estimator reference to robust estimator
     * @param progress progress of estimation expressed as a value between 0.0
     * and 1.0
     */
    public void onEstimateProgressChange(RobustEstimator<T> estimator, 
            float progress);
}
