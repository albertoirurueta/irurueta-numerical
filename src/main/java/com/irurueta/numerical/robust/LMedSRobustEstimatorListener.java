/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.LMedSRobustEstimatorListener
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date February 4, 2015
 */
package com.irurueta.numerical.robust;

import java.util.List;

/**
 * Listener to get data samples and residuals for LMedS method
 * @param <T> type of object to be estimated
 */
public interface LMedSRobustEstimatorListener<T> 
    extends RobustEstimatorListener<T>{
    /**
     * Returns total number of samples to be randomly processed
     * @return total number of samples
     */
    public int getTotalSamples();
    
    /**
     * Returns size of subsets of samples to be selected
     * @return size of subsets of samples to be selected
     */
    public int getSubsetSize();  
    
    /**
     * Estimates a list of possible preliminar solutions to be used during an
     * iteration of LMedS robust estimator
     * @param samplesIndices indices of random subset of samples that have been
     * picked
     * @param solutions list where possible preliminar solutions to be used 
     * during an iteration of LMedS robust estimator will be stored. Provided 
     * list will always be empty, and it is up to the implementor to fill it
     * with preliminar solutions based on provided sample indices
     */
    public void estimatePreliminarSolutions(
            int[] samplesIndices, List<T> solutions);    
    
    /**
     * Computes residual for sample located at i-th position using estimation
     * on current iteration.
     * @param currentEstimation a preliminar estimation that has been found for
     * current iteration
     * @param i position of sample to be checked
     * @return residual for i-th sample
     */
    public double computeResidual(T currentEstimation, int i);    
}
