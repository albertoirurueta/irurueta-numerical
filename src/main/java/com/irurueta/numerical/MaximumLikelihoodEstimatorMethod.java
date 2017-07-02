/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.MaximumLikelihoodEstimatorMethod
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 5, 2012
 */
package com.irurueta.numerical;

/**
 * Types of maximum likelihood estimation to determine the real value 
 * corresponing to a set of values
 */
public enum MaximumLikelihoodEstimatorMethod {
    /**
     * MLE method based on a histogram of all samples
     */
    HISTOGRAM_MAXIMUM_LIKELIHOOD_ESTIMATOR,
    
    /**
     * MLE method that refines the histogram method by using Gaussian 
     * interpolation
     */
    ACCURATE_MAXIMUM_LIKELIHOOD_ESTIMATOR
}
