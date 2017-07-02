/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fitting.Fitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 21, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.numerical.NotReadyException;

/**
 * Base class for function fitters used to estimate function parameters along
 * with their covariance matrix and chi square value
 */
public abstract class Fitter {
    
    /**
     * Indicates whether result has been estimated and is available for 
     * retrieval
     */
    protected boolean resultAvailable;
    
    /**
     * Returns boolean indicating whether result has been estimated and is 
     * available for retrieval
     * @return true if result has been estimated and is available for retrieval
     */
    public boolean isResultAvailable(){
        return resultAvailable;
    }
    
    /**
     * Indicates whether this instance is ready because enough input data has 
     * been provided to start the fitting process
     * @return true if this fitter is ready, false otherwise
     */
    public abstract boolean isReady();
    
    /**
     * Fits a function to provided data so that parameters associated to that
     * function can be estimated along with their covariance matrix and chi
     * square value
     * @throws FittingException if fitting fails
     * @throws NotReadyException if enough input data has not yet been provided
     */
    public abstract void fit() throws FittingException, NotReadyException;
}
