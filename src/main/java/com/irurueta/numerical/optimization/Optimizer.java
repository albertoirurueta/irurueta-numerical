/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.optimization.Optimizer
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical.optimization;

import com.irurueta.numerical.LockedException;
import com.irurueta.numerical.NotReadyException;

/**
 * Abstract class to find function minima.
 * Implementations will take into account whether the function is single or
 * multidimension, and will use different algorithms to find minima.
 */
public abstract class Optimizer {
    
    /**
     * Boolean indicating whether this instance is locked because computations
     * are being done.
     */
    protected boolean locked;
    
    /**
     * Empty constructor.
     */
    public Optimizer() {
        locked = false;
    }
    
    /**
     * Returns boolean indicating whether this instance is locked.
     * This instance will be locked while computations are being done.
     * Attempting to change any parameter while this instance is locked will
     * raise a LockedException.
     * @return True if this instance is locked, false otherwise.
     */
    public boolean isLocked() {
        return locked;
    }
    
    /**
     * This function estimates a function minimum.
     * Implementations of this class will usually search a local minimum within
     * a bracket of input values.
     * Because this is an abstract class, this method is meant to be overridden,
     * otherwise a NotReadyException will always be thrown.
     * @throws LockedException Raised if this instance is locked, because
     * estimation is being computed.
     * @throws NotReadyException Raised if this instance is not ready, usually
     * because listener has not yet been provided.
     * @throws OptimizationException Raised if the algorithm failed because of
     * lack of convergence or because function couldn't be evaluated.
     */
    public void minimize() throws LockedException, NotReadyException, 
            OptimizationException {
        throw new NotReadyException();
    }
    
    /**
     * Returns boolean indicating whether this instance is ready.
     * Usually an instance will be ready once its listener has been provided.
     * Because this is an abstract class, it will always return false;
     * @return True if this instance is ready, false otherwise.
     */
    public boolean isReady() {
        return false;
    }
}
