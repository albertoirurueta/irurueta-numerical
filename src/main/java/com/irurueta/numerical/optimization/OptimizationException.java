/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.optimization.OptimizationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical.optimization;

import com.irurueta.numerical.NumericalException;

/**
 * Raised when an optimizer cannot find a minimum on a function, usually because
 * of lack of convergence.
 */
public class OptimizationException extends NumericalException {
    
    /**
     * Constructor.
     */            
    public OptimizationException() {
        super();
    }
    
    /**
     * Constructor with String containing message.
     * @param message Message indicating the cause of the exception.
     */        
    public OptimizationException(String message) {
        super(message);
    }
    
    /**
     * Constructor with message and cause.
     * @param message Message describing the cause of the exception.
     * @param cause Instance containing the cause of the exception.
     */        
    public OptimizationException(String message, Throwable cause) {
        super(message, cause);
    }
    
    /**
     * Constructor with cause.
     * @param cause Instance containing the cause of the exception.
     */        
    public OptimizationException(Throwable cause) {
        super(cause);
    }
}
