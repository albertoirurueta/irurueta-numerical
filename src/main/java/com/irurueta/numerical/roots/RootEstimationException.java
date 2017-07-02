/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.roots.RootEstimationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 11, 2012
 */
package com.irurueta.numerical.roots;

import com.irurueta.numerical.NumericalException;

/**
 * Raised when a root estimator cannot determine a root of a polynomial, usually
 * because of lack of convergence
 */
public class RootEstimationException extends NumericalException{
    
    /**
     * Constructor
     */            
    public RootEstimationException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */        
    public RootEstimationException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */        
    public RootEstimationException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */        
    public RootEstimationException(Throwable cause){
        super(cause);
    }    
}
