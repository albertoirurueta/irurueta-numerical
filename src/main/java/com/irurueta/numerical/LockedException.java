/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.LockedException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Exception raised when an instance is locked
 */
public class LockedException extends NumericalException{
    
    /**
     * Constructor
     */    
    public LockedException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */
    public LockedException(String message){
        super(message);
    }

    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */
    public LockedException(String message, Throwable cause){
        super(message, cause);
    }

    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */
    public LockedException(Throwable cause){
        super(cause);
    }                
}
