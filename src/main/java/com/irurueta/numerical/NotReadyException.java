/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.NotReadyException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Raised when attempting to do a certain operation and not all parameters have
 * been provided or are correctly set.
 */
public class NotReadyException extends NumericalException{
    
    /**
     * Constructor
     */        
    public NotReadyException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */    
    public NotReadyException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */    
    public NotReadyException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */    
    public NotReadyException(Throwable cause){
        super(cause);
    }
}
