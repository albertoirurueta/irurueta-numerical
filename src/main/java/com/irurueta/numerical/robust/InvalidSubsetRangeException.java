/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.InvalidSubsetRangeException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

/**
 * Raised if provided range of samples to pick subsets from is invalid
 */
public class InvalidSubsetRangeException extends SubsetSelectorException{
    /**
     * Constructor
     */            
    public InvalidSubsetRangeException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */        
    public InvalidSubsetRangeException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */        
    public InvalidSubsetRangeException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */        
    public InvalidSubsetRangeException(Throwable cause){
        super(cause);
    }            
}
