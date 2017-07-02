/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.InvalidSubsetSizeException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

/**
 * Raised if an invalid subset size is requested on a subset selector
 */
public class InvalidSubsetSizeException extends SubsetSelectorException{
    /**
     * Constructor
     */            
    public InvalidSubsetSizeException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */        
    public InvalidSubsetSizeException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */        
    public InvalidSubsetSizeException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */        
    public InvalidSubsetSizeException(Throwable cause){
        super(cause);
    }        
}
