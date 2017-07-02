/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.NotEnoughSamplesException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

/**
 * Raised if there aren't enough samples to make a computation
 */
public class NotEnoughSamplesException  extends SubsetSelectorException{
    
    /**
     * Constructor
     */        
    public NotEnoughSamplesException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */    
    public NotEnoughSamplesException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */    
    public NotEnoughSamplesException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */    
    public NotEnoughSamplesException(Throwable cause){
        super(cause);
    }
}
