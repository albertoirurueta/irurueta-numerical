/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.EvaluationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Exception raised when function evaluation fails
 */
public class EvaluationException extends Exception{
    
    /**
     * Constructor
     */            
    public EvaluationException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */        
    public EvaluationException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */        
    public EvaluationException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */        
    public EvaluationException(Throwable cause){
        super(cause);
    }    
}
