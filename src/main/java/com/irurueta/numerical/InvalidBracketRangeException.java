/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.InvalidBracketRangeException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Exception raised when provided bracket of values is not valid
 */
public class InvalidBracketRangeException extends NumericalException{
    /**
     * Constructor
     */    
    public InvalidBracketRangeException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */    
    public InvalidBracketRangeException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */    
    public InvalidBracketRangeException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */    
    public InvalidBracketRangeException(Throwable cause){
        super(cause);
    }    
}
