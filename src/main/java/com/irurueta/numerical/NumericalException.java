/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.NumericalException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Base class for all the exceptions in this package
 */
public class NumericalException extends Exception{
    
    /**
     * Constructor
     */    
    public NumericalException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */    
    public NumericalException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */    
    public NumericalException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */    
    public NumericalException(Throwable cause){
        super(cause);
    }
}
