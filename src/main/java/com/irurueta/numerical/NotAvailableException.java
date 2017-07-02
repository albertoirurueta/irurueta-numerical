/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.NotAvailableException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 1, 2012
 */
package com.irurueta.numerical;

/**
 * Exception raised when some value cannot be retrieved, usually because it has
 * not yet been provided or computed
 */
public class NotAvailableException extends NumericalException{
    
    /**
     * Constructor
     */    
    public NotAvailableException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */    
    public NotAvailableException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */    
    public NotAvailableException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */    
    public NotAvailableException(Throwable cause){
        super(cause);
    }    
}
