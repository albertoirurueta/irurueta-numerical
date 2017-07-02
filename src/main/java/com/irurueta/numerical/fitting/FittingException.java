/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fitting.FittingException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 21, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.numerical.NumericalException;

/**
 * Raised when a fitter fails to fit a function to provided data
 */
public class FittingException extends NumericalException{
    
    /**
     * Constructor
     */
    public FittingException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message message indicating the cause of the exception
     */
    public FittingException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message message describing the cause of the exception
     * @param cause instance containing the cause of the exception
     */
    public FittingException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause instance containing the cause of the exception
     */
    public FittingException(Throwable cause){
        super(cause);
    }
}
