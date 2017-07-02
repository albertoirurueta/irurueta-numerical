/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.signal.processing.SignalProcessingException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date October 11, 2015
 */
package com.irurueta.numerical.signal.processing;

import com.irurueta.numerical.NumericalException;

/**
 * Raised when something fails during signal processing
 */
public class SignalProcessingException extends NumericalException{
    
    /**
     * Constructor.
     */
    public SignalProcessingException(){
        super();
    }
    
    /**
     * Constructor with String containing message.
     * @param message message indicating the cause of the exception
     */
    public SignalProcessingException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause.
     * @param message message describing the cause of the exception
     * @param cause instance containing the cause of the exception
     */
    public SignalProcessingException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause.
     * @param cause instance containing the cause of the exception
     */
    public SignalProcessingException(Throwable cause){
        super(cause);
    }
}
