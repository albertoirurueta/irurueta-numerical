/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.PolynomialsException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 8, 2016.
 */
package com.irurueta.numerical.polynomials;

import com.irurueta.numerical.NumericalException;

/**
 * Base exception for polynomials.
 */
public class PolynomialsException extends NumericalException {
    
    /**
     * Constructor.
     */
    public PolynomialsException() {
        super();
    }
    
    /**
     * Constructor with String containing message.
     * @param message message indicating the cause of the exception.
     */
    public PolynomialsException(String message) {
        super(message);
    }
    
    /**
     * Constructor with message and cause.
     * @param message message describing the cause of the exception.
     * @param cause instance containing the cause of the exception.
     */
    public PolynomialsException(String message, Throwable cause) {
        super(message, cause);
    }
    
    /**
     * Constructor with cause.
     * @param cause instance containing the cause of the exception.
     */
    public PolynomialsException(Throwable cause) {
        super(cause);
    }
}
