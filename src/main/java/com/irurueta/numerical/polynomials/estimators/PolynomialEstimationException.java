/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialEstimationException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 8, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

import com.irurueta.numerical.polynomials.PolynomialsException;

/**
 * Exception raised if polynomial estimation fails.
 */
public class PolynomialEstimationException extends PolynomialsException {
    
    /**
     * Constructor.
     */
    public PolynomialEstimationException() {
        super();
    }
    
    /**
     * Constructor with String containing message.
     * @param message message indicating the cause of the exception.
     */
    public PolynomialEstimationException(String message) {
        super(message);
    }
    
    /**
     * Constructor with message and cause.
     * @param message message describing the cause of the exception.
     * @param cause instance containing the cause of the exception.
     */
    public PolynomialEstimationException(String message, Throwable cause) {
        super(message, cause);
    }
    
    /**
     * Constructor with cause.
     * @param cause instance containing the cause of the exception.
     */
    public PolynomialEstimationException(Throwable cause) {
        super(cause);
    }
}
