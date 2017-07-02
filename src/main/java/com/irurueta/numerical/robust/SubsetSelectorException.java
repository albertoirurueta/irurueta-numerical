/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.SubsetSelectorException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

import com.irurueta.numerical.NumericalException;

/**
 * Raised if subset selection of samples fails
 */
public class SubsetSelectorException extends NumericalException{
    
    /**
     * Constructor
     */        
    public SubsetSelectorException(){
        super();
    }
    
    /**
     * Constructor with String containing message
     * @param message Message indicating the cause of the exception
     */    
    public SubsetSelectorException(String message){
        super(message);
    }
    
    /**
     * Constructor with message and cause
     * @param message Message describing the cause of the exception
     * @param cause Instance containing the cause of the exception
     */    
    public SubsetSelectorException(String message, Throwable cause){
        super(message, cause);
    }
    
    /**
     * Constructor with cause
     * @param cause Instance containing the cause of the exception
     */    
    public SubsetSelectorException(Throwable cause){
        super(cause);
    }
}
