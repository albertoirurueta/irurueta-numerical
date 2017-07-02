/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.signal.processing.SignalProcessingException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date October 11, 2015
 */
package com.irurueta.numerical.signal.processing;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class SignalProcessingExceptionTest {
    
    public SignalProcessingExceptionTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}
    
    @Test
    public void testConstructor(){
        SignalProcessingException ex;
        assertNotNull(ex = new SignalProcessingException());
        
        ex = null;
        assertNotNull(ex = new SignalProcessingException("message"));
        
        ex = null;
        assertNotNull(ex = new SignalProcessingException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new SignalProcessingException("message",
                new Exception()));
    }
}
