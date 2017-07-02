/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.robust.InvalidSubsetRangeException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 10, 2013
 */
package com.irurueta.numerical.robust;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class InvalidSubsetRangeExceptionTest {
    
    public InvalidSubsetRangeExceptionTest() {}
    
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
        InvalidSubsetRangeException ex;
        assertNotNull(ex = new InvalidSubsetRangeException());
        
        ex = null;
        assertNotNull(ex = new InvalidSubsetRangeException("message"));
        
        ex = null;
        assertNotNull(ex = new InvalidSubsetRangeException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new InvalidSubsetRangeException("message", 
                new Exception()));        
    }                          
}