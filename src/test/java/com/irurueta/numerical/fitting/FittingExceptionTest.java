/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.fitting.FittingException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 22, 2015
 */
package com.irurueta.numerical.fitting;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class FittingExceptionTest {
    
    public FittingExceptionTest() {}
    
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
        FittingException ex;
        assertNotNull(ex = new FittingException());
        
        ex = null;
        assertNotNull(ex = new FittingException("message"));
        
        ex = null;
        assertNotNull(ex = new FittingException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new FittingException("message", new Exception()));
    }
}
