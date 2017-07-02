/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.robust.SubsetSelectorException
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

public class SubsetSelectorExceptionTest {
    
    public SubsetSelectorExceptionTest() {}
    
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
        SubsetSelectorException ex;
        assertNotNull(ex = new SubsetSelectorException());
        
        ex = null;
        assertNotNull(ex = new SubsetSelectorException("message"));
        
        ex = null;
        assertNotNull(ex = new SubsetSelectorException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new SubsetSelectorException("message", 
                new Exception()));        
    }                  
}