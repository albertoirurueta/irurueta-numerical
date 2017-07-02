/**
 * @file
 * This file contains Unit Tests for
 * com.irurueta.numerical.InvalidBracketRangeException
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 23, 2012
 */
package com.irurueta.numerical;

import org.junit.*;
import static org.junit.Assert.*;

public class InvalidBracketRangeExceptionTest {
    
    public InvalidBracketRangeExceptionTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testConstructor(){
        InvalidBracketRangeException ex;
        assertNotNull(ex = new InvalidBracketRangeException());
        
        ex = null;
        assertNotNull(ex = new InvalidBracketRangeException("message"));
        
        ex = null;
        assertNotNull(ex = new InvalidBracketRangeException(new Exception()));
        
        ex = null;
        assertNotNull(ex = new InvalidBracketRangeException("message", 
                new Exception()));        
    }      
}
