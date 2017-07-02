/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.polynomials.estimators.WeightSelection
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 12, 2016.
 */
package com.irurueta.numerical.robust;

import com.irurueta.numerical.robust.WeightSelection;
import com.irurueta.sorting.SortingException;
import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class WeightSelectionTest {
    
    public static final int MIN_LENGTH = 6;
    public static final int MAX_LENGTH = 50;
    
    public static final double MIN_RANDOM_VALUE = 0.0;
    public static final double MAX_RANDOM_VALUE = 1.0;    
    
    public WeightSelectionTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }

    @Test
    public void testSelectWeights() throws SortingException{
        
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);
        
        double[] weights = new double[length];
        randomizer.fill(weights, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);
        
        //test with sorting disabled and num selected lower than length
        boolean sortWeights = false;        
        int maxPoints = length -1;
        WeightSelection selection = WeightSelection.selectWeights(weights, 
                sortWeights, maxPoints);
        assertEquals(selection.getNumSelected(), maxPoints);        
        //check first maxPoints values are true
        for(int i = 0; i < maxPoints; i++){
            assertTrue(selection.getSelected()[i]);
        }
        for(int i = maxPoints; i < length; i++){
            assertFalse(selection.getSelected()[i]);
        }
        
        //test with sorting disabled and num selected greater than length
       sortWeights = false;
       maxPoints = length + 1;
       selection = WeightSelection.selectWeights(weights, sortWeights, 
               maxPoints);
       assertEquals(selection.getNumSelected(), length);
       //check all values are true
       for(int i = 0; i < length; i++){
           assertTrue(selection.getSelected()[i]);
       }
       
       //test with sorting enabled
       sortWeights = true;
       maxPoints = length - 1;
       selection = WeightSelection.selectWeights(weights, sortWeights, 
               maxPoints);
       assertEquals(selection.getNumSelected(), maxPoints);
       //check all values are true
       int counter = 0;
       for(int i = 0; i < length; i++){
           if(selection.getSelected()[i]) counter++;
       }       
       assertEquals(counter, selection.getNumSelected());
    }
}
