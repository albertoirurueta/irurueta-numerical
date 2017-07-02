/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.robust.SubsetSelector
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date August 11, 2013
 */
package com.irurueta.numerical.robust;

import com.irurueta.statistics.UniformRandomizer;
import java.util.Random;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class SubsetSelectorTest {
    
    public static final int MIN_SAMPLES = 100;
    public static final int MAX_SAMPLES = 500;
    
    public static final int MIN_SUBSET_SIZE = 5;
    public static final int MAX_SUBSET_SIZE = 10;    
    
    public SubsetSelectorTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}
    
    @Test
    public void testCreate(){
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        
        //Test create with number of samples and type
        SubsetSelector selector = SubsetSelector.create(numSamples, 
                SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);
        assertEquals(selector.getType(), 
                SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);
        assertEquals(selector.getNumSamples(), numSamples);
        
        //Force IllegalArgumentException
        selector = null;
        try{
            selector = SubsetSelector.create(0, 
                    SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(selector);
        
        //Test create with number of samples
        selector = SubsetSelector.create(numSamples);
        assertEquals(selector.getType(), 
                SubsetSelector.DEFAULT_SUBSET_SELECTOR_TYPE);
        assertEquals(selector.getNumSamples(), numSamples);
        
        //Force IllegalArgumentException
        selector = null;
        try{
            selector = SubsetSelector.create(0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
        assertNull(selector);
    }
    
    @Test
    public void testGetType(){
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        
        SubsetSelector selector = SubsetSelector.create(numSamples);
        assertEquals(selector.getType(), 
                SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);        
    }
    
    @Test
    public void testGetSetNumSamples(){
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples1 = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        int numSamples2 = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        
        SubsetSelector selector = SubsetSelector.create(numSamples1);
        //check correctness
        assertEquals(selector.getNumSamples(), numSamples1);
        
        //set new value
        selector.setNumSamples(numSamples2);
        
        //check correctness
        assertEquals(selector.getNumSamples(), numSamples2);
        
        //Force IllegalArgumentException
        try{
            selector.setNumSamples(0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}
    }
    
    @Test
    public void testComputeRandomSubsets() throws NotEnoughSamplesException, 
            InvalidSubsetSizeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);        
        int subsetSize = randomizer.nextInt(MIN_SUBSET_SIZE, MAX_SUBSET_SIZE);
        
        SubsetSelector selector = SubsetSelector.create(numSamples);
        int[] result1 = new int[subsetSize];
        selector.computeRandomSubsets(subsetSize, result1);
        int[] result2 = selector.computeRandomSubsets(subsetSize);
        
        //check length of results
        assertEquals(result1.length, subsetSize);
        assertEquals(result2.length, subsetSize);
        
        //check that indices in results are valid
        for(int i = 0; i < subsetSize; i++){
            assertTrue(result1[i] >= 0 && result1[i] < numSamples);
            assertTrue(result2[i] >= 0 && result2[i] < numSamples);
        }
        
        //Force InvalidSubsetSizeException
        try{
            selector.computeRandomSubsets(0);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsets(0, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsets(numSamples + 1, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        
        //Force NotEnoughSamplesException
        try{
            selector.computeRandomSubsets(numSamples + 1);
            fail("NotEnoughSamplesException expected but not thrown");
        }catch(NotEnoughSamplesException e){}
        result2 = new int[numSamples + 1];
        try{
            selector.computeRandomSubsets(numSamples + 1, result2);
            fail("NotEnoughSamplesException expected but not thrown");
        }catch(NotEnoughSamplesException e){}
    }
    
    @Test
    public void testComputeRandomSubsetsInRange() 
            throws NotEnoughSamplesException, InvalidSubsetSizeException, 
            InvalidSubsetRangeException{
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);        
        int subsetSize = randomizer.nextInt(MIN_SUBSET_SIZE, MAX_SUBSET_SIZE);
        int minPos = randomizer.nextInt(0, MIN_SAMPLES - subsetSize - 1);
        int maxPos = randomizer.nextInt(MIN_SAMPLES - 1, numSamples);
        
        SubsetSelector selector = SubsetSelector.create(numSamples);        
        int[] result1 = new int[subsetSize];
        selector.computeRandomSubsetsInRange(minPos, maxPos, subsetSize, false,
                result1);
        int[] result2 = selector.computeRandomSubsetsInRange(minPos, maxPos, 
                subsetSize, false);
        
        //check length of results
        assertEquals(result1.length, subsetSize);
        assertEquals(result2.length, subsetSize);
        
        //check that indices in results are valid
        for(int i = 0; i < subsetSize; i++){
            assertTrue(result1[i] >= minPos && result1[i] < maxPos && 
                    result1[i] >= 0 && result1[i] < numSamples);
            assertTrue(result2[i] >= minPos && result2[i] < maxPos && 
                    result2[i] >= 0 && result2[i] < numSamples);
        }
        
        //Force InvalidSubsetSizeException
        try{
            selector.computeRandomSubsetsInRange(minPos, maxPos, 0, false);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, maxPos, 0, false, 
                    result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1,
                    false, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, 
                    minPos + subsetSize - 1, subsetSize, false);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, 
                    minPos + subsetSize - 1, subsetSize, false, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1,
                    false);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1,
                    false, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        }catch(InvalidSubsetSizeException e){}
        
        //Force NotEnoughSamplesException
        try{
            selector.computeRandomSubsetsInRange(minPos, numSamples + 1, 
                    subsetSize, false);
            fail("NotEnoughSamplesException expected but not thrown");
        }catch(NotEnoughSamplesException e){}
        try{
            selector.computeRandomSubsetsInRange(minPos, numSamples + 1, 
                    subsetSize, false, result2);
            fail("NotEnoughSamplesException expected but not thrown");
        }catch(NotEnoughSamplesException e){}
        
        //Force InvalidSubsetRangeException
        try{
            selector.computeRandomSubsetsInRange(maxPos, minPos, subsetSize, 
                    false);
            fail("InvalidSubsetRangeException expected but not thrown");
        }catch(InvalidSubsetRangeException e){}
        try{
            selector.computeRandomSubsetsInRange(maxPos, minPos, subsetSize, 
                    false, result2);
            fail("InvalidSubsetRangeException expected but not thrown");
        }catch(InvalidSubsetRangeException e){}
    }        
}