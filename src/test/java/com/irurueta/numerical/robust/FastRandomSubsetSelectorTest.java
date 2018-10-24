/*
 * Copyright (C) 2013 Alberto Irurueta Carro (alberto@irurueta.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.irurueta.numerical.robust;

import com.irurueta.statistics.UniformRandomizer;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.*;

public class FastRandomSubsetSelectorTest {
    
    private static final int MIN_SAMPLES = 100;
    private static final int MAX_SAMPLES = 500;
    
    private static final int MIN_SUBSET_SIZE = 5;
    private static final int MAX_SUBSET_SIZE = 10;
    
    public FastRandomSubsetSelectorTest() { }
    
    @BeforeClass
    public static void setUpClass() { }
    
    @AfterClass
    public static void tearDownClass() { }
    
    @Before
    public void setUp() { }
    
    @After
    public void tearDown() { }
    
    @Test
    public void testConstructor() {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        
        //Test constructor with number of samples
        FastRandomSubsetSelector selector = new FastRandomSubsetSelector(
                numSamples);
        assertNotNull(selector.getRandomizer());
        assertEquals(selector.getType(), 
                SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);
        assertEquals(selector.getNumSamples(), numSamples);
        
        //Force IllegalArgumentException
        selector = null;
        try {
            selector = new FastRandomSubsetSelector(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(selector);
        
        //Test constructor with number of samples and indicating whether 
        //randomizer is initialized with system timer
        selector = new FastRandomSubsetSelector(numSamples, false);
        assertNotNull(selector.getRandomizer());
        assertEquals(selector.getType(),
                SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);
        assertEquals(selector.getNumSamples(), numSamples);
        
        //Force IllegalArgumentException
        selector = null;
        try {
            selector = new FastRandomSubsetSelector(0, false);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
        assertNull(selector);
    }
    
    @Test
    public void testGetType() {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        
        FastRandomSubsetSelector selector = new FastRandomSubsetSelector(
                numSamples);
        assertEquals(selector.getType(), 
                SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR);        
    }
    
    @Test
    public void testGetSetNumSamples() {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples1 = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        int numSamples2 = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        
        FastRandomSubsetSelector selector = new FastRandomSubsetSelector(
                numSamples1);
        //check correctness
        assertEquals(selector.getNumSamples(), numSamples1);
        
        //set new value
        selector.setNumSamples(numSamples2);
        
        //check correctness
        assertEquals(selector.getNumSamples(), numSamples2);
        
        //Force IllegalArgumentException
        try {
            selector.setNumSamples(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) { }
    }
    
    @Test
    public void testComputeRandomSubsets() throws NotEnoughSamplesException, 
            InvalidSubsetSizeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);        
        int subsetSize = randomizer.nextInt(MIN_SUBSET_SIZE, MAX_SUBSET_SIZE);
        
        FastRandomSubsetSelector selector = new FastRandomSubsetSelector(
                numSamples);
        int[] result1 = new int[subsetSize];
        selector.computeRandomSubsets(subsetSize, result1);
        int[] result2 = selector.computeRandomSubsets(subsetSize);
        
        //check length of results
        assertEquals(result1.length, subsetSize);
        assertEquals(result2.length, subsetSize);
        
        //check that indices in results are valid
        for (int i = 0; i < subsetSize; i++) {
            assertTrue(result1[i] >= 0 && result1[i] < numSamples);
            assertTrue(result2[i] >= 0 && result2[i] < numSamples);
        }
        
        //Force InvalidSubsetSizeException
        try {
            selector.computeRandomSubsets(0);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsets(0, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsets(numSamples + 1, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        
        //Force NotEnoughSamplesException
        try {
            selector.computeRandomSubsets(numSamples + 1);
            fail("NotEnoughSamplesException expected but not thrown");
        } catch (NotEnoughSamplesException ignore) { }
        result2 = new int[numSamples + 1];
        try {
            selector.computeRandomSubsets(numSamples + 1, result2);
            fail("NotEnoughSamplesException expected but not thrown");
        } catch (NotEnoughSamplesException ignore) { }
    }
    
    @Test
    public void testComputeRandomSubsetsInRange() 
            throws NotEnoughSamplesException, InvalidSubsetSizeException, 
            InvalidSubsetRangeException {
        UniformRandomizer randomizer = new UniformRandomizer(new Random());
        int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);        
        int subsetSize = randomizer.nextInt(MIN_SUBSET_SIZE, MAX_SUBSET_SIZE);
        int minPos = randomizer.nextInt(0, MIN_SAMPLES - subsetSize - 1);
        int maxPos = randomizer.nextInt(MIN_SAMPLES - 1, numSamples);
        
        FastRandomSubsetSelector selector = new FastRandomSubsetSelector(
                numSamples);
        int[] result1 = new int[subsetSize];
        selector.computeRandomSubsetsInRange(minPos, maxPos, subsetSize, false,
                result1);
        int[] result2 = selector.computeRandomSubsetsInRange(minPos, maxPos, 
                subsetSize, false);
        
        //check length of results
        assertEquals(result1.length, subsetSize);
        assertEquals(result2.length, subsetSize);
        
        //check that indices in results are valid
        for (int i = 0; i < subsetSize; i++) {
            assertTrue(result1[i] >= minPos && result1[i] < maxPos && 
                    result1[i] >= 0 && result1[i] < numSamples);
            assertTrue(result2[i] >= minPos && result2[i] < maxPos && 
                    result2[i] >= 0 && result2[i] < numSamples);
        }
        
        //Force InvalidSubsetSizeException
        try {
            selector.computeRandomSubsetsInRange(minPos, maxPos, 0, false);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(minPos, maxPos, 0, false, 
                    result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1,
                    false, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try{
            selector.computeRandomSubsetsInRange(minPos, 
                    minPos + subsetSize - 1, subsetSize, false);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(minPos, 
                    minPos + subsetSize - 1, subsetSize, false, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1,
                    false);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1,
                    false, result2);
            fail("InvalidSubsetSizeException expected but not thrown");
        } catch (InvalidSubsetSizeException ignore) { }
        
        //Force NotEnoughSamplesException
        try {
            selector.computeRandomSubsetsInRange(minPos, numSamples + 1, 
                    subsetSize, false);
            fail("NotEnoughSamplesException expected but not thrown");
        } catch (NotEnoughSamplesException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(minPos, numSamples + 1, 
                    subsetSize, false, result2);
            fail("NotEnoughSamplesException expected but not thrown");
        } catch (NotEnoughSamplesException ignore) { }
        
        //Force InvalidSubsetRangeException
        try {
            selector.computeRandomSubsetsInRange(maxPos, minPos, subsetSize, 
                    false);
            fail("InvalidSubsetRangeException expected but not thrown");
        } catch (InvalidSubsetRangeException ignore) { }
        try {
            selector.computeRandomSubsetsInRange(maxPos, minPos, subsetSize, 
                    false, result2);
            fail("InvalidSubsetRangeException expected but not thrown");
        } catch (InvalidSubsetRangeException ignore) { }
    }    
}