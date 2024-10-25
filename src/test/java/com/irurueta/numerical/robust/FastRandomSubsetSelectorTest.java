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
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class FastRandomSubsetSelectorTest {

    private static final int MIN_SAMPLES = 100;
    private static final int MAX_SAMPLES = 500;

    private static final int MIN_SUBSET_SIZE = 5;
    private static final int MAX_SUBSET_SIZE = 10;

    @Test
    void testConstructor() {
        final var randomizer = new UniformRandomizer();
        final var numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);

        // Test constructor with number of samples
        FastRandomSubsetSelector selector = new FastRandomSubsetSelector(numSamples);
        assertNotNull(selector.getRandomizer());
        assertEquals(SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR, selector.getType());
        assertEquals(selector.getNumSamples(), numSamples);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new FastRandomSubsetSelector(0));

        // Test constructor with number of samples and indicating whether
        // randomizer is initialized with system timer
        selector = new FastRandomSubsetSelector(numSamples, false);
        assertNotNull(selector.getRandomizer());
        assertEquals(SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR, selector.getType());
        assertEquals(numSamples, selector.getNumSamples());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new FastRandomSubsetSelector(0,
                false));
    }

    @Test
    void testGetType() {
        final var randomizer = new UniformRandomizer();
        final int numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);

        final var selector = new FastRandomSubsetSelector(numSamples);
        assertEquals(SubsetSelectorType.FAST_RANDOM_SUBSET_SELECTOR, selector.getType());
    }

    @Test
    void testGetSetNumSamples() {
        final var randomizer = new UniformRandomizer();
        final var numSamples1 = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        final var numSamples2 = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);

        final var selector = new FastRandomSubsetSelector(numSamples1);
        // check correctness
        assertEquals(selector.getNumSamples(), numSamples1);

        // set new value
        selector.setNumSamples(numSamples2);

        // check correctness
        assertEquals(numSamples2, selector.getNumSamples());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> selector.setNumSamples(0));
    }

    @Test
    void testComputeRandomSubsets() throws NotEnoughSamplesException, InvalidSubsetSizeException {
        final var randomizer = new UniformRandomizer();
        final var numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        final var subsetSize = randomizer.nextInt(MIN_SUBSET_SIZE, MAX_SUBSET_SIZE);

        final var selector = new FastRandomSubsetSelector(numSamples);
        final var result1 = new int[subsetSize];
        selector.computeRandomSubsets(subsetSize, result1);
        var result2 = selector.computeRandomSubsets(subsetSize);

        // check length of results
        assertEquals(result1.length, subsetSize);
        assertEquals(result2.length, subsetSize);

        // check that indices in results are valid
        for (var i = 0; i < subsetSize; i++) {
            assertTrue(result1[i] >= 0 && result1[i] < numSamples);
            assertTrue(result2[i] >= 0 && result2[i] < numSamples);
        }

        // Force InvalidSubsetSizeException
        assertThrows(InvalidSubsetSizeException.class, () -> selector.computeRandomSubsets(0));
        final var finalResult = result2;
        assertThrows(InvalidSubsetSizeException.class, () -> selector.computeRandomSubsets(0, finalResult));
        assertThrows(InvalidSubsetSizeException.class, () -> selector.computeRandomSubsets(numSamples + 1,
                finalResult));

        // Force NotEnoughSamplesException
        assertThrows(NotEnoughSamplesException.class, () -> selector.computeRandomSubsets(numSamples + 1));
        final var finalResult2 = new int[numSamples + 1];
        assertThrows(NotEnoughSamplesException.class, () -> selector.computeRandomSubsets(numSamples + 1,
                finalResult2));
    }

    @Test
    void testComputeRandomSubsetsInRange() throws NotEnoughSamplesException, InvalidSubsetSizeException,
            InvalidSubsetRangeException {
        final var randomizer = new UniformRandomizer();
        final var numSamples = randomizer.nextInt(MIN_SAMPLES, MAX_SAMPLES);
        final var subsetSize = randomizer.nextInt(MIN_SUBSET_SIZE, MAX_SUBSET_SIZE);
        final var minPos = randomizer.nextInt(0, MIN_SAMPLES - subsetSize - 1);
        final var maxPos = randomizer.nextInt(MIN_SAMPLES - 1, numSamples);

        final var selector = new FastRandomSubsetSelector(numSamples);
        final var result1 = new int[subsetSize];
        selector.computeRandomSubsetsInRange(minPos, maxPos, subsetSize, false, result1);
        final var result2 = selector.computeRandomSubsetsInRange(minPos, maxPos, subsetSize, false);

        // check length of results
        assertEquals(result1.length, subsetSize);
        assertEquals(result2.length, subsetSize);

        // check that indices in results are valid
        for (var i = 0; i < subsetSize; i++) {
            assertTrue(result1[i] >= minPos && result1[i] < maxPos && result1[i] >= 0
                    && result1[i] < numSamples);
            assertTrue(result2[i] >= minPos && result2[i] < maxPos && result2[i] >= 0
                    && result2[i] < numSamples);
        }

        // Force InvalidSubsetSizeException
        assertThrows(InvalidSubsetSizeException.class,
                () -> selector.computeRandomSubsetsInRange(minPos, maxPos, 0, false));
        assertThrows(InvalidSubsetSizeException.class,
                () -> selector.computeRandomSubsetsInRange(minPos, maxPos, 0, false, result2));
        assertThrows(InvalidSubsetSizeException.class,
                () -> selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1, false,
                        result2));
        assertThrows(InvalidSubsetSizeException.class, () -> selector.computeRandomSubsetsInRange(minPos,
                minPos + subsetSize - 1, subsetSize, false));
        assertThrows(InvalidSubsetSizeException.class, () -> selector.computeRandomSubsetsInRange(minPos,
                minPos + subsetSize - 1, subsetSize, false, result2));
        assertThrows(InvalidSubsetSizeException.class,
                () -> selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1, false));
        assertThrows(InvalidSubsetSizeException.class,
                () -> selector.computeRandomSubsetsInRange(minPos, maxPos, numSamples + 1, false,
                        result2));

        // Force NotEnoughSamplesException
        assertThrows(NotEnoughSamplesException.class, () -> selector.computeRandomSubsetsInRange(minPos,
                numSamples + 1, subsetSize, false));
        assertThrows(NotEnoughSamplesException.class, () -> selector.computeRandomSubsetsInRange(minPos,
                numSamples + 1, subsetSize, false, result2));

        // Force InvalidSubsetRangeException
        assertThrows(InvalidSubsetRangeException.class, () -> selector.computeRandomSubsetsInRange(maxPos, minPos,
                subsetSize, false));
        assertThrows(InvalidSubsetRangeException.class, () -> selector.computeRandomSubsetsInRange(maxPos, minPos,
                subsetSize, false, result2));
    }
}