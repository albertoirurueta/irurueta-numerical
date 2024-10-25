/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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

import com.irurueta.sorting.SortingException;
import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class WeightSelectionTest {

    private static final int MIN_LENGTH = 6;
    private static final int MAX_LENGTH = 50;

    private static final double MIN_RANDOM_VALUE = 0.0;
    private static final double MAX_RANDOM_VALUE = 1.0;

    @Test
    void testSelectWeights() throws SortingException {

        final var randomizer = new UniformRandomizer();
        final var length = randomizer.nextInt(MIN_LENGTH, MAX_LENGTH);

        final var weights = new double[length];
        randomizer.fill(weights, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        // test with sorting disabled and num selected lower than length
        var maxPoints = length - 1;
        var selection = WeightSelection.selectWeights(weights, false, maxPoints);
        assertEquals(selection.getNumSelected(), maxPoints);
        // check first maxPoints values are true
        for (var i = 0; i < maxPoints; i++) {
            assertTrue(selection.getSelected()[i]);
        }
        for (var i = maxPoints; i < length; i++) {
            assertFalse(selection.getSelected()[i]);
        }

        // test with sorting disabled and num selected greater than length
        maxPoints = length + 1;
        selection = WeightSelection.selectWeights(weights, false, maxPoints);
        assertEquals(selection.getNumSelected(), length);
        // check all values are true
        for (var i = 0; i < length; i++) {
            assertTrue(selection.getSelected()[i]);
        }

        // test with sorting enabled
        maxPoints = length - 1;
        selection = WeightSelection.selectWeights(weights, true, maxPoints);
        assertEquals(selection.getNumSelected(), maxPoints);
        // check all values are true
        var counter = 0;
        for (var i = 0; i < length; i++) {
            if (selection.getSelected()[i]) {
                counter++;
            }
        }
        assertEquals(counter, selection.getNumSelected());
    }
}
