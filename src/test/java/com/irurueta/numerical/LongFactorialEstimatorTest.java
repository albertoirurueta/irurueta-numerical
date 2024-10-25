/*
 * Copyright (C) 2023 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class LongFactorialEstimatorTest {

    @Test
    void defaultConstructor_setsExpectedCache() {
        final var estimator = new LongFactorialEstimator();

        assertEquals(21, estimator.getCacheSize());
    }

    @Test
    void constructor_whenSmallValue_throwsException() {
        assertThrows(IllegalArgumentException.class, () -> new LongFactorialEstimator(0));
    }

    @Test
    void constructor_whenLargeValue_throwsException() {
        assertThrows(IllegalArgumentException.class, () -> new LongFactorialEstimator(22));
    }

    @Test
    void factorial_returnsExpectedValue() {
        final var estimator = new LongFactorialEstimator();

        final var cacheSize = estimator.getCacheSize();
        var previous = 0L;
        for (var i = 0; i < cacheSize; i++) {
            final var factorial = estimator.factorial(i);
            final var expected = fact(i);
            assertEquals(expected, factorial, 0.0);
            assertTrue(factorial >= previous);
            previous = factorial;
        }
    }

    @Test
    void factorial_whenSmallValue_throwsException() {
        final var estimator = new LongFactorialEstimator();

        assertThrows(IllegalArgumentException.class, () -> estimator.factorial(-1));
    }

    @Test
    void factorial_whenLargeValue_throwsException() {
        final var estimator = new LongFactorialEstimator();

        final var cacheSize = estimator.getCacheSize();
        assertThrows(IllegalArgumentException.class, () -> estimator.factorial(cacheSize));
    }

    private static long fact(long value) {
        var result = 1L;
        if (value > 0) {
            for (var i = 1; i <= value; i++) {
                result *= i;
            }
        }
        return result;
    }
}