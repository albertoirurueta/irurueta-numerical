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
package com.irurueta.numerical.interpolation;

import static org.junit.jupiter.api.Assertions.assertEquals;

import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class InverseMultiQuadricRadialBasisFunctionTest {

    @Test
    void evaluate_whenDefaultScale_returnsExpectedValue() {
        final var rbf = new InverseMultiQuadricRadialBasisFunction();

        final var randomizer = new UniformRandomizer();
        final var r = randomizer.nextDouble();
        final var expected = 1.0 / Math.sqrt((r * r) + 1.0);

        assertEquals(expected, rbf.evaluate(r), 0.0);
    }

    @Test
    void evaluate_whenRandomScale_returnsExpectedValue() {
        final var randomizer = new UniformRandomizer();
        final var scale = randomizer.nextDouble();

        final var rbf = new InverseMultiQuadricRadialBasisFunction(scale);

        final var r = randomizer.nextDouble();
        final var expected = 1.0 / Math.sqrt((r * r) + (scale * scale));

        assertEquals(expected, rbf.evaluate(r), 0.0);
    }
}