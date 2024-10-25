package com.irurueta.numerical.interpolation;

import static org.junit.jupiter.api.Assertions.*;

import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class MultiQuadricRadialBasisFunctionTest {

    @Test
    void evaluate_whenDefaultScale_returnsExpectedValue() {
        final var rbf = new MultiQuadricRadialBasisFunction();

        final var randomizer = new UniformRandomizer();
        final var r = randomizer.nextDouble();
        final var expected = Math.sqrt((r * r) + 1.0);

        assertEquals(expected, rbf.evaluate(r), 0.0);
    }

    @Test
    void evaluate_whenRandomScale_returnsExpectedValue() {
        final var randomizer = new UniformRandomizer();
        final var scale = randomizer.nextDouble();

        final var rbf = new MultiQuadricRadialBasisFunction(scale);

        final var r = randomizer.nextDouble();
        final var expected = Math.sqrt((r * r) + (scale * scale));

        assertEquals(expected, rbf.evaluate(r), 0.0);
    }
}