package com.irurueta.numerical.interpolation;

import static org.junit.jupiter.api.Assertions.*;

import com.irurueta.statistics.UniformRandomizer;

import org.junit.jupiter.api.Test;

class ThinPlateRadialBasisFunctionTest {

    @Test
    void evaluate_whenDefaultScale_returnsExpectedValue() {
        final var rbf = new ThinPlateRadialBasisFunction();

        final var randomizer = new UniformRandomizer();
        final var r = randomizer.nextDouble(-1.0, 1.0);
        final var expected = r <= 0.0 ? 0.0 : r * r * Math.log(r);

        assertEquals(expected, rbf.evaluate(r), 0.0);
    }

    @Test
    void evaluate_whenRandomScale_returnsExpectedValue() {
        final var randomizer = new UniformRandomizer();
        final var scale = randomizer.nextDouble();

        final var rbf = new ThinPlateRadialBasisFunction(scale);

        final var r = randomizer.nextDouble();
        final var expected = r <= 0.0 ? 0.0 : r * r * Math.log(r / scale);

        assertEquals(expected, rbf.evaluate(r), 0.0);
    }
}