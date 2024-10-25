/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.signal.processing;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

class MeasurementNoiseCovarianceEstimatorTest {

    private static final String ACCELERATION_NO_MOTION =
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-no-motion.dat";

    private static final double ABSOLUTE_ERROR = 1e-6;

    private static final double LARGE_ABSOLUTE_ERROR = 1e-3;

    @Test
    void testConstructor() throws SignalProcessingException {

        var estimator = new MeasurementNoiseCovarianceEstimator(3);

        // check correctness
        assertEquals(3, estimator.getMeasureParams());
        assertEquals(3, estimator.getMeasurementNoiseCov().getRows());
        assertEquals(3, estimator.getMeasurementNoiseCov().getColumns());
        assertEquals(3, estimator.getSampleAverage().length);
        assertEquals(0, estimator.getSampleCount());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> new MeasurementNoiseCovarianceEstimator(0));
    }

    @Test
    void testUpdate() throws SignalProcessingException, IOException, WrongSizeException {

        final var estimator = new MeasurementNoiseCovarianceEstimator(3);

        final var f = new File(ACCELERATION_NO_MOTION);
        final var data = AccelerationFileLoader.load(f);

        final var sample = new double[3];

        for (var i = 0; i < data.numSamples; i++) {
            sample[0] = data.accelerationX[i];
            sample[1] = data.accelerationY[i];
            sample[2] = data.accelerationZ[i];

            estimator.update(sample);
        }

        final var measurementNoiseCov = estimator.getMeasurementNoiseCov();
        final var sampleAverage = estimator.getSampleAverage();

        assertEquals(data.numSamples, estimator.getSampleCount());

        // compute sample average
        final var sampleAverage2 = new double[3];
        for (var i = 0; i < data.numSamples; i++) {
            sampleAverage2[0] += data.accelerationX[i] / data.numSamples;
            sampleAverage2[1] += data.accelerationY[i] / data.numSamples;
            sampleAverage2[2] += data.accelerationZ[i] / data.numSamples;
        }

        // check correctness of average
        assertArrayEquals(sampleAverage, sampleAverage2, ABSOLUTE_ERROR);

        // compute covariance

        // copy all acceleration samples into a matrix (each column contains a
        // sample)
        final var samples = new Matrix(3, data.numSamples);
        for (var i = 0; i < data.numSamples; i++) {
            samples.setElementAt(0, i, data.accelerationX[i]);
            samples.setElementAt(1, i, data.accelerationY[i]);
            samples.setElementAt(2, i, data.accelerationZ[i]);
        }

        final var transSamples = samples.transposeAndReturnNew();

        final var cov = samples.multiplyAndReturnNew(transSamples);
        cov.multiplyByScalar(1.0 / data.numSamples);

        // compare matrices
        for (var i = 0; i < 3; i++) {
            for (var j = 0; j < 3; j++) {
                assertEquals(measurementNoiseCov.getElementAt(i, j), cov.getElementAt(i, j), LARGE_ABSOLUTE_ERROR);
            }
        }
    }
}
