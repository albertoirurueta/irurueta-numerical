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
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class MeasurementNoiseCovarianceEstimatorTest {

    private static final String ACCELERATION_NO_MOTION =
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-no-motion.dat";

    private static final double ABSOLUTE_ERROR = 1e-6;

    private static final double RELATIVE_ERROR = 1e-9;

    private static final double TIGHT_ERROR = 1e-12;

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

        // compute centered batch population covariance
        final var samples = new double[data.numSamples][3];
        for (var i = 0; i < data.numSamples; i++) {
            samples[i][0] = data.accelerationX[i];
            samples[i][1] = data.accelerationY[i];
            samples[i][2] = data.accelerationZ[i];
        }
        final var cov = batchCovariance(samples, sampleAverage2);

        // compare matrices
        assertMatrixEquals(cov, measurementNoiseCov, RELATIVE_ERROR);
    }

    @Test
    void testUpdateKnownSamples() throws SignalProcessingException {
        final var estimator = new MeasurementNoiseCovarianceEstimator(1);

        final var expectedVariances = new double[]{0.0, 1.0, 8.0 / 3.0};
        final var expectedMeans = new double[]{0.0, 1.0, 2.0};
        final var values = new double[]{0.0, 2.0, 4.0};

        for (var i = 0; i < values.length; i++) {
            final var cov = estimator.update(new double[]{values[i]});

            assertEquals(i + 1, estimator.getSampleCount());
            assertEquals(expectedVariances[i], cov.getElementAt(0, 0), TIGHT_ERROR);
            assertEquals(expectedMeans[i], estimator.getSampleAverage()[0], TIGHT_ERROR);
        }
    }

    @Test
    void testUpdateMultivariate() throws SignalProcessingException, WrongSizeException {
        final var estimator = new MeasurementNoiseCovarianceEstimator(3);

        final var samples = new double[][]{
                {1.0, -2.0, 0.5},
                {3.0, 4.0, -1.5},
                {-2.0, 0.0, 2.5},
                {0.5, 1.5, 3.0},
                {4.0, -3.5, 0.0},
                {2.5, 2.0, -2.0}
        };

        for (var n = 1; n <= samples.length; n++) {
            final var cov = estimator.update(samples[n - 1]);

            final var subset = Arrays.copyOf(samples, n);
            final var mean = new double[3];
            for (final var sample : subset) {
                for (var j = 0; j < 3; j++) {
                    mean[j] += sample[j] / n;
                }
            }

            assertArrayEquals(mean, estimator.getSampleAverage(), TIGHT_ERROR);
            assertMatrixEquals(batchCovariance(subset, mean), cov, TIGHT_ERROR);

            // covariance must be symmetric
            for (var i = 0; i < 3; i++) {
                for (var j = i + 1; j < 3; j++) {
                    assertEquals(cov.getElementAt(i, j), cov.getElementAt(j, i), TIGHT_ERROR);
                }
            }
        }
    }

    @Test
    void testUpdateWrongSampleSize() throws SignalProcessingException {
        final var estimator = new MeasurementNoiseCovarianceEstimator(2);

        final var sample = new double[3];
        assertThrows(IllegalArgumentException.class, () -> estimator.update(sample));
    }

    private static Matrix batchCovariance(final double[][] samples, final double[] mean) throws WrongSizeException {
        final var dims = mean.length;
        final var cov = new Matrix(dims, dims);
        for (final var sample : samples) {
            for (var i = 0; i < dims; i++) {
                for (var j = 0; j < dims; j++) {
                    cov.setElementAt(i, j, cov.getElementAt(i, j)
                            + (sample[i] - mean[i]) * (sample[j] - mean[j]) / samples.length);
                }
            }
        }
        return cov;
    }

    private static void assertMatrixEquals(final Matrix expected, final Matrix actual, final double relativeError) {
        assertEquals(expected.getRows(), actual.getRows());
        assertEquals(expected.getColumns(), actual.getColumns());
        for (var i = 0; i < expected.getRows(); i++) {
            for (var j = 0; j < expected.getColumns(); j++) {
                final var e = expected.getElementAt(i, j);
                assertEquals(e, actual.getElementAt(i, j), relativeError * Math.max(1.0, Math.abs(e)));
            }
        }
    }
}
