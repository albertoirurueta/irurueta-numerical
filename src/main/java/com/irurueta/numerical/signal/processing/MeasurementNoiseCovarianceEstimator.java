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

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.ArrayUtils;
import com.irurueta.algebra.Matrix;

/**
 * Estimates noise covariance matrix for a given set of measures.
 * Covariance matrix is updated each time that measures are added, and it is
 * always the exact running population covariance of all the samples added so
 * far.
 * This class also computes average values of samples to estimate the bias
 * on given noise samples.
 * This class should only be used for samples obtained while system state is
 * held constant.
 */
public class MeasurementNoiseCovarianceEstimator {

    /**
     * Number of measurement vector dimensions (measure parameters).
     */
    private final int mp;

    /**
     * Estimated measurement noise covariance matrix.
     */
    private final Matrix measurementNoiseCov;

    /**
     * Estimated sample average.
     */
    private final double[] sampleAverage;

    /**
     * Number of samples used for estimation.
     */
    private long sampleCount;

    /**
     * A sample after removing the mean estimated before the sample was added
     * (i.e. the difference between the sample and the previous average). This
     * is used internally, and it is kept as an instance variable for reuse
     * purposes.
     */
    private final double[] sampleNoMean;

    /**
     * A sample expressed in matrix form. This is used internally, and it is kept
     * as an instance variable for reuse purposes.
     */
    private final Matrix sampleMatrix;

    /**
     * The transposed sample matrix. This is used internally, and it is kept as
     * an instance variable for reuse purposes.
     */
    private final Matrix transposedSampleMatrix;

    /**
     * A covariance matrix for a single sample. This is used internally, and it
     * is kept as an instance variable for reuse purposes.
     */
    private final Matrix singleCov;

    /**
     * Constructor.
     *
     * @param measureParams number of measurement parameters for each sample
     *                      (i.e. when sampling 3D acceleration samples, this value must be 3, since
     *                      each sample contains acceleration values for x, y, and z axes).
     * @throws IllegalArgumentException  if provided number of measure parameters
     *                                   is less than 1.
     * @throws SignalProcessingException if something fails.
     */
    public MeasurementNoiseCovarianceEstimator(final int measureParams) throws SignalProcessingException {

        if (measureParams < 1) {
            throw new IllegalArgumentException("Measure parameters must be greater than zero");
        }

        try {
            mp = measureParams;
            measurementNoiseCov = new Matrix(mp, mp);
            sampleAverage = new double[mp];

            sampleNoMean = new double[mp];
            sampleMatrix = new Matrix(mp, 1);
            transposedSampleMatrix = new Matrix(1, mp);

            singleCov = new Matrix(mp, mp);
        } catch (final AlgebraException ex) {
            throw new SignalProcessingException(ex);
        }
    }

    /**
     * Updates currently estimated covariance matrix by adding provided sample
     * data.
     * After each update, the returned covariance matrix is the exact running
     * population covariance of all samples added so far (i.e. it matches the
     * batch population covariance, with N as divisor, of those samples), and
     * the sample average is the exact mean of those samples.
     *
     * @param sample sample to be added to update covariance matrix.
     * @return covariance matrix after update.
     * @throws IllegalArgumentException  if provided sample length is not equal
     *                                   to the number of measure parameters set for this instance.
     * @throws SignalProcessingException if something fails.
     */
    public Matrix update(final double[] sample) throws SignalProcessingException {

        if (sample.length != mp) {
            throw new IllegalArgumentException("wrong sample size");
        }

        final var nextCount = sampleCount + 1;

        // compute sample without mean, using the mean prior to this update
        ArrayUtils.subtract(sample, sampleAverage, sampleNoMean);

        // update sample average
        for (int i = 0; i < mp; i++) {
            sampleAverage[i] = (sampleAverage[i] * sampleCount + sample[i]) / nextCount;
        }

        // copy sample without mean into matrix form
        sampleMatrix.setSubmatrix(0, 0, mp - 1, 0, sampleNoMean);
        // compute transpose of matrix form
        sampleMatrix.transpose(transposedSampleMatrix);

        try {
            // compute covariance matrix for a single sample
            sampleMatrix.multiply(transposedSampleMatrix, singleCov);

            // weight by N / (N + 1) to obtain (x - oldMean)(x - newMean)^T
            singleCov.multiplyByScalar((double) sampleCount / nextCount);

            // update covariance matrix
            measurementNoiseCov.multiplyByScalar(sampleCount);
            measurementNoiseCov.add(singleCov);
            measurementNoiseCov.multiplyByScalar(1.0 / nextCount);
        } catch (final AlgebraException ex) {
            throw new SignalProcessingException(ex);
        }

        sampleCount = nextCount;

        return measurementNoiseCov;
    }

    /**
     * Obtains the number of measurement vector dimensions (measure parameters).
     *
     * @return number of measurement vector dimensions (measure parameters).
     */
    public int getMeasureParams() {
        return mp;
    }

    /**
     * Obtains estimated measurement noise covariance matrix.
     *
     * @return estimated measurement noise covariance matrix.
     */
    public Matrix getMeasurementNoiseCov() {
        return measurementNoiseCov;
    }

    /**
     * Obtains estimated sample average.
     *
     * @return estimated sample average.
     */
    public double[] getSampleAverage() {
        return sampleAverage;
    }

    /**
     * Obtains number of samples used for estimation.
     *
     * @return number of samples used for estimation.
     */
    public long getSampleCount() {
        return sampleCount;
    }
}
