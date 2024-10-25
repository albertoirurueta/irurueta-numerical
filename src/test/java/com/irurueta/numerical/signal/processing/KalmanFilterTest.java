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
import com.irurueta.numerical.signal.processing.AccelerationFileLoader.Data;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

class KalmanFilterTest {

    private static final int N_SAMPLES = 500;
    private static final boolean WRITE_TO_CONSOLE = false;
    private static final boolean DO_NOT_SKIP = true;
    private static final int TIMES = 100;

    private static final String ACCELERATION_MOTION =
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-motion.dat";

    private static final String ACCELERATION_NO_MOTION =
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-no-motion.dat";

    private static final String CSV_FILE_MOTION =
            "./src/test/java/com/irurueta/numerical/signal/processing/motion.csv";

    private static final String CSV_FILE_NO_MOTION =
            "./src/test/java/com/irurueta/numerical/signal/processing/no-motion.csv";

    private static final String ACCELERATION_MOTION_FAST =
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-motion-fast.dat";

    private static final String ACCELERATION_NO_MOTION_FAST =
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-no-motion-fast.dat";

    private static final String CSV_FILE_MOTION_FAST =
            "./src/test/java/com/irurueta/numerical/signal/processing/motion-fast.csv";

    private static final String CSV_FILE_NO_MOTION_FAST =
            "./src/test/java/com/irurueta/numerical/signal/processing/no-motion-fast.csv";

    @Test
    void testConstructor() throws SignalProcessingException {
        var filter = new KalmanFilter(6, 9, -1);

        // check correctness
        assertEquals(6, filter.getDynamicParameters());
        assertEquals(9, filter.getMeasureParameters());
        // and because control parameters were set to negative value...
        assertEquals(6, filter.getControlParameters());

        assertEquals(6, filter.getStatePre().getRows());
        assertEquals(1, filter.getStatePre().getColumns());

        assertEquals(6, filter.getStatePost().getRows());
        assertEquals(1, filter.getStatePost().getColumns());

        assertEquals(6, filter.getTransitionMatrix().getRows());
        assertEquals(6, filter.getTransitionMatrix().getColumns());

        assertEquals(6, filter.getProcessNoiseCov().getRows());
        assertEquals(6, filter.getProcessNoiseCov().getColumns());

        assertEquals(9, filter.getMeasurementMatrix().getRows());
        assertEquals(6, filter.getMeasurementMatrix().getColumns());

        assertEquals(9, filter.getMeasurementNoiseCov().getRows());
        assertEquals(9, filter.getMeasurementNoiseCov().getColumns());

        assertEquals(6, filter.getErrorCovPre().getRows());
        assertEquals(6, filter.getErrorCovPre().getColumns());

        assertEquals(6, filter.getErrorCovPost().getRows());
        assertEquals(6, filter.getErrorCovPost().getColumns());

        assertEquals(6, filter.getGain().getRows());
        assertEquals(9, filter.getGain().getColumns());

        assertEquals(6, filter.getControlMatrix().getRows());
        assertEquals(6, filter.getControlMatrix().getColumns());

        // test with control parameters
        filter = new KalmanFilter(6, 9, 1);

        // check correctness
        assertEquals(6, filter.getDynamicParameters());
        assertEquals(9, filter.getMeasureParameters());
        assertEquals(1, filter.getControlParameters());

        assertEquals(6, filter.getStatePre().getRows());
        assertEquals(1, filter.getStatePre().getColumns());

        assertEquals(6, filter.getStatePost().getRows());
        assertEquals(1, filter.getStatePost().getColumns());

        assertEquals(6, filter.getTransitionMatrix().getRows());
        assertEquals(6, filter.getTransitionMatrix().getColumns());

        assertEquals(6, filter.getProcessNoiseCov().getRows());
        assertEquals(6, filter.getProcessNoiseCov().getColumns());

        assertEquals(9, filter.getMeasurementMatrix().getRows());
        assertEquals(6, filter.getMeasurementMatrix().getColumns());

        assertEquals(9, filter.getMeasurementNoiseCov().getRows());
        assertEquals(9, filter.getMeasurementNoiseCov().getColumns());

        assertEquals(6, filter.getErrorCovPre().getRows());
        assertEquals(6, filter.getErrorCovPre().getColumns());

        assertEquals(6, filter.getErrorCovPost().getRows());
        assertEquals(6, filter.getErrorCovPost().getColumns());

        assertEquals(6, filter.getGain().getRows());
        assertEquals(9, filter.getGain().getColumns());

        assertEquals(6, filter.getControlMatrix().getRows());
        assertEquals(1, filter.getControlMatrix().getColumns());

        // test without control parameters
        filter = new KalmanFilter(6, 9);

        // check correctness
        assertEquals(6, filter.getDynamicParameters());
        assertEquals(9, filter.getMeasureParameters());
        assertEquals(0, filter.getControlParameters());

        assertEquals(6, filter.getStatePre().getRows());
        assertEquals(1, filter.getStatePre().getColumns());

        assertEquals(6, filter.getStatePost().getRows());
        assertEquals(1, filter.getStatePost().getColumns());

        assertEquals(6, filter.getTransitionMatrix().getRows());
        assertEquals(6, filter.getTransitionMatrix().getColumns());

        assertEquals(6, filter.getProcessNoiseCov().getRows());
        assertEquals(6, filter.getProcessNoiseCov().getColumns());

        assertEquals(9, filter.getMeasurementMatrix().getRows());
        assertEquals(6, filter.getMeasurementMatrix().getColumns());

        assertEquals(9, filter.getMeasurementNoiseCov().getRows());
        assertEquals(9, filter.getMeasurementNoiseCov().getColumns());

        assertEquals(6, filter.getErrorCovPre().getRows());
        assertEquals(6, filter.getErrorCovPre().getColumns());

        assertEquals(6, filter.getErrorCovPost().getRows());
        assertEquals(6, filter.getErrorCovPost().getColumns());

        assertEquals(6, filter.getGain().getRows());
        assertEquals(9, filter.getGain().getColumns());

        assertNull(filter.getControlMatrix());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class,
                () -> new KalmanFilter(-1, 9, 1));
        assertThrows(IllegalArgumentException.class,
                () -> new KalmanFilter(6, -1, 1));
    }

    @Test
    void testGetSetMeasureParameters() throws SignalProcessingException {
        final var filter = new KalmanFilter(6, 9, -1);

        assertEquals(9, filter.getMeasureParameters());

        // set new value
        filter.setMeasureParameters(10);

        // check correctness
        assertEquals(10, filter.getMeasureParameters());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> filter.setMeasureParameters(0));
    }

    @Test
    void testPredictAndCorrectAcceleration() throws SignalProcessingException, WrongSizeException {
        var numValid = 0;
        for (var t = 0; t < TIMES; t++) {

            // using a Kalman filter we take noisy measures of acceleration (for
            // simplicity we only consider one dimension).
            // Using acceleration samples we want to obtain the state of speed and
            // position
            final var kalman = new KalmanFilter(1, 1);

            final var rand = new Random();
            // constant acceleration
            final var acceleration = rand.nextDouble();
            double sampleNoise;
            var correctedAcceleration = 0.0;

            // setup Kalman filter
            //[ acceleration]
            Matrix correctedState;

            // measurement (acceleration)
            final var measurement = new Matrix(1, 1);
            measurement.setElementAtIndex(0, acceleration);

            // transitions for acceleration
            // [1]
            final var transitionMatrix = new Matrix(1, 1);
            transitionMatrix.setElementAtIndex(0, 1.0);
            kalman.setTransitionMatrix(transitionMatrix);

            if (WRITE_TO_CONSOLE) {
                System.out.println("no;truePosition;trueSpeed;trueAcceleration;" +
                        "predictedPosition;predictedSpeed;predictedAcceleration;" +
                        "correctedPosition;correctedSpeed;correctedAcceleration;");
            }

            // assume each sample happens every second
            for (var i = 0; i < N_SAMPLES; i++) {
                sampleNoise = 1e-3 * rand.nextGaussian();

                // take measure with noise
                measurement.setElementAtIndex(0, measurement.getElementAtIndex(0) + sampleNoise);

                // correct
                correctedState = kalman.correct(measurement);

                correctedAcceleration = correctedState.getElementAtIndex(0);

                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";;;" + acceleration + ";" +
                            //";;" + predictedAcceleration + ";" +
                            ";;" + correctedAcceleration);
                }
            }

            if (Math.abs(acceleration - correctedAcceleration) > KalmanFilter.DEFAULT_MEASUREMENT_NOISE_VARIANCE) {
                continue;
            }
            assertEquals(acceleration, correctedAcceleration, KalmanFilter.DEFAULT_MEASUREMENT_NOISE_VARIANCE);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    void testPredictAndCorrectPositionSpeedAndAcceleration() throws SignalProcessingException, WrongSizeException {

        // using a Kalman filter we take noisy measures of acceleration (for
        // simplicity we only consider one dimension).
        // Using acceleration samples we want to obtain the state of speed and
        // position
        final var kalman = new KalmanFilter(3, 1);

        final var rand = new Random();
        // constant acceleration
        final var acceleration = rand.nextDouble();
        var trueSpeed = 0.0;
        var truePosition = 0.0;
        var rawSpeed = 0.0;
        var rawPosition = 0.0;
        double rawAcceleration;
        double sampleNoise;
        var correctedAcceleration = 0.0;

        // setup Kalman filter
        Matrix predictedState; //[position, speed, acceleration]
        Matrix correctedState; //[position, speed, acceleration]

        final var measurement = new Matrix(1, 1); //measurement (acceleration)
        measurement.setElementAtIndex(0, acceleration);

        // assuming delta time between samples of 1 second
        final var deltaTime = 1.0;
        // transitions for position, speed and acceleration
        // [1,   deltaTime,    deltaTime^2/2]
        // [0,   1,            deltaTime]
        // [0,   0,            1]
        final var transitionMatrix = new Matrix(3, 3);
        // column order
        transitionMatrix.setSubmatrix(0, 0, 2, 2, new double[]
                {
                        1.0, 0.0, 0.0,
                        deltaTime, 1.0, 0.0,
                        0.5 * deltaTime * deltaTime, deltaTime, 1.0
                }, true);
        kalman.setTransitionMatrix(transitionMatrix);

        final var processNoiseStd = Math.sqrt(1e-6);
        final var processNoiseVariance = processNoiseStd * processNoiseStd;

        // process noise covariance =
        // processNoiseVariance * [deltaTime^4/4,   deltaTime^3/2,  deltaTime^2/2]
        //                        [deltaTime^3/2,   deltaTime^2,    deltaTime]
        //                        [deltaTime^2/2,   deltaTime,      1]
        final var processNoiseCov = new Matrix(3, 3);
        // column order
        processNoiseCov.setSubmatrix(0, 0, 2, 2, new double[]
                {
                        processNoiseVariance * Math.pow(deltaTime, 4.0) * 0.25,
                        processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 2.0),
                        processNoiseVariance * deltaTime,
                        processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                        processNoiseVariance * deltaTime,
                        processNoiseVariance
                }, true);
        kalman.setProcessNoiseCov(processNoiseCov);

        final var measureNoiseStd = Math.sqrt(1e-1);
        final var measureNoiseVariance = measureNoiseStd * measureNoiseStd;

        // measurement noise covariance = [measureNoiseVariance]
        final var measurementNoiseCov = new Matrix(1, 1);
        measurementNoiseCov.setElementAtIndex(0, measureNoiseVariance);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        final var measurementMatrix = new Matrix(1, 3);
        measurementMatrix.setSubmatrix(0, 0, 0, 2,
                new double[]{0.0, 0.0, 1.0}, true);
        kalman.setMeasurementMatrix(measurementMatrix);

        if (WRITE_TO_CONSOLE) {
            System.out.println("no;truePosition;trueSpeed;trueAcceleration;" +
                    "predictedPosition;predictedSpeed;predictedAcceleration;" +
                    "correctedPosition;correctedSpeed;correctedAcceleration;" +
                    "rawPosition;rawSpeed;rawAcceleration");
        }

        // assume each sample happens every second
        for (var i = 0; i < N_SAMPLES; i++) {
            trueSpeed += acceleration;
            truePosition += trueSpeed;

            sampleNoise = measureNoiseStd * rand.nextGaussian();

            rawAcceleration = acceleration + sampleNoise;
            rawSpeed += rawAcceleration;
            rawPosition += rawSpeed;

            predictedState = kalman.predict();
            assertNotNull(predictedState);

            // correct state for the 75% of samples if  DO_NOT_SKIP is enabled
            if (rand.nextDouble() <= 0.75 || DO_NOT_SKIP) {
                // take measure with noise
                measurement.setElementAtIndex(0, rawAcceleration);

                // correct
                correctedState = kalman.correct(measurement);

                correctedAcceleration = correctedState.getElementAtIndex(2);

                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";" + truePosition + ";" + trueSpeed + ";" + acceleration + ";" +
                            correctedAcceleration + ";" + rawPosition + ";" + rawSpeed + ";" + rawAcceleration);
                }

            } else {
                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";" + truePosition + ";" + trueSpeed + ";" + acceleration + ";" +
                            ";;;" + rawPosition + ";" + rawSpeed + ";" + rawAcceleration);
                }
            }
        }

        // check that the filter state converged to the true acceleration
        assertEquals(acceleration, correctedAcceleration, measureNoiseVariance);
    }

    @Test
    void testPredictAndCorrectPositionSpeedAndAcceleration3D() throws SignalProcessingException, WrongSizeException {

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final var kalman = new KalmanFilter(9, 3);

        final var rand = new Random();
        // constant acceleration
        final var accelerationX = rand.nextDouble();
        final var accelerationY = rand.nextDouble();
        final var accelerationZ = rand.nextDouble();
        var trueSpeedX = 0.0;
        var trueSpeedY = 0.0;
        var trueSpeedZ = 0.0;
        var truePositionX = 0.0;
        var truePositionY = 0.0;
        var truePositionZ = 0.0;
        var rawSpeedX = 0.0;
        var rawSpeedY = 0.0;
        var rawSpeedZ = 0.0;
        var rawPositionX = 0.0;
        var rawPositionY = 0.0;
        var rawPositionZ = 0.0;
        double rawAccelerationX;
        double rawAccelerationY;
        double rawAccelerationZ;
        double sampleNoiseX;
        double sampleNoiseY;
        double sampleNoiseZ;
        var correctedAccelerationX = 0.0;
        var correctedAccelerationY = 0.0;
        var correctedAccelerationZ = 0.0;

        // setup Kalman filter
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix predictedState;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix correctedState;

        // measurement (accelerationX, accelerationY, accelerationZ)
        final var measurement = new Matrix(3, 1);
        measurement.setElementAtIndex(0, accelerationX);
        measurement.setElementAtIndex(1, accelerationY);
        measurement.setElementAtIndex(2, accelerationZ);

        // assuming delta time between samples of 1 second
        final var deltaTime = 1.0;
        // transitions for position, speed and acceleration
        // [1,   deltaTime,    deltaTime^2/2,    0,  0,          0,              0,  0,          0]
        // [0,   1,            deltaTime,        0,  0,          0,              0,  0,          0]
        // [0,   0,            1,                0,  0,          0,              0,  0,          0]
        // [0,   0,            0,                1,  deltaTime,  deltaTime^2/2,  0,  0,          0]
        // [0,   0,            0,                0,  1,          deltaTime,      0,  0,          0]
        // [0,   0,            0,                0,  0,          1,              0,  0,          0]
        // [0,   0,            0,                0,  0,          0,              1,  deltaTime,  deltaTime^2/2]
        // [0,   0,            0,                0,  0,          0,              0,  1,          deltaTime]
        // [0,   0,            0,                0,  0,          0,              0,  0,          1]
        final var transitionMatrix = new Matrix(9, 9);
        // column order
        transitionMatrix.setSubmatrix(0, 0, 8, 8, new double[]
                {
                        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        deltaTime, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.5 * deltaTime * deltaTime, deltaTime, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, deltaTime, 1.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.5 * deltaTime * deltaTime, deltaTime, 1.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, deltaTime, 1.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5 * deltaTime * deltaTime, deltaTime, 1.0
                }, true);
        kalman.setTransitionMatrix(transitionMatrix);

        final var processNoiseStd = Math.sqrt(1e-6);
        final var processNoiseVariance = processNoiseStd * processNoiseStd;

        // process noise covariance is a symmetric matrix obtained as a result of:
        // processNoiseVariance * [deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1]'*[deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1]
        // which results in a symmetric matrix containing 9 times whe following block:
        // processNoiseVariance * [deltaTime^4/4,   deltaTime^3/2,  deltaTime^2/2]
        //                        [deltaTime^3/2,   deltaTime^2,    deltaTime]
        //                        [deltaTime^2/2,   deltaTime,      1]
        final var block = new Matrix(3, 3);
        // column order
        block.setSubmatrix(0, 0, 2, 2, new double[]
                {
                        processNoiseVariance * Math.pow(deltaTime, 4.0) * 0.25,
                        processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 2.0),
                        processNoiseVariance * deltaTime,
                        processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                        processNoiseVariance * deltaTime,
                        processNoiseVariance
                }, true);
        var processNoiseCov = new Matrix(9, 9);
        for (var i = 0; i < 9; i += 3) {
            for (var j = 0; j < 9; j += 3) {
                processNoiseCov.setSubmatrix(i, j, i + 2, j + 2, block);
            }
        }
        kalman.setProcessNoiseCov(processNoiseCov);

        final var measureNoiseStd = Math.sqrt(1e-1);
        final var measureNoiseVariance = measureNoiseStd * measureNoiseStd;

        // measurement noise covariance = measureNoiseVariance * I
        final var measurementNoiseCov = Matrix.identity(3, 3);
        measurementNoiseCov.multiplyByScalar(measureNoiseVariance);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final var measurementMatrix = new Matrix(3, 9);
        // column order
        measurementMatrix.setSubmatrix(0, 0, 2, 8, new double[]
                {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 1.0
                }, true);
        kalman.setMeasurementMatrix(measurementMatrix);

        if (WRITE_TO_CONSOLE) {
            System.out.println("no;truePositionX;trueSpeedX;trueAccelerationX;" +
                    "predictedPositionX;predictedSpeedX;predictedAccelerationX;" +
                    "correctedPositionX;correctedSpeedX;correctedAccelerationX;" +
                    "rawPositionX;rawSpeedX;rawAccelerationX;" +
                    "truePositionY;trueSpeedY;trueAccelerationY;" +
                    "predictedPositionY;predictedSpeedY;predictedAccelerationY;" +
                    "correctedPositionY;correctedSpeedY;correctedAccelerationY;" +
                    "rawPositionY;rawSpeedY;rawAccelerationY;" +
                    "truePositionZ;trueSpeedZ;trueAccelerationZ;" +
                    "predictedPositionZ;predictedSpeedZ;predictedAccelerationZ;" +
                    "correctedPositionZ;correctedSpeedZ;correctedAccelerationZ;" +
                    "rawPositionZ;rawSpeedZ;rawAccelerationZ;");
        }

        // assume each sample happens every second
        for (var i = 0; i < N_SAMPLES; i++) {
            trueSpeedX += accelerationX;
            trueSpeedY += accelerationY;
            trueSpeedZ += accelerationZ;
            truePositionX += trueSpeedX;
            truePositionY += trueSpeedY;
            truePositionZ += trueSpeedZ;

            sampleNoiseX = measureNoiseStd * rand.nextGaussian();
            sampleNoiseY = measureNoiseStd * rand.nextGaussian();
            sampleNoiseZ = measureNoiseStd * rand.nextGaussian();

            rawAccelerationX = accelerationX + sampleNoiseX;
            rawAccelerationY = accelerationY + sampleNoiseY;
            rawAccelerationZ = accelerationZ + sampleNoiseZ;
            rawSpeedX += rawAccelerationX;
            rawSpeedY += rawAccelerationY;
            rawSpeedZ += rawAccelerationZ;
            rawPositionX += rawSpeedX;
            rawPositionY += rawSpeedY;
            rawPositionZ += rawSpeedZ;

            predictedState = kalman.predict();
            assertNotNull(predictedState);

            // correct state for the 75% of samples if  DO_NOT_SKIP is enabled
            if (rand.nextDouble() <= 0.75 || DO_NOT_SKIP) {
                // take measure with noise
                measurement.setElementAtIndex(0, rawAccelerationX);
                measurement.setElementAtIndex(1, rawAccelerationY);
                measurement.setElementAtIndex(2, rawAccelerationZ);

                // correct
                correctedState = kalman.correct(measurement);

                correctedAccelerationX = correctedState.getElementAtIndex(2);
                correctedAccelerationY = correctedState.getElementAtIndex(5);
                correctedAccelerationZ = correctedState.getElementAtIndex(8);

                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";" + truePositionX + ";" + trueSpeedX + ";" + accelerationX + ";"
                            + correctedAccelerationX + ";" + rawPositionX + ";" + rawSpeedX + ";"
                            + rawAccelerationX + ";" + truePositionY + ";" + trueSpeedY + ";" + accelerationY + ";"
                            + correctedAccelerationY + ";" + rawPositionY + ";" + rawSpeedY + ";"
                            + rawAccelerationY + ";" + truePositionZ + ";" + trueSpeedZ + ";" + accelerationZ + ";"
                            + correctedAccelerationZ + ";" + rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
                }

            } else {
                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";" + truePositionX + ";" + trueSpeedX + ";" + accelerationX + ";" + ";;;"
                            + rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" + truePositionY + ";"
                            + trueSpeedY + ";" + accelerationY + ";" + ";;;" + rawPositionY + ";" + rawSpeedY + ";"
                            + rawAccelerationY + ";" + truePositionZ + ";" + trueSpeedZ + ";" + accelerationZ + ";"
                            + ";;;" + rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
                }
            }
        }

        // check that the filter state converged to the true acceleration
        assertEquals(accelerationX, correctedAccelerationX, measureNoiseVariance);
        assertEquals(accelerationY, correctedAccelerationY, measureNoiseVariance);
        assertEquals(accelerationZ, correctedAccelerationZ, measureNoiseVariance);
    }

    @Test
    void testPredictAndCorrectRealDataNoMotion() throws IOException, SignalProcessingException, WrongSizeException {

        var f = new File(ACCELERATION_NO_MOTION);
        final var data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final var kalman = new KalmanFilter(9, 3);

        var rawSpeedX = 0.0;
        var rawSpeedY = 0.0;
        var rawSpeedZ = 0.0;
        var rawPositionX = 0.0;
        var rawPositionY = 0.0;
        var rawPositionZ = 0.0;
        double rawAccelerationX;
        double rawAccelerationY;
        double rawAccelerationZ;
        double predictedPositionX;
        double predictedPositionY;
        double predictedPositionZ;
        double predictedSpeedX;
        double predictedSpeedY;
        double predictedSpeedZ;
        double predictedAccelerationX;
        double predictedAccelerationY;
        double predictedAccelerationZ;
        double correctedPositionX;
        double correctedPositionY;
        double correctedPositionZ;
        double correctedSpeedX;
        double correctedSpeedY;
        double correctedSpeedZ;
        double correctedAccelerationX;
        double correctedAccelerationY;
        double correctedAccelerationZ;

        // setup Kalman filter
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix predictedState;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix correctedState;

        // measurement (accelerationX, accelerationY, accelerationZ)
        final var measurement = new Matrix(3, 1);

        final var transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final var processNoiseStd = Math.sqrt(1e-6);
        final var processNoiseVariance = processNoiseStd * processNoiseStd;

        final var block = new Matrix(3, 3);
        final var processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block, processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final var measurementNoiseCov = noiseCovarianceMatrix(data);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final var measurementMatrix = new Matrix(3, 9);
        // column order
        measurementMatrix.setSubmatrix(0, 0, 2, 8, new double[]
                {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 1.0
                }, true);
        kalman.setMeasurementMatrix(measurementMatrix);

        f = new File(CSV_FILE_NO_MOTION);
        final var stream = new PrintStream(f);

        stream.println("time;truePositionX;trueSpeedX;trueAccelerationX;" +
                "predictedPositionX;predictedSpeedX;predictedAccelerationX;" +
                "correctedPositionX;correctedSpeedX;correctedAccelerationX;" +
                "rawPositionX;rawSpeedX;rawAccelerationX;" +
                "truePositionY;trueSpeedY;trueAccelerationY;" +
                "predictedPositionY;predictedSpeedY;predictedAccelerationY;" +
                "correctedPositionY;correctedSpeedY;correctedAccelerationY;" +
                "rawPositionY;rawSpeedY;rawAccelerationY;" +
                "truePositionZ;trueSpeedZ;trueAccelerationZ;" +
                "predictedPositionZ;predictedSpeedZ;predictedAccelerationZ;" +
                "correctedPositionZ;correctedSpeedZ;correctedAccelerationZ;" +
                "rawPositionZ;rawSpeedZ;rawAccelerationZ;");

        var deltaTime = 20e-3;
        long prevTimestampNanos;
        var timestampNanos = 0L;
        var time = 0.0;
        for (var i = 0; i < data.numSamples; i++) {
            prevTimestampNanos = timestampNanos;
            timestampNanos = data.timestamp[i];
            if (i > 0) {
                deltaTime = (timestampNanos - prevTimestampNanos) * 1e-9;
                time += deltaTime;
            }

            rawAccelerationX = data.accelerationX[i];
            rawAccelerationY = data.accelerationY[i];
            rawAccelerationZ = data.accelerationZ[i];
            rawSpeedX += rawAccelerationX;
            rawSpeedY += rawAccelerationY;
            rawSpeedZ += rawAccelerationZ;
            rawPositionX += rawSpeedX;
            rawPositionY += rawSpeedY;
            rawPositionZ += rawSpeedZ;

            updateTransitionMatrix(deltaTime, transitionMatrix);
            kalman.setTransitionMatrix(transitionMatrix);

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block, processNoiseCov);
            kalman.setProcessNoiseCov(processNoiseCov);

            predictedState = kalman.predict();

            predictedPositionX = predictedState.getElementAtIndex(0);
            predictedSpeedX = predictedState.getElementAtIndex(1);
            predictedAccelerationX = predictedState.getElementAtIndex(2);

            predictedPositionY = predictedState.getElementAtIndex(3);
            predictedSpeedY = predictedState.getElementAtIndex(4);
            predictedAccelerationY = predictedState.getElementAtIndex(5);

            predictedPositionZ = predictedState.getElementAtIndex(6);
            predictedSpeedZ = predictedState.getElementAtIndex(7);
            predictedAccelerationZ = predictedState.getElementAtIndex(8);

            // correct system using new acceleration measures
            measurement.setElementAtIndex(0, rawAccelerationX);
            measurement.setElementAtIndex(1, rawAccelerationY);
            measurement.setElementAtIndex(2, rawAccelerationZ);

            // correct
            correctedState = kalman.correct(measurement);

            correctedPositionX = correctedState.getElementAtIndex(0);
            correctedSpeedX = correctedState.getElementAtIndex(1);
            correctedAccelerationX = correctedState.getElementAtIndex(2);

            correctedPositionY = correctedState.getElementAtIndex(3);
            correctedSpeedY = correctedState.getElementAtIndex(4);
            correctedAccelerationY = correctedState.getElementAtIndex(5);

            correctedPositionZ = correctedState.getElementAtIndex(6);
            correctedSpeedZ = correctedState.getElementAtIndex(7);
            correctedAccelerationZ = correctedState.getElementAtIndex(8);

            stream.println(time + ";;;;" +
                    predictedPositionX + ";" + predictedSpeedX + ";" + predictedAccelerationX + ";" +
                    correctedPositionX + ";" + correctedSpeedX + ";" + correctedAccelerationX + ";" +
                    rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" +
                    ";;;" +
                    predictedPositionY + ";" + predictedSpeedY + ";" + predictedAccelerationY + ";" +
                    correctedPositionY + ";" + correctedSpeedY + ";" + correctedAccelerationY + ";" +
                    rawPositionY + ";" + rawSpeedY + ";" + rawAccelerationY + ";" +
                    ";;;" +
                    predictedPositionZ + ";" + predictedSpeedZ + ";" + predictedAccelerationZ + ";" +
                    correctedPositionZ + ";" + correctedSpeedZ + ";" + correctedAccelerationZ + ";" +
                    rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
        }

        stream.close();

        assertTrue(f.exists());
    }

    @Test
    void testPredictAndCorrectRealDataMotion() throws IOException, SignalProcessingException, WrongSizeException {

        var f = new File(ACCELERATION_NO_MOTION);
        final var dataNoMotion = AccelerationFileLoader.load(f);

        f = new File(ACCELERATION_MOTION);
        final var data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final var kalman = new KalmanFilter(9, 3);

        var rawSpeedX = 0.0;
        var rawSpeedY = 0.0;
        var rawSpeedZ = 0.0;
        var rawPositionX = 0.0;
        var rawPositionY = 0.0;
        var rawPositionZ = 0.0;
        double rawAccelerationX;
        double rawAccelerationY;
        double rawAccelerationZ;
        double predictedPositionX;
        double predictedPositionY;
        double predictedPositionZ;
        double predictedSpeedX;
        double predictedSpeedY;
        double predictedSpeedZ;
        double predictedAccelerationX;
        double predictedAccelerationY;
        double predictedAccelerationZ;
        double correctedPositionX;
        double correctedPositionY;
        double correctedPositionZ;
        double correctedSpeedX;
        double correctedSpeedY;
        double correctedSpeedZ;
        double correctedAccelerationX;
        double correctedAccelerationY;
        double correctedAccelerationZ;

        // setup Kalman filter
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix predictedState;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix correctedState;

        // measurement (accelerationX, accelerationY, accelerationZ)
        final var measurement = new Matrix(3, 1);

        final var transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final var processNoiseStd = Math.sqrt(1e-6);
        final var processNoiseVariance = processNoiseStd * processNoiseStd;

        final var block = new Matrix(3, 3);
        final var processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block, processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final var measurementNoiseCov = noiseCovarianceMatrix(dataNoMotion);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final var measurementMatrix = new Matrix(3, 9);
        // column order
        measurementMatrix.setSubmatrix(0, 0, 2, 8, new double[]
                {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 1.0
                }, true);
        kalman.setMeasurementMatrix(measurementMatrix);

        f = new File(CSV_FILE_MOTION);
        final var stream = new PrintStream(f);

        stream.println("time;truePositionX;trueSpeedX;trueAccelerationX;" +
                "predictedPositionX;predictedSpeedX;predictedAccelerationX;" +
                "correctedPositionX;correctedSpeedX;correctedAccelerationX;" +
                "rawPositionX;rawSpeedX;rawAccelerationX;" +
                "truePositionY;trueSpeedY;trueAccelerationY;" +
                "predictedPositionY;predictedSpeedY;predictedAccelerationY;" +
                "correctedPositionY;correctedSpeedY;correctedAccelerationY;" +
                "rawPositionY;rawSpeedY;rawAccelerationY;" +
                "truePositionZ;trueSpeedZ;trueAccelerationZ;" +
                "predictedPositionZ;predictedSpeedZ;predictedAccelerationZ;" +
                "correctedPositionZ;correctedSpeedZ;correctedAccelerationZ;" +
                "rawPositionZ;rawSpeedZ;rawAccelerationZ;");

        var deltaTime = 20e-3;
        long prevTimestampNanos;
        var timestampNanos = 0L;
        var time = 0.0;
        for (var i = 0; i < data.numSamples; i++) {
            prevTimestampNanos = timestampNanos;
            timestampNanos = data.timestamp[i];
            if (i > 0) {
                deltaTime = (timestampNanos - prevTimestampNanos) * 1e-9;
                time += deltaTime;
            }

            rawAccelerationX = data.accelerationX[i];
            rawAccelerationY = data.accelerationY[i];
            rawAccelerationZ = data.accelerationZ[i];
            rawSpeedX += rawAccelerationX;
            rawSpeedY += rawAccelerationY;
            rawSpeedZ += rawAccelerationZ;
            rawPositionX += rawSpeedX;
            rawPositionY += rawSpeedY;
            rawPositionZ += rawSpeedZ;

            updateTransitionMatrix(deltaTime, transitionMatrix);
            kalman.setTransitionMatrix(transitionMatrix);

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block, processNoiseCov);
            kalman.setProcessNoiseCov(processNoiseCov);

            predictedState = kalman.predict();

            predictedPositionX = predictedState.getElementAtIndex(0);
            predictedSpeedX = predictedState.getElementAtIndex(1);
            predictedAccelerationX = predictedState.getElementAtIndex(2);

            predictedPositionY = predictedState.getElementAtIndex(3);
            predictedSpeedY = predictedState.getElementAtIndex(4);
            predictedAccelerationY = predictedState.getElementAtIndex(5);

            predictedPositionZ = predictedState.getElementAtIndex(6);
            predictedSpeedZ = predictedState.getElementAtIndex(7);
            predictedAccelerationZ = predictedState.getElementAtIndex(8);

            // correct system using new acceleration measures
            measurement.setElementAtIndex(0, rawAccelerationX);
            measurement.setElementAtIndex(1, rawAccelerationY);
            measurement.setElementAtIndex(2, rawAccelerationZ);

            // correct
            correctedState = kalman.correct(measurement);

            correctedPositionX = correctedState.getElementAtIndex(0);
            correctedSpeedX = correctedState.getElementAtIndex(1);
            correctedAccelerationX = correctedState.getElementAtIndex(2);

            correctedPositionY = correctedState.getElementAtIndex(3);
            correctedSpeedY = correctedState.getElementAtIndex(4);
            correctedAccelerationY = correctedState.getElementAtIndex(5);

            correctedPositionZ = correctedState.getElementAtIndex(6);
            correctedSpeedZ = correctedState.getElementAtIndex(7);
            correctedAccelerationZ = correctedState.getElementAtIndex(8);

            stream.println(time + ";;;;" +
                    predictedPositionX + ";" + predictedSpeedX + ";" + predictedAccelerationX + ";" +
                    correctedPositionX + ";" + correctedSpeedX + ";" + correctedAccelerationX + ";" +
                    rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" +
                    ";;;" +
                    predictedPositionY + ";" + predictedSpeedY + ";" + predictedAccelerationY + ";" +
                    correctedPositionY + ";" + correctedSpeedY + ";" + correctedAccelerationY + ";" +
                    rawPositionY + ";" + rawSpeedY + ";" + rawAccelerationY + ";" +
                    ";;;" +
                    predictedPositionZ + ";" + predictedSpeedZ + ";" + predictedAccelerationZ + ";" +
                    correctedPositionZ + ";" + correctedSpeedZ + ";" + correctedAccelerationZ + ";" +
                    rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
        }

        stream.close();

        assertTrue(f.exists());
    }

    @Test
    void testPredictAndCorrectRealDataNoMotionFast() throws IOException, SignalProcessingException, WrongSizeException {

        var f = new File(ACCELERATION_NO_MOTION_FAST);
        final var data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final var kalman = new KalmanFilter(9, 3);

        var rawSpeedX = 0.0;
        var rawSpeedY = 0.0;
        var rawSpeedZ = 0.0;
        var rawPositionX = 0.0;
        var rawPositionY = 0.0;
        var rawPositionZ = 0.0;
        double rawAccelerationX;
        double rawAccelerationY;
        double rawAccelerationZ;
        double predictedPositionX;
        double predictedPositionY;
        double predictedPositionZ;
        double predictedSpeedX;
        double predictedSpeedY;
        double predictedSpeedZ;
        double predictedAccelerationX;
        double predictedAccelerationY;
        double predictedAccelerationZ;
        double correctedPositionX;
        double correctedPositionY;
        double correctedPositionZ;
        double correctedSpeedX;
        double correctedSpeedY;
        double correctedSpeedZ;
        double correctedAccelerationX;
        double correctedAccelerationY;
        double correctedAccelerationZ;

        // setup Kalman filter
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix predictedState;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix correctedState;

        // measurement (accelerationX, accelerationY, accelerationZ)
        final var measurement = new Matrix(3, 1);

        final var transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final var processNoiseStd = Math.sqrt(1e-6);
        final var processNoiseVariance = processNoiseStd * processNoiseStd;

        final var block = new Matrix(3, 3);
        final var processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block, processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final var measurementNoiseCov = noiseCovarianceMatrix(data);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final var measurementMatrix = new Matrix(3, 9);
        // column order
        measurementMatrix.setSubmatrix(0, 0, 2, 8, new double[]
                {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 1.0
                }, true);
        kalman.setMeasurementMatrix(measurementMatrix);

        f = new File(CSV_FILE_NO_MOTION_FAST);
        final var stream = new PrintStream(f);

        stream.println("time;truePositionX;trueSpeedX;trueAccelerationX;" +
                "predictedPositionX;predictedSpeedX;predictedAccelerationX;" +
                "correctedPositionX;correctedSpeedX;correctedAccelerationX;" +
                "rawPositionX;rawSpeedX;rawAccelerationX;" +
                "truePositionY;trueSpeedY;trueAccelerationY;" +
                "predictedPositionY;predictedSpeedY;predictedAccelerationY;" +
                "correctedPositionY;correctedSpeedY;correctedAccelerationY;" +
                "rawPositionY;rawSpeedY;rawAccelerationY;" +
                "truePositionZ;trueSpeedZ;trueAccelerationZ;" +
                "predictedPositionZ;predictedSpeedZ;predictedAccelerationZ;" +
                "correctedPositionZ;correctedSpeedZ;correctedAccelerationZ;" +
                "rawPositionZ;rawSpeedZ;rawAccelerationZ;");

        var deltaTime = 4e-3;
        long prevTimestampNanos;
        var timestampNanos = 0L;
        var time = 0.0;
        for (var i = 0; i < data.numSamples; i++) {
            prevTimestampNanos = timestampNanos;
            timestampNanos = data.timestamp[i];
            if (i > 0) {
                deltaTime = (timestampNanos - prevTimestampNanos) * 1e-9;
                time += deltaTime;
            }

            rawAccelerationX = data.accelerationX[i];
            rawAccelerationY = data.accelerationY[i];
            rawAccelerationZ = data.accelerationZ[i];
            rawSpeedX += rawAccelerationX;
            rawSpeedY += rawAccelerationY;
            rawSpeedZ += rawAccelerationZ;
            rawPositionX += rawSpeedX;
            rawPositionY += rawSpeedY;
            rawPositionZ += rawSpeedZ;

            updateTransitionMatrix(deltaTime, transitionMatrix);
            kalman.setTransitionMatrix(transitionMatrix);

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block, processNoiseCov);
            kalman.setProcessNoiseCov(processNoiseCov);

            predictedState = kalman.predict();

            predictedPositionX = predictedState.getElementAtIndex(0);
            predictedSpeedX = predictedState.getElementAtIndex(1);
            predictedAccelerationX = predictedState.getElementAtIndex(2);

            predictedPositionY = predictedState.getElementAtIndex(3);
            predictedSpeedY = predictedState.getElementAtIndex(4);
            predictedAccelerationY = predictedState.getElementAtIndex(5);

            predictedPositionZ = predictedState.getElementAtIndex(6);
            predictedSpeedZ = predictedState.getElementAtIndex(7);
            predictedAccelerationZ = predictedState.getElementAtIndex(8);

            // correct system using new acceleration measures
            measurement.setElementAtIndex(0, rawAccelerationX);
            measurement.setElementAtIndex(1, rawAccelerationY);
            measurement.setElementAtIndex(2, rawAccelerationZ);

            // correct
            correctedState = kalman.correct(measurement);

            correctedPositionX = correctedState.getElementAtIndex(0);
            correctedSpeedX = correctedState.getElementAtIndex(1);
            correctedAccelerationX = correctedState.getElementAtIndex(2);

            correctedPositionY = correctedState.getElementAtIndex(3);
            correctedSpeedY = correctedState.getElementAtIndex(4);
            correctedAccelerationY = correctedState.getElementAtIndex(5);

            correctedPositionZ = correctedState.getElementAtIndex(6);
            correctedSpeedZ = correctedState.getElementAtIndex(7);
            correctedAccelerationZ = correctedState.getElementAtIndex(8);

            stream.println(time + ";;;;" +
                    predictedPositionX + ";" + predictedSpeedX + ";" + predictedAccelerationX + ";" +
                    correctedPositionX + ";" + correctedSpeedX + ";" + correctedAccelerationX + ";" +
                    rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" +
                    ";;;" +
                    predictedPositionY + ";" + predictedSpeedY + ";" + predictedAccelerationY + ";" +
                    correctedPositionY + ";" + correctedSpeedY + ";" + correctedAccelerationY + ";" +
                    rawPositionY + ";" + rawSpeedY + ";" + rawAccelerationY + ";" +
                    ";;;" +
                    predictedPositionZ + ";" + predictedSpeedZ + ";" + predictedAccelerationZ + ";" +
                    correctedPositionZ + ";" + correctedSpeedZ + ";" + correctedAccelerationZ + ";" +
                    rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
        }

        stream.close();

        assertTrue(f.exists());
    }

    @Test
    void testPredictAndCorrectRealDataMotionFast() throws IOException, SignalProcessingException, WrongSizeException {

        var f = new File(ACCELERATION_NO_MOTION_FAST);
        final var dataNoMotion = AccelerationFileLoader.load(f);

        f = new File(ACCELERATION_MOTION_FAST);
        final var data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final var kalman = new KalmanFilter(9, 3);

        var rawSpeedX = 0.0;
        var rawSpeedY = 0.0;
        var rawSpeedZ = 0.0;
        var rawPositionX = 0.0;
        var rawPositionY = 0.0;
        var rawPositionZ = 0.0;
        double rawAccelerationX;
        double rawAccelerationY;
        double rawAccelerationZ;
        double predictedPositionX;
        double predictedPositionY;
        double predictedPositionZ;
        double predictedSpeedX;
        double predictedSpeedY;
        double predictedSpeedZ;
        double predictedAccelerationX;
        double predictedAccelerationY;
        double predictedAccelerationZ;
        double correctedPositionX;
        double correctedPositionY;
        double correctedPositionZ;
        double correctedSpeedX;
        double correctedSpeedY;
        double correctedSpeedZ;
        double correctedAccelerationX;
        double correctedAccelerationY;
        double correctedAccelerationZ;

        // setup Kalman filter
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix predictedState;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix correctedState;

        // measurement (accelerationX, accelerationY, accelerationZ)
        final var measurement = new Matrix(3, 1);

        final var transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final var processNoiseStd = Math.sqrt(1e-6);
        final var processNoiseVariance = processNoiseStd * processNoiseStd;

        final var block = new Matrix(3, 3);
        final var processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block, processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final var measurementNoiseCov = noiseCovarianceMatrix(dataNoMotion);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final var measurementMatrix = new Matrix(3, 9);
        // column order
        measurementMatrix.setSubmatrix(0, 0, 2, 8, new double[]
                {
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 1.0
                }, true);
        kalman.setMeasurementMatrix(measurementMatrix);

        f = new File(CSV_FILE_MOTION_FAST);
        final var stream = new PrintStream(f);

        stream.println("time;truePositionX;trueSpeedX;trueAccelerationX;" +
                "predictedPositionX;predictedSpeedX;predictedAccelerationX;" +
                "correctedPositionX;correctedSpeedX;correctedAccelerationX;" +
                "rawPositionX;rawSpeedX;rawAccelerationX;" +
                "truePositionY;trueSpeedY;trueAccelerationY;" +
                "predictedPositionY;predictedSpeedY;predictedAccelerationY;" +
                "correctedPositionY;correctedSpeedY;correctedAccelerationY;" +
                "rawPositionY;rawSpeedY;rawAccelerationY;" +
                "truePositionZ;trueSpeedZ;trueAccelerationZ;" +
                "predictedPositionZ;predictedSpeedZ;predictedAccelerationZ;" +
                "correctedPositionZ;correctedSpeedZ;correctedAccelerationZ;" +
                "rawPositionZ;rawSpeedZ;rawAccelerationZ;");

        var deltaTime = 4e-3;
        long prevTimestampNanos;
        var timestampNanos = 0L;
        var time = 0.0;
        for (var i = 0; i < data.numSamples; i++) {
            prevTimestampNanos = timestampNanos;
            timestampNanos = data.timestamp[i];
            if (i > 0) {
                deltaTime = (timestampNanos - prevTimestampNanos) * 1e-9;
                time += deltaTime;
            }

            rawAccelerationX = data.accelerationX[i];
            rawAccelerationY = data.accelerationY[i];
            rawAccelerationZ = data.accelerationZ[i];
            rawSpeedX += rawAccelerationX;
            rawSpeedY += rawAccelerationY;
            rawSpeedZ += rawAccelerationZ;
            rawPositionX += rawSpeedX;
            rawPositionY += rawSpeedY;
            rawPositionZ += rawSpeedZ;

            updateTransitionMatrix(deltaTime, transitionMatrix);
            kalman.setTransitionMatrix(transitionMatrix);

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block, processNoiseCov);
            kalman.setProcessNoiseCov(processNoiseCov);

            predictedState = kalman.predict();

            predictedPositionX = predictedState.getElementAtIndex(0);
            predictedSpeedX = predictedState.getElementAtIndex(1);
            predictedAccelerationX = predictedState.getElementAtIndex(2);

            predictedPositionY = predictedState.getElementAtIndex(3);
            predictedSpeedY = predictedState.getElementAtIndex(4);
            predictedAccelerationY = predictedState.getElementAtIndex(5);

            predictedPositionZ = predictedState.getElementAtIndex(6);
            predictedSpeedZ = predictedState.getElementAtIndex(7);
            predictedAccelerationZ = predictedState.getElementAtIndex(8);

            // correct system using new acceleration measures
            measurement.setElementAtIndex(0, rawAccelerationX);
            measurement.setElementAtIndex(1, rawAccelerationY);
            measurement.setElementAtIndex(2, rawAccelerationZ);

            // correct
            correctedState = kalman.correct(measurement);

            correctedPositionX = correctedState.getElementAtIndex(0);
            correctedSpeedX = correctedState.getElementAtIndex(1);
            correctedAccelerationX = correctedState.getElementAtIndex(2);

            correctedPositionY = correctedState.getElementAtIndex(3);
            correctedSpeedY = correctedState.getElementAtIndex(4);
            correctedAccelerationY = correctedState.getElementAtIndex(5);

            correctedPositionZ = correctedState.getElementAtIndex(6);
            correctedSpeedZ = correctedState.getElementAtIndex(7);
            correctedAccelerationZ = correctedState.getElementAtIndex(8);

            stream.println(time + ";;;;" +
                    predictedPositionX + ";" + predictedSpeedX + ";" + predictedAccelerationX + ";" +
                    correctedPositionX + ";" + correctedSpeedX + ";" + correctedAccelerationX + ";" +
                    rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" +
                    ";;;" +
                    predictedPositionY + ";" + predictedSpeedY + ";" + predictedAccelerationY + ";" +
                    correctedPositionY + ";" + correctedSpeedY + ";" + correctedAccelerationY + ";" +
                    rawPositionY + ";" + rawSpeedY + ";" + rawAccelerationY + ";" +
                    ";;;" +
                    predictedPositionZ + ";" + predictedSpeedZ + ";" + predictedAccelerationZ + ";" +
                    correctedPositionZ + ";" + correctedSpeedZ + ";" + correctedAccelerationZ + ";" +
                    rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
        }

        stream.close();

        assertTrue(f.exists());
    }

    @Test
    void testGetSetStatePre() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        // check default value
        final var statePre = filter.getStatePre();

        assertEquals(6, statePre.getRows());
        assertEquals(1, statePre.getColumns());

        // set new value
        final var statePre2 = new Matrix(6, 1);
        filter.setStatePre(statePre2);

        // check correctness
        assertSame(statePre2, filter.getStatePre());

        // Force IllegalArgumentException
        final var wrongStatePre = new Matrix(5, 1);
        assertThrows(IllegalArgumentException.class, () -> filter.setStatePre(wrongStatePre));

        final var wrongStatePre2 = new Matrix(6, 2);
        assertThrows(IllegalArgumentException.class, () -> filter.setStatePre(wrongStatePre2));
    }

    @Test
    void testGetSetStatePost() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        // check default value
        final var statePost = filter.getStatePost();

        assertEquals(6, statePost.getRows());
        assertEquals(1, statePost.getColumns());

        // set new value
        final var statePost2 = new Matrix(6, 1);
        filter.setStatePost(statePost2);

        // check correctness
        assertSame(statePost2, filter.getStatePost());

        // Force IllegalArgumentException
        final var wrongStatePost = new Matrix(5, 1);
        assertThrows(IllegalArgumentException.class, () -> filter.setStatePost(wrongStatePost));

        final var wrongStatePost2 = new Matrix(6, 2);
        assertThrows(IllegalArgumentException.class, () -> filter.setStatePost(wrongStatePost2));
    }

    @Test
    void testGetSetTransitionMatrix() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        // check default value
        final var transitionMatrix = filter.getTransitionMatrix();

        assertEquals(6, transitionMatrix.getRows());
        assertEquals(6, transitionMatrix.getColumns());

        // set new value
        final var transitionMatrix2 = new Matrix(6, 6);
        filter.setTransitionMatrix(transitionMatrix2);

        // check correctness
        assertSame(transitionMatrix2, filter.getTransitionMatrix());

        // Force IllegalArgumentException
        final var wrongTransitionMatrix = new Matrix(5, 6);
        assertThrows(IllegalArgumentException.class, () -> filter.setTransitionMatrix(wrongTransitionMatrix));

        final var wrongTransitionMatrix2 = new Matrix(6, 5);
        assertThrows(IllegalArgumentException.class, () -> filter.setTransitionMatrix(wrongTransitionMatrix2));
    }

    @Test
    void testGetSetControlMatrix() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9, 1);

        // check default value
        final var controlMatrix = filter.getControlMatrix();

        assertEquals(6, controlMatrix.getRows());
        assertEquals(1, controlMatrix.getColumns());

        // set new value
        final var controlMatrix2 = new Matrix(6, 1);
        filter.setControlMatrix(controlMatrix2);

        // check correctness
        assertSame(controlMatrix2, filter.getControlMatrix());

        // Force IllegalArgumentException
        final var wrongControlMatrix = new Matrix(5, 1);
        assertThrows(IllegalArgumentException.class, () -> filter.setControlMatrix(wrongControlMatrix));

        final var wrongControlMatrix2 = new Matrix(6, 2);
        assertThrows(IllegalArgumentException.class, () -> filter.setControlMatrix(wrongControlMatrix2));
    }

    @Test
    void testGetSetMeasurementMatrix() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        // check default value
        final var measurementMatrix = filter.getMeasurementMatrix();

        assertEquals(9, measurementMatrix.getRows());
        assertEquals(6, measurementMatrix.getColumns());

        // set new value
        final var measurementMatrix2 = new Matrix(9, 6);
        filter.setMeasurementMatrix(measurementMatrix2);

        // check correctness
        assertSame(measurementMatrix2, filter.getMeasurementMatrix());

        // Force IllegalArgumentException
        final var wrongMeasurementMatrix = new Matrix(8, 6);
        assertThrows(IllegalArgumentException.class, () -> filter.setMeasurementMatrix(wrongMeasurementMatrix));

        final var wrongMeasurementMatrix2 = new Matrix(9, 5);
        assertThrows(IllegalArgumentException.class, () -> filter.setMeasurementMatrix(wrongMeasurementMatrix2));
    }

    @Test
    void testGetSetProcessNoiseCov() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        // check default value
        final var processNoiseCov = filter.getProcessNoiseCov();

        assertEquals(6, processNoiseCov.getRows());
        assertEquals(6, processNoiseCov.getColumns());

        // set new value
        final var processNoiseCov2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        filter.setProcessNoiseCov(processNoiseCov2);

        // check correctness
        assertSame(processNoiseCov2, filter.getProcessNoiseCov());

        // Force IllegalArgumentException
        final var wrongProcessNoiseCov = new Matrix(5, 6);
        assertThrows(IllegalArgumentException.class, () -> filter.setProcessNoiseCov(wrongProcessNoiseCov));

        final var wrongProcessNoiseCov2 = new Matrix(6, 5);
        assertThrows(IllegalArgumentException.class, () -> filter.setProcessNoiseCov(wrongProcessNoiseCov2));

        final var wrongProcessNoiseCov3 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        wrongProcessNoiseCov3.setElementAt(0, 1, 1.0);
        assertThrows(IllegalArgumentException.class, () -> filter.setProcessNoiseCov(wrongProcessNoiseCov3));
    }

    @Test
    void testGetSetMeasurementNoiseCov() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        // check default value
        final var measurementNoiseCov = filter.getMeasurementNoiseCov();

        assertEquals(9, measurementNoiseCov.getRows());
        assertEquals(9, measurementNoiseCov.getColumns());

        // set new value
        final var measurementNoiseCov2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9});
        filter.setMeasurementNoiseCov(measurementNoiseCov2);

        // check correctness
        assertSame(measurementNoiseCov2, filter.getMeasurementNoiseCov());

        // Force IllegalArgumentException
        final var wrongMeasurementNoiseCov = new Matrix(8, 9);
        assertThrows(IllegalArgumentException.class, () -> filter.setMeasurementNoiseCov(wrongMeasurementNoiseCov));

        final var wrongMeasurementNoiseCov2 = new Matrix(9, 8);
        assertThrows(IllegalArgumentException.class, () -> filter.setMeasurementNoiseCov(wrongMeasurementNoiseCov2));

        final var wrongMeasurementNoiseCov3 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9});
        wrongMeasurementNoiseCov3.setElementAt(0, 1, 1.0);
        assertThrows(IllegalArgumentException.class, () -> filter.setMeasurementNoiseCov(wrongMeasurementNoiseCov3));
    }

    @Test
    void testGetSetErrorCovPre() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        final var errorCovPre = filter.getErrorCovPre();

        assertEquals(6, errorCovPre.getRows());
        assertEquals(6, errorCovPre.getColumns());

        // set new value
        final var errorCovPre2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        filter.setErrorCovPre(errorCovPre2);

        // check correctness
        assertSame(errorCovPre2, filter.getErrorCovPre());

        // Force IllegalArgumentException
        final var wrongErrorCovPre = new Matrix(5, 6);
        assertThrows(IllegalArgumentException.class, () -> filter.setErrorCovPre(wrongErrorCovPre));

        final var wrongErrorCovPre2 = new Matrix(6, 5);
        assertThrows(IllegalArgumentException.class, () -> filter.setErrorCovPre(wrongErrorCovPre2));

        final var wrongErrorCovPre3 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        wrongErrorCovPre3.setElementAt(0, 1, 1.0);
        assertThrows(IllegalArgumentException.class, () -> filter.setErrorCovPre(wrongErrorCovPre3));
    }

    @Test
    void testGetSetGain() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        final var gain = filter.getGain();

        assertEquals(6, gain.getRows());
        assertEquals(9, gain.getColumns());

        // set new value
        final var gain2 = new Matrix(6, 9);
        filter.setGain(gain2);

        // check correctness
        assertSame(gain2, filter.getGain());

        // Force IllegalArgumentException
        final var wrongGain = new Matrix(5, 9);
        assertThrows(IllegalArgumentException.class, () -> filter.setGain(wrongGain));

        final var wrongGain2 = new Matrix(6, 8);
        assertThrows(IllegalArgumentException.class, () -> filter.setGain(wrongGain2));
    }

    @Test
    void testGetSetErrorCovPost() throws SignalProcessingException, WrongSizeException {
        final var filter = new KalmanFilter(6, 9);

        final var errorCovPost = filter.getErrorCovPost();

        assertEquals(6, errorCovPost.getRows());
        assertEquals(6, errorCovPost.getColumns());

        // set new value
        final var errorCovPost2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        filter.setErrorCovPost(errorCovPost2);

        // check correctness
        assertSame(errorCovPost2, filter.getErrorCovPost());

        // Force IllegalArgumentException
        final var wrongErrorCovPost = new Matrix(5, 6);
        assertThrows(IllegalArgumentException.class, () -> filter.setErrorCovPost(wrongErrorCovPost));

        final var wrongErrorCovPost2 = new Matrix(6, 5);
        assertThrows(IllegalArgumentException.class, () -> filter.setErrorCovPost(wrongErrorCovPost2));

        final var wrongErrorCovPost3 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        wrongErrorCovPost3.setElementAt(0, 1, 1.0);
        assertThrows(IllegalArgumentException.class, () -> filter.setErrorCovPost(wrongErrorCovPost3));
    }

    private static void updateProcessNoiseCov(
            final double processNoiseVariance, final double deltaTime, final Matrix block, final Matrix result) {

        // block matrix has the following form:
        // processNoiseVariance * [deltaTime^4/4,   deltaTime^3/2,  deltaTime^2/2]
        //                        [deltaTime^3/2,   deltaTime^2,    deltaTime]
        //                        [deltaTime^2/2,   deltaTime,      1]

        // set elements in column order
        block.setElementAtIndex(0, processNoiseVariance * Math.pow(deltaTime, 4.0) * 0.25,
                true);
        block.setElementAtIndex(1, processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5,
                true);
        block.setElementAtIndex(2, processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                true);

        block.setElementAtIndex(3, processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5,
                true);
        block.setElementAtIndex(4, processNoiseVariance * Math.pow(deltaTime, 2.0), true);
        block.setElementAtIndex(5, processNoiseVariance * deltaTime, true);

        block.setElementAtIndex(6, processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                true);
        block.setElementAtIndex(7, processNoiseVariance * deltaTime, true);
        block.setElementAtIndex(8, processNoiseVariance);

        // repeat block matrix throughout result matrix (9 times)
        for (var i = 0; i < 9; i += 3) {
            for (var j = 0; j < 9; j += 3) {
                result.setSubmatrix(i, j, i + 2, j + 2, block);
            }
        }
    }

    private static void updateTransitionMatrix(final double deltaTime, final Matrix result) {
        //transitions for position, speed and acceleration
        //[1,   deltaTime,    deltaTime^2/2,    0,  0,          0,              0,  0,          0]
        //[0,   1,            deltaTime,        0,  0,          0,              0,  0,          0]
        //[0,   0,            1,                0,  0,          0,              0,  0,          0]
        //[0,   0,            0,                1,  deltaTime,  deltaTime^2/2,  0,  0,          0]
        //[0,   0,            0,                0,  1,          deltaTime,      0,  0,          0]
        //[0,   0,            0,                0,  0,          1,              0,  0,          0]
        //[0,   0,            0,                0,  0,          0,              1,  deltaTime,  deltaTime^2/2]
        //[0,   0,            0,                0,  0,          0,              0,  1,          deltaTime]
        //[0,   0,            0,                0,  0,          0,              0,  0,          1]

        // set elements in column order
        result.setElementAtIndex(0, 1.0, true);
        result.setElementAtIndex(1, 0.0, true);
        result.setElementAtIndex(2, 0.0, true);
        result.setElementAtIndex(3, 0.0, true);
        result.setElementAtIndex(4, 0.0, true);
        result.setElementAtIndex(5, 0.0, true);
        result.setElementAtIndex(6, 0.0, true);
        result.setElementAtIndex(7, 0.0, true);
        result.setElementAtIndex(8, 0.0, true);

        result.setElementAtIndex(9, deltaTime, true);
        result.setElementAtIndex(10, 1.0, true);
        result.setElementAtIndex(11, 0.0, true);
        result.setElementAtIndex(12, 0.0, true);
        result.setElementAtIndex(13, 0.0, true);
        result.setElementAtIndex(14, 0.0, true);
        result.setElementAtIndex(15, 0.0, true);
        result.setElementAtIndex(16, 0.0, true);
        result.setElementAtIndex(17, 0.0, true);

        result.setElementAtIndex(18, 0.5 * deltaTime * deltaTime, true);
        result.setElementAtIndex(19, deltaTime, true);
        result.setElementAtIndex(20, 1.0, true);
        result.setElementAtIndex(21, 0.0, true);
        result.setElementAtIndex(22, 0.0, true);
        result.setElementAtIndex(23, 0.0, true);
        result.setElementAtIndex(24, 0.0, true);
        result.setElementAtIndex(25, 0.0, true);
        result.setElementAtIndex(26, 0.0, true);

        result.setElementAtIndex(27, 0.0, true);
        result.setElementAtIndex(28, 0.0, true);
        result.setElementAtIndex(29, 0.0, true);
        result.setElementAtIndex(30, 1.0, true);
        result.setElementAtIndex(31, 0.0, true);
        result.setElementAtIndex(32, 0.0, true);
        result.setElementAtIndex(33, 0.0, true);
        result.setElementAtIndex(34, 0.0, true);
        result.setElementAtIndex(35, 0.0, true);

        result.setElementAtIndex(36, 0.0, true);
        result.setElementAtIndex(37, 0.0, true);
        result.setElementAtIndex(38, 0.0, true);
        result.setElementAtIndex(39, deltaTime, true);
        result.setElementAtIndex(40, 1.0, true);
        result.setElementAtIndex(41, 0.0, true);
        result.setElementAtIndex(42, 0.0, true);
        result.setElementAtIndex(43, 0.0, true);
        result.setElementAtIndex(44, 0.0, true);

        result.setElementAtIndex(45, 0.0, true);
        result.setElementAtIndex(46, 0.0, true);
        result.setElementAtIndex(47, 0.0, true);
        result.setElementAtIndex(48, 0.5 * deltaTime * deltaTime, true);
        result.setElementAtIndex(49, deltaTime, true);
        result.setElementAtIndex(50, 1.0, true);
        result.setElementAtIndex(51, 0.0, true);
        result.setElementAtIndex(52, 0.0, true);
        result.setElementAtIndex(53, 0.0, true);

        result.setElementAtIndex(54, 0.0, true);
        result.setElementAtIndex(55, 0.0, true);
        result.setElementAtIndex(56, 0.0, true);
        result.setElementAtIndex(57, 0.0, true);
        result.setElementAtIndex(58, 0.0, true);
        result.setElementAtIndex(59, 0.0, true);
        result.setElementAtIndex(60, 1.0, true);
        result.setElementAtIndex(61, 0.0, true);
        result.setElementAtIndex(62, 0.0, true);

        result.setElementAtIndex(63, 0.0, true);
        result.setElementAtIndex(64, 0.0, true);
        result.setElementAtIndex(65, 0.0, true);
        result.setElementAtIndex(66, 0.0, true);
        result.setElementAtIndex(67, 0.0, true);
        result.setElementAtIndex(68, 0.0, true);
        result.setElementAtIndex(69, deltaTime, true);
        result.setElementAtIndex(70, 1.0, true);
        result.setElementAtIndex(71, 0.0, true);

        result.setElementAtIndex(72, 0.0, true);
        result.setElementAtIndex(73, 0.0, true);
        result.setElementAtIndex(74, 0.0, true);
        result.setElementAtIndex(75, 0.0, true);
        result.setElementAtIndex(76, 0.0, true);
        result.setElementAtIndex(77, 0.0, true);
        result.setElementAtIndex(78, 0.5 * deltaTime * deltaTime, true);
        result.setElementAtIndex(79, deltaTime, true);
        result.setElementAtIndex(80, 1.0, true);
    }

    private static Matrix noiseCovarianceMatrix(final Data data) throws SignalProcessingException {

        final var estimator = new MeasurementNoiseCovarianceEstimator(3);

        final var sample = new double[3];

        for (var i = 0; i < data.numSamples; i++) {
            sample[0] = data.accelerationX[i];
            sample[1] = data.accelerationY[i];
            sample[2] = data.accelerationZ[i];

            estimator.update(sample);
        }

        return estimator.getMeasurementNoiseCov();
    }
}
