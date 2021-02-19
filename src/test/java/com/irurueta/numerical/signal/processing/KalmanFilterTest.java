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
import org.junit.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;

import static org.junit.Assert.*;

public class KalmanFilterTest {

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
    public void testConstructor() throws SignalProcessingException {
        KalmanFilter filter = new KalmanFilter(6, 9, -1);

        // check correctness
        assertEquals(filter.getDynamicParameters(), 6);
        assertEquals(filter.getMeasureParameters(), 9);
        // and because control parameters were set to negative value...
        assertEquals(filter.getControlParameters(), 6);

        assertEquals(filter.getStatePre().getRows(), 6);
        assertEquals(filter.getStatePre().getColumns(), 1);

        assertEquals(filter.getStatePost().getRows(), 6);
        assertEquals(filter.getStatePost().getColumns(), 1);

        assertEquals(filter.getTransitionMatrix().getRows(), 6);
        assertEquals(filter.getTransitionMatrix().getColumns(), 6);

        assertEquals(filter.getProcessNoiseCov().getRows(), 6);
        assertEquals(filter.getProcessNoiseCov().getColumns(), 6);

        assertEquals(filter.getMeasurementMatrix().getRows(), 9);
        assertEquals(filter.getMeasurementMatrix().getColumns(), 6);

        assertEquals(filter.getMeasurementNoiseCov().getRows(), 9);
        assertEquals(filter.getMeasurementNoiseCov().getColumns(), 9);

        assertEquals(filter.getErrorCovPre().getRows(), 6);
        assertEquals(filter.getErrorCovPre().getColumns(), 6);

        assertEquals(filter.getErrorCovPost().getRows(), 6);
        assertEquals(filter.getErrorCovPost().getColumns(), 6);

        assertEquals(filter.getGain().getRows(), 6);
        assertEquals(filter.getGain().getColumns(), 9);

        assertEquals(filter.getControlMatrix().getRows(), 6);
        assertEquals(filter.getControlMatrix().getColumns(), 6);


        // test with control parameters
        filter = new KalmanFilter(6, 9, 1);

        // check correctness
        assertEquals(filter.getDynamicParameters(), 6);
        assertEquals(filter.getMeasureParameters(), 9);
        assertEquals(filter.getControlParameters(), 1);

        assertEquals(filter.getStatePre().getRows(), 6);
        assertEquals(filter.getStatePre().getColumns(), 1);

        assertEquals(filter.getStatePost().getRows(), 6);
        assertEquals(filter.getStatePost().getColumns(), 1);

        assertEquals(filter.getTransitionMatrix().getRows(), 6);
        assertEquals(filter.getTransitionMatrix().getColumns(), 6);

        assertEquals(filter.getProcessNoiseCov().getRows(), 6);
        assertEquals(filter.getProcessNoiseCov().getColumns(), 6);

        assertEquals(filter.getMeasurementMatrix().getRows(), 9);
        assertEquals(filter.getMeasurementMatrix().getColumns(), 6);

        assertEquals(filter.getMeasurementNoiseCov().getRows(), 9);
        assertEquals(filter.getMeasurementNoiseCov().getColumns(), 9);

        assertEquals(filter.getErrorCovPre().getRows(), 6);
        assertEquals(filter.getErrorCovPre().getColumns(), 6);

        assertEquals(filter.getErrorCovPost().getRows(), 6);
        assertEquals(filter.getErrorCovPost().getColumns(), 6);

        assertEquals(filter.getGain().getRows(), 6);
        assertEquals(filter.getGain().getColumns(), 9);

        assertEquals(filter.getControlMatrix().getRows(), 6);
        assertEquals(filter.getControlMatrix().getColumns(), 1);


        // test without control parameters
        filter = new KalmanFilter(6, 9);

        // check correctness
        assertEquals(filter.getDynamicParameters(), 6);
        assertEquals(filter.getMeasureParameters(), 9);
        assertEquals(filter.getControlParameters(), 0);

        assertEquals(filter.getStatePre().getRows(), 6);
        assertEquals(filter.getStatePre().getColumns(), 1);

        assertEquals(filter.getStatePost().getRows(), 6);
        assertEquals(filter.getStatePost().getColumns(), 1);

        assertEquals(filter.getTransitionMatrix().getRows(), 6);
        assertEquals(filter.getTransitionMatrix().getColumns(), 6);

        assertEquals(filter.getProcessNoiseCov().getRows(), 6);
        assertEquals(filter.getProcessNoiseCov().getColumns(), 6);

        assertEquals(filter.getMeasurementMatrix().getRows(), 9);
        assertEquals(filter.getMeasurementMatrix().getColumns(), 6);

        assertEquals(filter.getMeasurementNoiseCov().getRows(), 9);
        assertEquals(filter.getMeasurementNoiseCov().getColumns(), 9);

        assertEquals(filter.getErrorCovPre().getRows(), 6);
        assertEquals(filter.getErrorCovPre().getColumns(), 6);

        assertEquals(filter.getErrorCovPost().getRows(), 6);
        assertEquals(filter.getErrorCovPost().getColumns(), 6);

        assertEquals(filter.getGain().getRows(), 6);
        assertEquals(filter.getGain().getColumns(), 9);

        assertNull(filter.getControlMatrix());


        // Force IllegalArgumentException
        filter = null;
        try {
            filter = new KalmanFilter(-1, 9, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }

        try {
            filter = new KalmanFilter(6, -1, 1);
            fail("IllegalArgumentException expected but not thrown");
        } catch (IllegalArgumentException ignore) {
        }

        assertNull(filter);
    }

    @Test
    public void testGetSetMeasureParameters() throws SignalProcessingException {
        final KalmanFilter filter = new KalmanFilter(6, 9, -1);

        assertEquals(filter.getMeasureParameters(), 9);

        // set new value
        filter.setMeasureParameters(10);

        // check correctness
        assertEquals(filter.getMeasureParameters(), 10);

        // Force IllegalArgumentException
        try {
            filter.setMeasureParameters(0);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testPredictAndCorrectAcceleration()
            throws SignalProcessingException, WrongSizeException {
        int numValid = 0;
        for (int t = 0; t < TIMES; t++) {

            // using a Kalman filter we take noisy measures of acceleration (for
            // simplicity we only consider one dimension).
            // Using acceleration samples we want to obtain the state of speed and
            // position
            final KalmanFilter kalman = new KalmanFilter(1, 1);

            final Random rand = new Random();
            // constant acceleration
            final double acceleration = rand.nextDouble();
            double sampleNoise;
            // double predictedAcceleration = 0.0;
            double correctedAcceleration = 0.0;


            // setup Kalman filter;
            //[ acceleration]
            Matrix correctedState;

            // measurement (aceleration)
            final Matrix measurement = new Matrix(1, 1);
            measurement.setElementAtIndex(0, acceleration);

            // transitions for acceleration
            // [1]
            final Matrix transitionMatrix = new Matrix(1, 1);
            transitionMatrix.setElementAtIndex(0, 1.0);
            kalman.setTransitionMatrix(transitionMatrix);

            if (WRITE_TO_CONSOLE) {
                System.out.println("no;truePosition;trueSpeed;trueAcceleration;" +
                        "predictedPosition;predictedSpeed;predictedAcceleration;" +
                        "correctedPosition;correctedSpeed;correctedAcceleration;");
            }

            // assume each sample happens every second
            for (int i = 0; i < N_SAMPLES; i++) {
                sampleNoise = 1e-3 * rand.nextGaussian();

                // take measure with noise
                measurement.setElementAtIndex(0,
                        measurement.getElementAtIndex(0) + sampleNoise);

                // correct
                correctedState = kalman.correct(measurement);

                correctedAcceleration = correctedState.getElementAtIndex(0);

                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";;;" + acceleration + ";" +
                            //";;" + predictedAcceleration + ";" +
                            ";;" + correctedAcceleration);
                }
            }

            if (Math.abs(acceleration - correctedAcceleration) >
                    KalmanFilter.DEFAULT_MEASUREMENT_NOISE_VARIANCE) {
                continue;
            }
            assertEquals(acceleration, correctedAcceleration,
                    KalmanFilter.DEFAULT_MEASUREMENT_NOISE_VARIANCE);

            numValid++;
            break;
        }

        assertTrue(numValid > 0);
    }

    @Test
    public void testPredictAndCorrectPositionSpeedAndAcceleration()
            throws SignalProcessingException, WrongSizeException {

        // using a Kalman filter we take noisy measures of acceleration (for
        // simplicity we only consider one dimension).
        // Using acceleration samples we want to obtain the state of speed and
        // position
        final KalmanFilter kalman = new KalmanFilter(3, 1);

        final Random rand = new Random();
        // constant acceleration
        final double acceleration = rand.nextDouble();
        double trueSpeed = 0.0;
        double truePosition = 0.0;
        double rawSpeed = 0.0;
        double rawPosition = 0.0;
        double rawAcceleration;
        double sampleNoise;
        double correctedAcceleration = 0.0;

        // setup Kalman filter;
        Matrix predictedState; //[position, speed, acceleration]
        Matrix correctedState; //[position, speed, acceleration]

        final Matrix measurement = new Matrix(1, 1); //measurement (acceleration)
        measurement.setElementAtIndex(0, acceleration);

        // assuming delta time between samples of 1 second
        final double deltaTime = 1.0;
        // transitions for position, speed and acceleration
        // [1,   deltaTime,    deltaTime^2/2]
        // [0,   1,            deltaTime]
        // [0,   0,            1]
        final Matrix transitionMatrix = new Matrix(3, 3);
        // column order
        transitionMatrix.setSubmatrix(0, 0, 2, 2, new double[]
                {
                        1.0, 0.0, 0.0,
                        deltaTime, 1.0, 0.0,
                        0.5 * deltaTime * deltaTime, deltaTime, 1.0
                }, true);
        kalman.setTransitionMatrix(transitionMatrix);

        final double processNoiseStd = Math.sqrt(1e-6);
        final double processNoiseVariance = processNoiseStd * processNoiseStd;

        // process noise covariance =
        // processNoiseVariance * [deltaTime^4/4,   deltaTime^3/2,  deltaTime^2/2]
        //                        [deltaTime^3/2,   deltaTime^2,    deltaTime]
        //                        [deltaTime^2/2,   deltaTime,      1]
        final Matrix processNoiseCov = new Matrix(3, 3);
        // column order
        processNoiseCov.setSubmatrix(0, 0, 2, 2, new double[]
                {
                        processNoiseVariance * Math.pow(deltaTime, 4.0) * 0.25, processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5, processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5, processNoiseVariance * Math.pow(deltaTime, 2.0), processNoiseVariance * deltaTime,
                        processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5, processNoiseVariance * deltaTime, processNoiseVariance
                }, true);
        kalman.setProcessNoiseCov(processNoiseCov);

        final double measureNoiseStd = Math.sqrt(1e-1);
        final double measureNoiseVariance = measureNoiseStd * measureNoiseStd;

        // measurement noise covariance = [measureNoiseVariance]
        final Matrix measurementNoiseCov = new Matrix(1, 1);
        measurementNoiseCov.setElementAtIndex(0, measureNoiseVariance);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);


        final Matrix measurementMatrix = new Matrix(1, 3);
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
        for (int i = 0; i < N_SAMPLES; i++) {
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
    public void testPredictAndCorrectPositionSpeedAndAcceleration3D()
            throws SignalProcessingException, WrongSizeException {

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final KalmanFilter kalman = new KalmanFilter(9, 3);

        final Random rand = new Random();
        // constant acceleration
        final double accelerationX = rand.nextDouble();
        final double accelerationY = rand.nextDouble();
        final double accelerationZ = rand.nextDouble();
        double trueSpeedX = 0.0;
        double trueSpeedY = 0.0;
        double trueSpeedZ = 0.0;
        double truePositionX = 0.0;
        double truePositionY = 0.0;
        double truePositionZ = 0.0;
        double rawSpeedX = 0.0;
        double rawSpeedY = 0.0;
        double rawSpeedZ = 0.0;
        double rawPositionX = 0.0;
        double rawPositionY = 0.0;
        double rawPositionZ = 0.0;
        double rawAccelerationX;
        double rawAccelerationY;
        double rawAccelerationZ;
        double sampleNoiseX;
        double sampleNoiseY;
        double sampleNoiseZ;
        double correctedAccelerationX = 0.0;
        double correctedAccelerationY = 0.0;
        double correctedAccelerationZ = 0.0;

        // setup Kalman filter;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix predictedState;
        // [positionX, speedX, accelerationX, positionY, speedY, accelerationY, positionZ, speedZ, accelerationZ]
        Matrix correctedState;

        // measurement (accelerationX, accelerationY, accelerationZ)
        final Matrix measurement = new Matrix(3, 1);
        measurement.setElementAtIndex(0, accelerationX);
        measurement.setElementAtIndex(1, accelerationY);
        measurement.setElementAtIndex(2, accelerationZ);

        // assuming delta time betweeen samples of 1 second
        final double deltaTime = 1.0;
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
        final Matrix transitionMatrix = new Matrix(9, 9);
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

        final double processNoiseStd = Math.sqrt(1e-6);
        final double processNoiseVariance = processNoiseStd * processNoiseStd;

        // process noise covariance is a symmetric matrix obtained as a result of:
        // processNoiseVariance * [deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1]'*[deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1, deltaTime^2/2, deltaTime, 1]
        // which results in a symmetric matrix containing 9 times whe following block:
        // processNoiseVariance * [deltaTime^4/4,   deltaTime^3/2,  deltaTime^2/2]
        //                        [deltaTime^3/2,   deltaTime^2,    deltaTime]
        //                        [deltaTime^2/2,   deltaTime,      1]
        final Matrix block = new Matrix(3, 3);
        // column order
        block.setSubmatrix(0, 0, 2, 2, new double[]
                {
                        processNoiseVariance * Math.pow(deltaTime, 4.0) * 0.25, processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5, processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5,
                        processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5, processNoiseVariance * Math.pow(deltaTime, 2.0), processNoiseVariance * deltaTime,
                        processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5, processNoiseVariance * deltaTime, processNoiseVariance
                }, true);
        Matrix processNoiseCov = new Matrix(9, 9);
        for (int i = 0; i < 9; i += 3) {
            for (int j = 0; j < 9; j += 3) {
                processNoiseCov.setSubmatrix(i, j, i + 2, j + 2, block);
            }
        }
        kalman.setProcessNoiseCov(processNoiseCov);

        final double measureNoiseStd = Math.sqrt(1e-1);
        final double measureNoiseVariance = measureNoiseStd * measureNoiseStd;

        // measurement noise covariance = measureNoiseVariance * I
        final Matrix measurementNoiseCov = Matrix.identity(3, 3);
        measurementNoiseCov.multiplyByScalar(measureNoiseVariance);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final Matrix measurementMatrix = new Matrix(3, 9);
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
        for (int i = 0; i < N_SAMPLES; i++) {
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
                    System.out.println(i + ";" + truePositionX + ";" + trueSpeedX + ";" + accelerationX + ";" +
                            correctedAccelerationX + ";" +
                            rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" +
                            truePositionY + ";" + trueSpeedY + ";" + accelerationY + ";" +
                            correctedAccelerationY + ";" +
                            rawPositionY + ";" + rawSpeedY + ";" + rawAccelerationY + ";" +
                            truePositionZ + ";" + trueSpeedZ + ";" + accelerationZ + ";" +
                            correctedAccelerationZ + ";" +
                            rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
                }

            } else {
                if (WRITE_TO_CONSOLE) {
                    System.out.println(i + ";" + truePositionX + ";" + trueSpeedX + ";" + accelerationX + ";" +
                            ";;;" +
                            rawPositionX + ";" + rawSpeedX + ";" + rawAccelerationX + ";" +
                            truePositionY + ";" + trueSpeedY + ";" + accelerationY + ";" +
                            ";;;" +
                            rawPositionY + ";" + rawSpeedY + ";" + rawAccelerationY + ";" +
                            truePositionZ + ";" + trueSpeedZ + ";" + accelerationZ + ";" +
                            ";;;" +
                            rawPositionZ + ";" + rawSpeedZ + ";" + rawAccelerationZ);
                }
            }
        }

        // check that the filter state converged to the true acceleration
        assertEquals(accelerationX, correctedAccelerationX, measureNoiseVariance);
        assertEquals(accelerationY, correctedAccelerationY, measureNoiseVariance);
        assertEquals(accelerationZ, correctedAccelerationZ, measureNoiseVariance);
    }

    @Test
    public void testPredictAndCorrectRealDataNoMotion() throws IOException,
            SignalProcessingException, WrongSizeException {

        File f = new File(ACCELERATION_NO_MOTION);
        final Data data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final KalmanFilter kalman = new KalmanFilter(9, 3);

        double rawSpeedX = 0.0;
        double rawSpeedY = 0.0;
        double rawSpeedZ = 0.0;
        double rawPositionX = 0.0;
        double rawPositionY = 0.0;
        double rawPositionZ = 0.0;
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
        final Matrix measurement = new Matrix(3, 1);

        final Matrix transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final double processNoiseStd = Math.sqrt(1e-6);
        final double processNoiseVariance = processNoiseStd * processNoiseStd;

        final Matrix block = new Matrix(3, 3);
        final Matrix processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block,
                processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final Matrix measurementNoiseCov = noiseCovarianceMatrix(data);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final Matrix measurementMatrix = new Matrix(3, 9);
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
        final PrintStream stream = new PrintStream(f);

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

        double deltaTime = 20e-3;
        long prevTimestampNanos;
        long timestampNanos = 0;
        double time = 0.0;
        for (int i = 0; i < data.numSamples; i++) {
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

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block,
                    processNoiseCov);
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
    public void testPredictAndCorrectRealDataMotion() throws IOException,
            SignalProcessingException, WrongSizeException {

        File f = new File(ACCELERATION_NO_MOTION);
        final Data dataNoMotion = AccelerationFileLoader.load(f);

        f = new File(ACCELERATION_MOTION);
        final Data data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final KalmanFilter kalman = new KalmanFilter(9, 3);

        double rawSpeedX = 0.0;
        double rawSpeedY = 0.0;
        double rawSpeedZ = 0.0;
        double rawPositionX = 0.0;
        double rawPositionY = 0.0;
        double rawPositionZ = 0.0;
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
        final Matrix measurement = new Matrix(3, 1);

        final Matrix transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final double processNoiseStd = Math.sqrt(1e-6);
        final double processNoiseVariance = processNoiseStd * processNoiseStd;

        final Matrix block = new Matrix(3, 3);
        final Matrix processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block,
                processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final Matrix measurementNoiseCov = noiseCovarianceMatrix(dataNoMotion);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final Matrix measurementMatrix = new Matrix(3, 9);
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
        final PrintStream stream = new PrintStream(f);

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

        double deltaTime = 20e-3;
        long prevTimestampNanos;
        long timestampNanos = 0;
        double time = 0.0;
        for (int i = 0; i < data.numSamples; i++) {
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

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block,
                    processNoiseCov);
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
    public void testPredictAndCorrectRealDataNoMotionFast() throws IOException,
            SignalProcessingException, WrongSizeException {

        File f = new File(ACCELERATION_NO_MOTION_FAST);
        final Data data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final KalmanFilter kalman = new KalmanFilter(9, 3);

        double rawSpeedX = 0.0;
        double rawSpeedY = 0.0;
        double rawSpeedZ = 0.0;
        double rawPositionX = 0.0;
        double rawPositionY = 0.0;
        double rawPositionZ = 0.0;
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
        final Matrix measurement = new Matrix(3, 1);

        final Matrix transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final double processNoiseStd = Math.sqrt(1e-6);
        final double processNoiseVariance = processNoiseStd * processNoiseStd;

        final Matrix block = new Matrix(3, 3);
        final Matrix processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block,
                processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final Matrix measurementNoiseCov = noiseCovarianceMatrix(data);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final Matrix measurementMatrix = new Matrix(3, 9);
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
        final PrintStream stream = new PrintStream(f);


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

        double deltaTime = 4e-3;
        long prevTimestampNanos;
        long timestampNanos = 0;
        double time = 0.0;
        for (int i = 0; i < data.numSamples; i++) {
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

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block,
                    processNoiseCov);
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
    public void testPredictAndCorrectRealDataMotionFast() throws IOException,
            SignalProcessingException, WrongSizeException {

        File f = new File(ACCELERATION_NO_MOTION_FAST);
        final Data dataNoMotion = AccelerationFileLoader.load(f);

        f = new File(ACCELERATION_MOTION_FAST);
        final Data data = AccelerationFileLoader.load(f);

        // using a Kalman filter we take noisy measures of acceleration in 3D.
        // Using acceleration samples we want to obtain the state of speed and
        // position in 3D
        final KalmanFilter kalman = new KalmanFilter(9, 3);

        double rawSpeedX = 0.0;
        double rawSpeedY = 0.0;
        double rawSpeedZ = 0.0;
        double rawPositionX = 0.0;
        double rawPositionY = 0.0;
        double rawPositionZ = 0.0;
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
        final Matrix measurement = new Matrix(3, 1);

        final Matrix transitionMatrix = new Matrix(9, 9);
        updateTransitionMatrix(1.0, transitionMatrix);
        kalman.setTransitionMatrix(transitionMatrix);

        final double processNoiseStd = Math.sqrt(1e-6);
        final double processNoiseVariance = processNoiseStd * processNoiseStd;

        final Matrix block = new Matrix(3, 3);
        final Matrix processNoiseCov = new Matrix(9, 9);
        updateProcessNoiseCov(processNoiseVariance, 1.0, block,
                processNoiseCov);
        kalman.setProcessNoiseCov(processNoiseCov);

        // noise covariance matrix is obtained by sampling data when system state
        // is held constant (no motion). Ideally noise should be statistically
        // independent (i.e. a diagonal noise matrix)
        final Matrix measurementNoiseCov = noiseCovarianceMatrix(dataNoMotion);
        kalman.setMeasurementNoiseCov(measurementNoiseCov);

        // measurement matrix relates system state to measured data. This matrix
        // for the acceleration case is constant and has the following form:
        // H = [0, 0, 1, 0, 0, 0, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 1, 0, 0, 0]
        //     [0, 0, 0, 0, 0, 0, 0, 0, 1]
        final Matrix measurementMatrix = new Matrix(3, 9);
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
        final PrintStream stream = new PrintStream(f);

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

        double deltaTime = 4e-3;
        long prevTimestampNanos;
        long timestampNanos = 0;
        double time = 0.0;
        for (int i = 0; i < data.numSamples; i++) {
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

            updateProcessNoiseCov(processNoiseVariance, deltaTime, block,
                    processNoiseCov);
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
    public void testGetSetStatePre() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        // check default value
        final Matrix statePre = filter.getStatePre();

        assertEquals(statePre.getRows(), 6);
        assertEquals(statePre.getColumns(), 1);

        // set new value
        final Matrix statePre2 = new Matrix(6, 1);
        filter.setStatePre(statePre2);

        // check correctness
        assertSame(filter.getStatePre(), statePre2);

        // Force IllegalArgumentException
        Matrix wrongStatePre = new Matrix(5, 1);
        try {
            filter.setStatePre(wrongStatePre);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongStatePre = new Matrix(6, 2);
        try {
            filter.setStatePre(wrongStatePre);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetStatePost() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        // check default value
        final Matrix statePost = filter.getStatePost();

        assertEquals(statePost.getRows(), 6);
        assertEquals(statePost.getColumns(), 1);

        // set new value
        final Matrix statePost2 = new Matrix(6, 1);
        filter.setStatePost(statePost2);

        // check correctness
        assertSame(filter.getStatePost(), statePost2);

        // Force IllegalArgumentException
        Matrix wrongStatePost = new Matrix(5, 1);
        try {
            filter.setStatePost(wrongStatePost);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongStatePost = new Matrix(6, 2);
        try {
            filter.setStatePost(wrongStatePost);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetTransitionMatrix() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        // check default value
        final Matrix transitionMatrix = filter.getTransitionMatrix();

        assertEquals(transitionMatrix.getRows(), 6);
        assertEquals(transitionMatrix.getColumns(), 6);

        // set new value
        final Matrix transitionMatrix2 = new Matrix(6, 6);
        filter.setTransitionMatrix(transitionMatrix2);

        // check correctness
        assertSame(filter.getTransitionMatrix(), transitionMatrix2);

        // Force IllegalArgumentException
        Matrix wrongTransitionMatrix = new Matrix(5, 6);
        try {
            filter.setTransitionMatrix(wrongTransitionMatrix);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongTransitionMatrix = new Matrix(6, 5);
        try {
            filter.setTransitionMatrix(wrongTransitionMatrix);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetControlMatrix() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9, 1);

        // check default value
        final Matrix controlMatrix = filter.getControlMatrix();

        assertEquals(controlMatrix.getRows(), 6);
        assertEquals(controlMatrix.getColumns(), 1);

        // set new value
        final Matrix controlMatrix2 = new Matrix(6, 1);
        filter.setControlMatrix(controlMatrix2);

        // check correctness
        assertSame(filter.getControlMatrix(), controlMatrix2);

        // Force IllegalArgumentException
        Matrix wrongControlMatrix = new Matrix(5, 1);
        try {
            filter.setControlMatrix(wrongControlMatrix);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongControlMatrix = new Matrix(6, 2);
        try {
            filter.setControlMatrix(wrongControlMatrix);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetMeasurementMatrix() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        // check default value
        final Matrix measurementMatrix = filter.getMeasurementMatrix();

        assertEquals(measurementMatrix.getRows(), 9);
        assertEquals(measurementMatrix.getColumns(), 6);

        // set new value
        final Matrix measurementMatrix2 = new Matrix(9, 6);
        filter.setMeasurementMatrix(measurementMatrix2);

        // check correctness
        assertSame(filter.getMeasurementMatrix(), measurementMatrix2);

        // Force IllegalArgumentException
        Matrix wrongMeasurementMatrix = new Matrix(8, 6);
        try {
            filter.setMeasurementMatrix(wrongMeasurementMatrix);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongMeasurementMatrix = new Matrix(9, 5);
        try {
            filter.setMeasurementMatrix(wrongMeasurementMatrix);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetProcessNoiseCov() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        // check default value
        final Matrix processNoiseCov = filter.getProcessNoiseCov();

        assertEquals(processNoiseCov.getRows(), 6);
        assertEquals(processNoiseCov.getColumns(), 6);

        // set new value
        final Matrix processNoiseCov2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        filter.setProcessNoiseCov(processNoiseCov2);

        // check correctness
        assertSame(filter.getProcessNoiseCov(), processNoiseCov2);

        // Force IllegalArgumentException
        Matrix wrongProcessNoiseCov = new Matrix(5, 6);
        try {
            filter.setProcessNoiseCov(wrongProcessNoiseCov);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongProcessNoiseCov = new Matrix(6, 5);
        try {
            filter.setProcessNoiseCov(wrongProcessNoiseCov);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongProcessNoiseCov = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        wrongProcessNoiseCov.setElementAt(0, 1, 1.0);
        try {
            filter.setProcessNoiseCov(wrongProcessNoiseCov);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetMeasurementNoiseCov()
            throws SignalProcessingException, WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        // check default value
        final Matrix measurementNoiseCov = filter.getMeasurementNoiseCov();

        assertEquals(measurementNoiseCov.getRows(), 9);
        assertEquals(measurementNoiseCov.getColumns(), 9);

        // set new value
        final Matrix measurementNoiseCov2 = Matrix.diagonal(
                new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9});
        filter.setMeasurementNoiseCov(measurementNoiseCov2);

        // check correctness
        assertSame(filter.getMeasurementNoiseCov(), measurementNoiseCov2);

        // Force IllegalArgumentException
        Matrix wrongMeasurmentNoiseCov = new Matrix(8, 9);
        try {
            filter.setMeasurementNoiseCov(wrongMeasurmentNoiseCov);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongMeasurmentNoiseCov = new Matrix(9, 8);
        try {
            filter.setMeasurementNoiseCov(wrongMeasurmentNoiseCov);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongMeasurmentNoiseCov = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9});
        wrongMeasurmentNoiseCov.setElementAt(0, 1, 1.0);
        try {
            filter.setMeasurementNoiseCov(wrongMeasurmentNoiseCov);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetErrorCovPre() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        final Matrix errorCovPre = filter.getErrorCovPre();

        assertEquals(errorCovPre.getRows(), 6);
        assertEquals(errorCovPre.getColumns(), 6);

        // set new value
        final Matrix errorCovPre2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        filter.setErrorCovPre(errorCovPre2);

        // check correctness
        assertSame(filter.getErrorCovPre(), errorCovPre2);

        // Force IllegalArgumentException
        Matrix wrongErrorCovPre = new Matrix(5, 6);
        try {
            filter.setErrorCovPre(wrongErrorCovPre);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongErrorCovPre = new Matrix(6, 5);
        try {
            filter.setErrorCovPre(wrongErrorCovPre);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongErrorCovPre = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        wrongErrorCovPre.setElementAt(0, 1, 1.0);
        try {
            filter.setErrorCovPre(wrongErrorCovPre);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetGain() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        final Matrix gain = filter.getGain();

        assertEquals(gain.getRows(), 6);
        assertEquals(gain.getColumns(), 9);

        // set new value
        final Matrix gain2 = new Matrix(6, 9);
        filter.setGain(gain2);

        // check correctness
        assertSame(filter.getGain(), gain2);

        // Force IllegalArgumentException
        Matrix wrongGain = new Matrix(5, 9);
        try {
            filter.setGain(wrongGain);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongGain = new Matrix(6, 8);
        try {
            filter.setGain(wrongGain);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    @Test
    public void testGetSetErrorCovPost() throws SignalProcessingException,
            WrongSizeException {
        final KalmanFilter filter = new KalmanFilter(6, 9);

        final Matrix errorCovPost = filter.getErrorCovPost();

        assertEquals(errorCovPost.getRows(), 6);
        assertEquals(errorCovPost.getColumns(), 6);

        // set new value
        final Matrix errorCovPost2 = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        filter.setErrorCovPost(errorCovPost2);

        // check correctness
        assertSame(filter.getErrorCovPost(), errorCovPost2);

        // Force IllegalArgumentException
        Matrix wrongErrorCovPost = new Matrix(5, 6);
        try {
            filter.setErrorCovPost(wrongErrorCovPost);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongErrorCovPost = new Matrix(6, 5);
        try {
            filter.setErrorCovPost(wrongErrorCovPost);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }

        wrongErrorCovPost = Matrix.diagonal(new double[]{1, 2, 3, 4, 5, 6});
        wrongErrorCovPost.setElementAt(0, 1, 1.0);
        try {
            filter.setErrorCovPost(wrongErrorCovPost);
            fail("IllegalArgumentException expected but not thrown");
        } catch (final IllegalArgumentException ignore) {
        }
    }

    private void updateProcessNoiseCov(final double processNoiseVariance,
                                       final double deltaTime, final Matrix block, final Matrix result) {

        // block matrix has the following form:
        // processNoiseVariance * [deltaTime^4/4,   deltaTime^3/2,  deltaTime^2/2]
        //                        [deltaTime^3/2,   deltaTime^2,    deltaTime]
        //                        [deltaTime^2/2,   deltaTime,      1]

        // set elements in column order
        block.setElementAtIndex(0,
                processNoiseVariance * Math.pow(deltaTime, 4.0) * 0.25, true);
        block.setElementAtIndex(1,
                processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5, true);
        block.setElementAtIndex(2,
                processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5, true);

        block.setElementAtIndex(3,
                processNoiseVariance * Math.pow(deltaTime, 3.0) * 0.5, true);
        block.setElementAtIndex(4,
                processNoiseVariance * Math.pow(deltaTime, 2.0), true);
        block.setElementAtIndex(5,
                processNoiseVariance * deltaTime, true);

        block.setElementAtIndex(6,
                processNoiseVariance * Math.pow(deltaTime, 2.0) * 0.5, true);
        block.setElementAtIndex(7,
                processNoiseVariance * deltaTime, true);
        block.setElementAtIndex(8, processNoiseVariance);

        // repeat block matrix throughout result matrix (9 times)
        for (int i = 0; i < 9; i += 3) {
            for (int j = 0; j < 9; j += 3) {
                result.setSubmatrix(i, j, i + 2, j + 2, block);
            }
        }
    }

    private void updateTransitionMatrix(final double deltaTime, final Matrix result) {
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

    private Matrix noiseCovarianceMatrix(final Data data)
            throws SignalProcessingException {

        final MeasurementNoiseCovarianceEstimator estimator =
                new MeasurementNoiseCovarianceEstimator(3);

        final double[] sample = new double[3];

        for (int i = 0; i < data.numSamples; i++) {
            sample[0] = data.accelerationX[i];
            sample[1] = data.accelerationY[i];
            sample[2] = data.accelerationZ[i];

            estimator.update(sample);
        }

        return estimator.getMeasurementNoiseCov();
    }
}
