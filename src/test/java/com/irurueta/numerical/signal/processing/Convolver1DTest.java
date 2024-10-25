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
package com.irurueta.numerical.signal.processing;

import com.irurueta.statistics.UniformRandomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class Convolver1DTest implements Convolver1D.Convolver1DListener {

    public static final double MIN_RANDOM_VALUE = -10.0;
    public static final double MAX_RANDOM_VALUE = 10.0;

    private int startConvolution;
    private int finishConvolution;
    private int convolveProgressChange;

    @Test
    void testConstructor() {
        // empty constructor
        var convolver = new Convolver1D();

        // check default values
        assertNull(convolver.getSignal());
        assertNull(convolver.getKernel());
        assertEquals(0, convolver.getKernelCenter());
        assertEquals(ConvolverEdgeMethod.ZERO_EDGE, convolver.getEdgeMethod());
        assertEquals(0.0, convolver.getConstantValue(), 0.0);
        assertNull(convolver.getListener());
        assertFalse(convolver.isReady());

        // constructor with signal and kernel
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        final double[] kernel = {1.0, 2.0, 1.0};
        convolver = new Convolver1D(signal, kernel);

        // check default values
        assertSame(signal, convolver.getSignal());
        assertSame(kernel, convolver.getKernel());
        assertEquals(0, convolver.getKernelCenter());
        assertEquals(ConvolverEdgeMethod.ZERO_EDGE, convolver.getEdgeMethod());
        assertEquals(0.0, convolver.getConstantValue(), 0.0);
        assertNull(convolver.getListener());
        assertTrue(convolver.isReady());
    }

    @Test
    void testGetSetSignal() {
        final var convolver = new Convolver1D();

        // check default value
        assertNull(convolver.getSignal());

        // set new value
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        convolver.setSignal(signal);

        // check correctness
        assertSame(convolver.getSignal(), signal);
    }

    @Test
    void testGetSetKernel() {
        final var convolver = new Convolver1D();

        // check default value
        assertNull(convolver.getKernel());

        // set new value
        final double[] kernel = {1.0, 2.0, 1.0};
        convolver.setKernel(kernel);

        // check correctness
        assertSame(convolver.getKernel(), kernel);
    }

    @Test
    void testGetSetKernelCenter() {
        final var convolver = new Convolver1D();

        // check default value
        assertEquals(0, convolver.getKernelCenter());

        // set new value
        convolver.setKernelCenter(2);

        // check correctness
        assertEquals(2, convolver.getKernelCenter());

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> convolver.setKernelCenter(-1));
    }

    @Test
    void testGetSetEdgeMethod() {
        final var convolver = new Convolver1D();

        // check default value
        assertEquals(ConvolverEdgeMethod.ZERO_EDGE, convolver.getEdgeMethod());

        // set new value
        convolver.setEdgeMethod(ConvolverEdgeMethod.REPEAT_EDGE);

        // check correctness
        assertEquals(ConvolverEdgeMethod.REPEAT_EDGE, convolver.getEdgeMethod());
    }

    @Test
    void testGetSetConstantValue() {
        final var randomizer = new UniformRandomizer();

        final var constantValue = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final var convolver = new Convolver1D();

        // default value
        assertEquals(0.0, convolver.getConstantValue(), 0.0);

        // set new value
        convolver.setConstantValue(constantValue);

        // check correctness
        assertEquals(constantValue, convolver.getConstantValue(), 0.0);
    }

    @Test
    void testGetSetListener() {
        final var convolver = new Convolver1D();

        // check default value
        assertNull(convolver.getListener());

        // new value
        convolver.setListener(this);

        // check correctness
        assertSame(this, convolver.getListener());
    }

    @Test
    void testIsReady() {
        final var convolver = new Convolver1D();

        // check default value
        assertFalse(convolver.isReady());

        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        convolver.setSignal(signal);

        // check correctness
        assertFalse(convolver.isReady());

        final double[] kernel = {1.0, 2.0, 1.0};
        convolver.setKernel(kernel);

        // check correctness
        assertTrue(convolver.isReady());

        convolver.setKernelCenter(3);

        // check correctness
        assertFalse(convolver.isReady());

        convolver.setKernelCenter(2);

        // check correctness
        assertTrue(convolver.isReady());
    }

    @Test
    void testConvolveZero() {

        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        final double[] kernel = {1.0, 2.0, 1.0};

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);

        // test: center = 0, zero edge
        final var result = new double[7];
        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE, 0.0, result,
                this);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 0,
                ConvolverEdgeMethod.ZERO_EDGE, 0.0, new double[1], this));
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 3,
                ConvolverEdgeMethod.ZERO_EDGE, 0.0, result, this));

        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE, 0.0, result);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        var result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE,
                0.0, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE, 0.0);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);

        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE, result, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE, result);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.ZERO_EDGE);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        Convolver1D.convolve(signal, kernel, 0, result, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 0, result);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, 0, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 0);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        Convolver1D.convolve(signal, kernel, result, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, result);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, this);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        final Convolver1D convolver = new Convolver1D(signal, kernel);
        convolver.setListener(this);
        convolver.convolve(result);

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = convolver.convolve();

        // check correctness
        assertArrayEquals(new double[]{1.0, 4.0, 8.0, 12.0, 16.0, 14.0, 5.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();


        // test: center = 1, zero edge
        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, 0.0, result,
                this);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, 0.0, result);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0},
                result, 0.0);


        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, 0.0,
                this);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, 0.0);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result2, 0.0);

        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, result, this);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, result);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE, this);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.ZERO_EDGE);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result2, 0.0);

        Convolver1D.convolve(signal, kernel, 1, result, this);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 1, result);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, 1, this);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 1);

        // check correctness
        assertArrayEquals(new double[]{0.0, 1.0, 4.0, 8.0, 12.0, 16.0, 14.0}, result2, 0.0);
    }

    @Test
    void testConvolveConstant() {

        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        final double[] kernel = {1.0, 2.0, 1.0};
        final double constantValue = 10.0;

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);

        // test: center = 0, constant edge
        final var result = new double[7];
        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue, result,
                this);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        // check correctness
        assertArrayEquals(new double[]{31.0, 14.0, 8.0, 12.0, 16.0, 24.0, 35.0}, result, 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 0,
                ConvolverEdgeMethod.CONSTANT_EDGE, constantValue, new double[1], this));
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 3,
                ConvolverEdgeMethod.CONSTANT_EDGE, constantValue, result, this));

        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue, result);

        // check correctness
        assertArrayEquals(new double[]{31.0, 14.0, 8.0, 12.0, 16.0, 24.0, 35.0}, result, 0.0);

        var result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.CONSTANT_EDGE,
                constantValue, this);

        // check correctness
        assertArrayEquals(new double[]{31.0, 14.0, 8.0, 12.0, 16.0, 24.0, 35.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue);

        // check correctness
        assertArrayEquals(new double[]{31.0, 14.0, 8.0, 12.0, 16.0, 24.0, 35.0}, result2, 0.0);

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);


        // test: center = 1, constant edge
        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue, result,
                this);

        // check correctness
        assertArrayEquals(new double[]{40.0, 31.0, 14.0, 8.0, 12.0, 16.0, 24.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue, result);

        // check correctness
        assertArrayEquals(new double[]{40.0, 31.0, 14.0, 8.0, 12.0, 16.0, 24.0}, result, 0.0);


        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue,
                this);

        // check correctness
        assertArrayEquals(new double[]{40.0, 31.0, 14.0, 8.0, 12.0, 16.0, 24.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.CONSTANT_EDGE, constantValue);

        // check correctness
        assertArrayEquals(new double[]{40.0, 31.0, 14.0, 8.0, 12.0, 16.0, 24.0}, result2, 0.0);
    }

    @Test
    void testConvolveRepeat() {
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        final double[] kernel = {1.0, 2.0, 1.0};

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);

        // test: center = 0, repeat edge
        final var result = new double[7];
        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.REPEAT_EDGE, 0.0, result,
                this);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        // check correctness
        assertArrayEquals(new double[]{15.0, 9.0, 8.0, 12.0, 16.0, 15.0, 9.0}, result, 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 0,
                ConvolverEdgeMethod.REPEAT_EDGE, 0.0, new double[1], this));
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 3,
                ConvolverEdgeMethod.REPEAT_EDGE, 0.0, result, this));

        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.REPEAT_EDGE, 0.0, result);

        // check correctness
        assertArrayEquals(new double[]{15.0, 9.0, 8.0, 12.0, 16.0, 15.0, 9.0}, result, 0.0);

        var result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.REPEAT_EDGE,
                0.0, this);

        // check correctness
        assertArrayEquals(new double[]{15.0, 9.0, 8.0, 12.0, 16.0, 15.0, 9.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.REPEAT_EDGE,
                0.0);

        // check correctness
        assertArrayEquals(new double[]{15.0, 9.0, 8.0, 12.0, 16.0, 15.0, 9.0}, result2, 0.0);

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);

        // test: center = 1, repeat edge
        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.REPEAT_EDGE, 0.0, result,
                this);

        // check correctness
        assertArrayEquals(new double[]{16.0, 15.0, 9.0, 8.0, 12.0, 16.0, 15.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.REPEAT_EDGE, 0.0, result);

        // check correctness
        assertArrayEquals(new double[]{16.0, 15.0, 9.0, 8.0, 12.0, 16.0, 15.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.REPEAT_EDGE, 0.0,
                this);

        // check correctness
        assertArrayEquals(new double[]{16.0, 15.0, 9.0, 8.0, 12.0, 16.0, 15.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.REPEAT_EDGE,
                0.0);

        // check correctness
        assertArrayEquals(new double[]{16.0, 15.0, 9.0, 8.0, 12.0, 16.0, 15.0}, result2, 0.0);
    }

    @Test
    void testConvolveMirror() {
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};
        final double[] kernel = {1.0, 2.0, 1.0};

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);

        // test: center = 0, mirror edge
        final var result = new double[7];
        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.MIRROR_EDGE, 0.0, result,
                this);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        // check correctness
        assertArrayEquals(new double[]{5.0, 5.0, 8.0, 12.0, 16.0, 19.0, 19.0}, result, 0.0);

        // Force IllegalArgumentException
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 0,
                ConvolverEdgeMethod.MIRROR_EDGE, 0.0, new double[1], this));
        assertThrows(IllegalArgumentException.class, () -> Convolver1D.convolve(signal, kernel, 3,
                ConvolverEdgeMethod.MIRROR_EDGE, 0.0, result, this));

        Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.MIRROR_EDGE, 0.0, result);

        // check correctness
        assertArrayEquals(new double[]{5.0, 5.0, 8.0, 12.0, 16.0, 19.0, 19.0}, result, 0.0);

        var result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.MIRROR_EDGE,
                0.0, this);

        // check correctness
        assertArrayEquals(new double[]{5.0, 5.0, 8.0, 12.0, 16.0, 19.0, 19.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        result2 = Convolver1D.convolve(signal, kernel, 0, ConvolverEdgeMethod.MIRROR_EDGE, 0.0);

        // check correctness
        assertArrayEquals(new double[]{5.0, 5.0, 8.0, 12.0, 16.0, 19.0, 19.0}, result2, 0.0);

        reset();

        assertEquals(0, startConvolution);
        assertEquals(0, finishConvolution);
        assertEquals(0, convolveProgressChange);


        // test: center = 1, mirror edge
        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.MIRROR_EDGE, 0.0, result,
                this);

        // check correctness
        assertArrayEquals(new double[]{8.0, 5.0, 5.0, 8.0, 12.0, 16.0, 19.0}, result, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.MIRROR_EDGE, 0.0, result);

        // check correctness
        assertArrayEquals(new double[]{8.0, 5.0, 5.0, 8.0, 12.0, 16.0, 19.0}, result, 0.0);

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.MIRROR_EDGE, 0.0,
                this);

        // check correctness
        assertArrayEquals(new double[]{8.0, 5.0, 5.0, 8.0, 12.0, 16.0, 19.0}, result2, 0.0);

        assertEquals(1, startConvolution);
        assertEquals(1, finishConvolution);
        assertEquals(7, convolveProgressChange);

        reset();

        result2 = Convolver1D.convolve(signal, kernel, 1, ConvolverEdgeMethod.MIRROR_EDGE,
                0.0);

        // check correctness
        assertArrayEquals(new double[]{8.0, 5.0, 5.0, 8.0, 12.0, 16.0, 19.0}, result2, 0.0);
    }

    @Test
    void testGetSignalValueZero() {
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};

        assertEquals(0.0, Convolver1D.getSignalValueZero(signal, -1), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueZero(signal, 0), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueZero(signal, 1), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueZero(signal, 2), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueZero(signal, 3), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueZero(signal, 4), 0.0);
        assertEquals(0.0, Convolver1D.getSignalValueZero(signal, 5), 0.0);
    }

    @Test
    void testGetSignalValueConstant() {
        final var randomizer = new UniformRandomizer();
        final var constantValue = randomizer.nextDouble(MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};

        assertEquals(Convolver1D.getSignalValueConstant(signal, -1, constantValue), constantValue, 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueConstant(signal, 0, constantValue), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueConstant(signal, 1, constantValue), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueConstant(signal, 2, constantValue), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueConstant(signal, 3, constantValue), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueConstant(signal, 4, constantValue), 0.0);
        assertEquals(Convolver1D.getSignalValueConstant(signal, 5, constantValue), constantValue, 0.0);
    }

    @Test
    void testGetSignalValueRepeat() {
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};

        assertEquals(4.0, Convolver1D.getSignalValueRepeat(signal, -12), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueRepeat(signal, -11), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueRepeat(signal, -10), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueRepeat(signal, -9), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueRepeat(signal, -8), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueRepeat(signal, -7), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueRepeat(signal, -6), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueRepeat(signal, -5), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueRepeat(signal, -4), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueRepeat(signal, -3), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueRepeat(signal, -2), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueRepeat(signal, -1), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueRepeat(signal, 0), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueRepeat(signal, 1), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueRepeat(signal, 2), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueRepeat(signal, 3), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueRepeat(signal, 4), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueRepeat(signal, 5), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueRepeat(signal, 6), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueRepeat(signal, 7), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueRepeat(signal, 8), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueRepeat(signal, 9), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueRepeat(signal, 10), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueRepeat(signal, 11), 0.0);
    }

    @Test
    void testGetSignalValueMirror() {
        final double[] signal = {1.0, 2.0, 3.0, 4.0, 5.0};

        assertEquals(2.0, Convolver1D.getSignalValueMirror(signal, -12), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueMirror(signal, -11), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueMirror(signal, -10), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueMirror(signal, -9), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueMirror(signal, -8), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueMirror(signal, -7), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueMirror(signal, -6), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueMirror(signal, -5), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueMirror(signal, -4), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueMirror(signal, -3), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueMirror(signal, -2), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueMirror(signal, -1), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueMirror(signal, 0), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueMirror(signal, 1), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueMirror(signal, 2), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueMirror(signal, 3), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueMirror(signal, 4), 0.0);
        assertEquals(5.0, Convolver1D.getSignalValueMirror(signal, 5), 0.0);
        assertEquals(4.0, Convolver1D.getSignalValueMirror(signal, 6), 0.0);
        assertEquals(3.0, Convolver1D.getSignalValueMirror(signal, 7), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueMirror(signal, 8), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueMirror(signal, 9), 0.0);
        assertEquals(1.0, Convolver1D.getSignalValueMirror(signal, 10), 0.0);
        assertEquals(2.0, Convolver1D.getSignalValueMirror(signal, 11), 0.0);
    }

    private void reset() {
        startConvolution = finishConvolution = convolveProgressChange = 0;
    }

    @Override
    public void onStartConvolution() {
        startConvolution++;
    }

    @Override
    public void onFinishConvolution() {
        finishConvolution++;
    }

    @Override
    public void onConvolveProgressChange(final float progress) {
        convolveProgressChange++;
    }
}
