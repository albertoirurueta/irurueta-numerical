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

/**
 * This enumerator indicates how edges should be treated during convolution.
 * Repeat and mirror edge extension methods can be used to prevent certain
 * artifacts when reaching edges of signal.
 */
public enum ConvolverEdgeMethod {
    /**
     * When convolution kernel reaches edge of signal being convoluted, it is
     * assumed that the signal value is zero.
     */
    ZERO_EDGE,

    /**
     * When convolution kernel reaches edge of signal being convoluted, it is
     * assumed that the signal has a constant value (with a value that can be
     * setup).
     */
    CONSTANT_EDGE,

    /**
     * When convolution kernel reaches edge of signal being convoluted, it is
     * assumed that the signal is repeated indefinitely.
     */
    REPEAT_EDGE,

    /**
     * When convolution kernel reaches edge of signal being convoluted, it is
     * assumed that the signal is mirrored.
     */
    MIRROR_EDGE
}
