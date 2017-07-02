/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.signal.processing.ConvolverEdgeMethod
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
