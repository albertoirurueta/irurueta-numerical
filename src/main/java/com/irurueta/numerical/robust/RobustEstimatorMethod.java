/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.RobustEstimatorMethod
 */
package com.irurueta.numerical.robust;

/**
 * Enumerator containing different robust estimation algorithms
 */
public enum RobustEstimatorMethod {
    /**
     * Random Sample Consensus
     */
    RANSAC,
    
    /**
     * Least Median of Squares
     */
    LMedS,    
    
    /**
     * Median Sample Consensus
     */
    MSAC,
    
    /**
     * Progressive Sample Consensus
     */
    PROSAC,
    
    /**
     * Progressive Median of Squares
     */
    PROMedS
}
