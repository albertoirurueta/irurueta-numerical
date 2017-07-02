/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.PolynomialEvaluation
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016
 */
package com.irurueta.numerical.polynomials.estimators;

import java.io.Serializable;

/**
 * Contains an evaluation of a polynomial and the point where the
 * polynomial has been evaluated.
 */
public abstract class PolynomialEvaluation implements Serializable{    
    /**
     * Evaluation of polynomial at point x.
     */
    private double mEvaluation;
    
    /**
     * Constructor.
     */
    public PolynomialEvaluation() {}
    
    /**
     * Constructor.
     * @param evaluation evaluation of polynomial at a point x.
     */
    public PolynomialEvaluation(double evaluation) {
        mEvaluation = evaluation;
    }
    
    /**
     * Constructor.
     * @param x point where polynomial has been evaluated.
     * @param evaluation evaluation of polynomial at point x.
     */
    public PolynomialEvaluation(double x, double evaluation) {
        mEvaluation = evaluation;
    }
        
    /**
     * Gets evaluation of polynomial at point x.
     * @return evaluation of polynomial at point x.
     */
    public double getEvaluation() {
        return mEvaluation;
    }
    
    /**
     * Sets evaluation of polynomial at point x.
     * @param evaluation evaluation of polynomial at point x.
     */
    public void setEvaluation(double evaluation) {
        mEvaluation = evaluation;
    }
    
    /**
     * Gets type of polynomial evaluation.
     * @return type of polynomial evaluation.
     */
    public abstract PolynomialEvaluationType getType();
}
