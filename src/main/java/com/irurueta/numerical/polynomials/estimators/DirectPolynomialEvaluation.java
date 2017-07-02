/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.DirectPolynomialEvaluation
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

/**
 * Contains an evaluation of a polynomial and the point where the polynomial
 * has been evaluated.
 */
public class DirectPolynomialEvaluation extends PolynomialEvaluation {
    /**
     * Point where polynomial has been evaluated.
     */
    private double mX;    
    
    /**
     * Constructor.
     */
    public DirectPolynomialEvaluation() {
        super();
    }
    
    /**
     * Constructor.
     * @param x point where polynomial has been evaluated.
     * @param evaluation evaluation of polynomial at point x.
     */
    public DirectPolynomialEvaluation(double x, double evaluation) {
        super(evaluation);
        mX = x;
    }
    
    /**
     * Gets point where polynomial has been evaluated.
     * @return point where polynomial has been evaluated.
     */
    public double getX() {
        return mX;
    }
    
    /**
     * Sets point where polynomial has been evaluated.
     * @param x point where polynomial has been evaluated.
     */
    public void setX(double x) {
        mX = x;
    }    
    
    /**
     * Gets type of polynomial evaluation.
     * @return type of polynomial evaluation.
     */
    @Override
    public PolynomialEvaluationType getType() {
        return PolynomialEvaluationType.DIRECT_EVALUATION;
    }
}
