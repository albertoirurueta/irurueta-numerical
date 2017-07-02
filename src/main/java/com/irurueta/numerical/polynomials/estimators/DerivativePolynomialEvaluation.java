/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.polynomials.estimators.DerivativePolynomialEvaluation
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 6, 2016.
 */
package com.irurueta.numerical.polynomials.estimators;

/**
 * Contains an evaluation of the derivative of a given order of a polynomial and
 * the point where such derivative has been evaluated.
 */
public class DerivativePolynomialEvaluation extends PolynomialEvaluation {
    
    /**
     * Minimum allowed derivative order.
     */
    public static final int MIN_DERIVATIVE_ORDER = 1;
    
    /**
     * Point where derivative of a given order has been evaluated.
     */
    private double mX;
    
    /**
     * Order of derivative.
     */
    private int mDerivativeOrder;
    
    /**
     * Constructor.
     */
    public DerivativePolynomialEvaluation() {
        super();
        setDerivativeOrder(MIN_DERIVATIVE_ORDER);
    }
    
    /**
     * Constructor.
     * @param x point where derivative of polynomial has been evaluated.
     * @param evaluation evaluation of nth-derivative of polynomial at point x.
     * @param derivativeOrder order of derivative.
     * @throws IllegalArgumentException if order of derivative is less than 1.
     */
    public DerivativePolynomialEvaluation(double x, double evaluation, 
            int derivativeOrder) throws IllegalArgumentException {
        super(evaluation);
        mX = x;
        setDerivativeOrder(derivativeOrder);
    }
    
    /**
     * @param x point where derivative of polynomial has been evaluated.
     * @param evaluation evaluation of nth-derivative of polynomial at point x.
     */
    public DerivativePolynomialEvaluation(double x, double evaluation) {
        this(x, evaluation, MIN_DERIVATIVE_ORDER);
    }
    
    /**
     * Gets point where polynomial derivative has been evaluated.
     * @return point where polynomial derivative has been evaluated.
     */
    public double getX() {
        return mX;
    }
    
    /**
     * Sets point where polynomial derivative has been evaluated.
     * @param x point where polynomial derivative has been evaluated.
     */
    public void setX(double x) {
        mX = x;
    }
    
    /**
     * Gets order of derivative.
     * @return order of derivative.
     */
    public int getDerivativeOrder() {
        return mDerivativeOrder;
    }
    
    /**
     * Sets order of derivative.
     * @param derivativeOrder order of derivative.
     * @throws IllegalArgumentException if order of derivative is less than 1.
     */
    public final void setDerivativeOrder(int derivativeOrder) 
            throws IllegalArgumentException {
        if(derivativeOrder < MIN_DERIVATIVE_ORDER) {
            throw new IllegalArgumentException(
                    "derivative order must be at least 1");
        }
        
        mDerivativeOrder = derivativeOrder;
    }

    /**
     * Gets type of polynomial evaluation.
     * @return type of polynomial evaluation.
     */
    @Override
    public PolynomialEvaluationType getType() {
        return PolynomialEvaluationType.DERIVATIVE_EVALUATION;
    }
}
