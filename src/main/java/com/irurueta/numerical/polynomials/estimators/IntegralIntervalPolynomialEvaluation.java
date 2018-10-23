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
package com.irurueta.numerical.polynomials.estimators;

/**
 * Contains an evaluation of an interval of the nth-integral of a polynomial.
 */
@SuppressWarnings("WeakerAccess")
public class IntegralIntervalPolynomialEvaluation extends PolynomialEvaluation {
    
    /**
     * Minimum allowed integral order.
     */
    public static final int MIN_INTEGRAL_ORDER = 1;
    
    /**
     * Start point of interval being integrated.
     */
    private double mStartX;
    
    /**
     * End point of interval being integrated.
     */
    private double mEndX;
    
    /**
     * Order of integral.
     */
    private int mIntegralOrder;
    
    /**
     * Constant terms of integral.
     */
    private double[] mConstants;
    
    /**
     * Constructor.
     */
    public IntegralIntervalPolynomialEvaluation() {
        super();
        setIntegralOrder(MIN_INTEGRAL_ORDER);
    }
    
    /**
     * Constructor.
     * @param startX start point of interval being integrated.
     * @param endX end point of interval being integrated.
     * @param evaluation evaluation of nth-integral of polynomial between startX 
     * and endX.
     * @param integralOrder order of integral.
     * @throws IllegalArgumentException if order of integral is less than 1.
     */
    public IntegralIntervalPolynomialEvaluation(double startX, double endX, 
            double evaluation, int integralOrder) 
            throws IllegalArgumentException {
        super(evaluation);
        mStartX = startX;
        mEndX = endX;
        setIntegralOrder(integralOrder);
    }
    
    /**
     * Constructor.
     * @param startX start point of interval being integrated.
     * @param endX end point of interval being integrated.
     * @param evaluation evaluation of nth-integral of polynomial between startX
     * and endX.
     */
    public IntegralIntervalPolynomialEvaluation(double startX, double endX,
            double evaluation) {
        this(startX, endX, evaluation, MIN_INTEGRAL_ORDER);
    }
    
    /**
     * Gets start point of interval being integrated.
     * @return start point of interval being integrated.
     */
    public double getStartX() {
        return mStartX;
    }
    
    /**
     * Sets start point of interval being integrated.
     * @param startX start point of interval being integrated.
     */
    public void setStartX(double startX) {
        mStartX = startX;
    }

    /**
     * Gets end point of interval being integrated.
     * @return end point of interval being integrated.
     */
    public double getEndX() {
        return mEndX;
    }
    
    /**
     * Sets end point of interval being integrated.
     * @param endX end point of interval being integrated.
     */
    public void setEndX(double endX) {
        mEndX = endX;
    }

    /**
     * Gets integral order.
     * @return integral order.
     */
    public int getIntegralOrder() {
        return mIntegralOrder;
    }
    
    /**
     * Sets integral order.
     * @param integralOrder integral order.
     * @throws IllegalArgumentException if integral order is less than 1.
     */
    public final void setIntegralOrder(int integralOrder)
            throws IllegalArgumentException {
        if(integralOrder < MIN_INTEGRAL_ORDER) {
            throw new IllegalArgumentException(
                    "integral order must be at least 1");
        }
        
        mIntegralOrder = integralOrder;
    }

    /**
     * Gets constant terms of integral.
     * @return constant terms of integral.
     */
    public double[] getConstants() {
        return mConstants;
    }
    
    /**
     * Sets constant terms of integral.
     * @param constants constant terms of integral.
     */
    public void setConstants(double[] constants) {
        mConstants = constants;
    }    
        
    /**
     * Gets type of polynomial evaluation.
     * @return type of polynomial evaluation.
     */
    @Override
    public PolynomialEvaluationType getType() {
        return PolynomialEvaluationType.INTEGRAL_INTERVAL;
    }    
}
