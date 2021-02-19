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
 * Contains an evaluation of the nth-integral of a polynomial and the
 * point where such integral has been evaluated.
 */
@SuppressWarnings("WeakerAccess")
public class IntegralPolynomialEvaluation extends PolynomialEvaluation {

    /**
     * Minimum allowed integral order.
     */
    public static final int MIN_INTEGRAL_ORDER = 1;

    /**
     * Point where integral of a given order has been evaluated.
     */
    private double mX;

    /**
     * Constant terms of integral.
     */
    private double[] mConstants;

    /**
     * Order of integral.
     */
    private int mIntegralOrder;

    /**
     * Constructor.
     */
    public IntegralPolynomialEvaluation() {
        super();
        setIntegralOrder(MIN_INTEGRAL_ORDER);
    }

    /**
     * Constructor.
     *
     * @param x             point where integral of polynomial has been evaluated.
     * @param evaluation    evaluation of nth-integral of polynomial at point x.
     * @param constants     constant terms of nth-integral.
     * @param integralOrder order of integral.
     * @throws IllegalArgumentException if order of integral is less than 1.
     */
    public IntegralPolynomialEvaluation(
            final double x, final double evaluation,
            final double[] constants, final int integralOrder) {
        super(evaluation);
        mX = x;
        mConstants = constants;
        setIntegralOrder(integralOrder);
    }

    /**
     * Constructor.
     *
     * @param x          point where integral of polynomial has been evaluated.
     * @param evaluation evaluation of nth-integral of polynomial at point x.
     * @param constants  constant terms of nth-integral.
     */
    public IntegralPolynomialEvaluation(
            final double x, final double evaluation,
            final double[] constants) {
        this(x, evaluation, constants, MIN_INTEGRAL_ORDER);
    }

    /**
     * Constructor.
     *
     * @param x             point where integral of polynomial has been evaluated.
     * @param evaluation    evaluation of nth-integral of polynomial at point x.
     * @param integralOrder order of integral.
     * @throws IllegalArgumentException if order of integral is less than 1.
     */
    public IntegralPolynomialEvaluation(
            final double x, final double evaluation,
            final int integralOrder) {
        this(x, evaluation, null, integralOrder);
    }

    /**
     * Constructor.
     *
     * @param x          point where integral of polynomial has been evaluated.
     * @param evaluation evaluation of nth-integral of polynomial at point x.
     */
    public IntegralPolynomialEvaluation(
            final double x, final double evaluation) {
        this(x, evaluation, null, MIN_INTEGRAL_ORDER);
    }

    /**
     * Gets point where nth-polynomial integral has been evaluated.
     *
     * @return point where nth-polynomial integral has been evaluated.
     */
    public double getX() {
        return mX;
    }

    /**
     * Sets point where nth-polynomial integral has been evaluated.
     *
     * @param x point where nth-polynomial integral has been evaluated.
     */
    public void setX(final double x) {
        mX = x;
    }

    /**
     * Gets constant terms of integral.
     *
     * @return constant terms of integral.
     */
    public double[] getConstants() {
        return mConstants;
    }

    /**
     * Sets constant terms of integral.
     *
     * @param constants constant terms of integral.
     */
    public void setConstant(final double[] constants) {
        mConstants = constants;
    }

    /**
     * Gets integral order.
     *
     * @return integral order.
     */
    public int getIntegralOrder() {
        return mIntegralOrder;
    }

    /**
     * Sets integral order.
     *
     * @param integralOrder integral order.
     * @throws IllegalArgumentException if integral order is less than 1.
     */
    public final void setIntegralOrder(final int integralOrder) {
        if (integralOrder < MIN_INTEGRAL_ORDER) {
            throw new IllegalArgumentException(
                    "integral order must be at least 1");
        }

        mIntegralOrder = integralOrder;
    }

    /**
     * Gets type of polynomial evaluation.
     *
     * @return type of polynomial evaluation.
     */
    @Override
    public PolynomialEvaluationType getType() {
        return PolynomialEvaluationType.INTEGRAL_EVALUATION;
    }
}
