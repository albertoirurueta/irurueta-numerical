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
public class IntegralIntervalPolynomialEvaluation extends PolynomialEvaluation {

    /**
     * Minimum allowed integral order.
     */
    public static final int MIN_INTEGRAL_ORDER = 1;

    /**
     * Start point of interval being integrated.
     */
    private double startX;

    /**
     * End point of interval being integrated.
     */
    private double endX;

    /**
     * Order of integral.
     */
    private int integralOrder;

    /**
     * Constant terms of integral.
     */
    private double[] constants;

    /**
     * Constructor.
     */
    public IntegralIntervalPolynomialEvaluation() {
        super();
        setIntegralOrder(MIN_INTEGRAL_ORDER);
    }

    /**
     * Constructor.
     *
     * @param startX        start point of interval being integrated.
     * @param endX          end point of interval being integrated.
     * @param evaluation    evaluation of nth-integral of polynomial between startX
     *                      and endX.
     * @param integralOrder order of integral.
     * @throws IllegalArgumentException if order of integral is less than 1.
     */
    public IntegralIntervalPolynomialEvaluation(
            final double startX, final double endX, final double evaluation, final int integralOrder) {
        super(evaluation);
        this.startX = startX;
        this.endX = endX;
        setIntegralOrder(integralOrder);
    }

    /**
     * Constructor.
     *
     * @param startX     start point of interval being integrated.
     * @param endX       end point of interval being integrated.
     * @param evaluation evaluation of nth-integral of polynomial between startX
     *                   and endX.
     */
    public IntegralIntervalPolynomialEvaluation(final double startX, final double endX, final double evaluation) {
        this(startX, endX, evaluation, MIN_INTEGRAL_ORDER);
    }

    /**
     * Gets start point of interval being integrated.
     *
     * @return start point of interval being integrated.
     */
    public double getStartX() {
        return startX;
    }

    /**
     * Sets start point of interval being integrated.
     *
     * @param startX start point of interval being integrated.
     */
    public void setStartX(final double startX) {
        this.startX = startX;
    }

    /**
     * Gets end point of interval being integrated.
     *
     * @return end point of interval being integrated.
     */
    public double getEndX() {
        return endX;
    }

    /**
     * Sets end point of interval being integrated.
     *
     * @param endX end point of interval being integrated.
     */
    public void setEndX(final double endX) {
        this.endX = endX;
    }

    /**
     * Gets integral order.
     *
     * @return integral order.
     */
    public int getIntegralOrder() {
        return integralOrder;
    }

    /**
     * Sets integral order.
     *
     * @param integralOrder integral order.
     * @throws IllegalArgumentException if integral order is less than 1.
     */
    public final void setIntegralOrder(final int integralOrder) {
        if (integralOrder < MIN_INTEGRAL_ORDER) {
            throw new IllegalArgumentException("integral order must be at least 1");
        }

        this.integralOrder = integralOrder;
    }

    /**
     * Gets constant terms of integral.
     *
     * @return constant terms of integral.
     */
    public double[] getConstants() {
        return constants;
    }

    /**
     * Sets constant terms of integral.
     *
     * @param constants constant terms of integral.
     */
    public void setConstants(final double[] constants) {
        this.constants = constants;
    }

    /**
     * Gets type of polynomial evaluation.
     *
     * @return type of polynomial evaluation.
     */
    @Override
    public PolynomialEvaluationType getType() {
        return PolynomialEvaluationType.INTEGRAL_INTERVAL;
    }
}
