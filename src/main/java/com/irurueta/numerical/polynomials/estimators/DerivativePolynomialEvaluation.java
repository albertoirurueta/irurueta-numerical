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
    private double x;

    /**
     * Order of derivative.
     */
    private int derivativeOrder;

    /**
     * Constructor.
     */
    public DerivativePolynomialEvaluation() {
        super();
        setDerivativeOrder(MIN_DERIVATIVE_ORDER);
    }

    /**
     * Constructor.
     *
     * @param x               point where derivative of polynomial has been evaluated.
     * @param evaluation      evaluation of nth-derivative of polynomial at point x.
     * @param derivativeOrder order of derivative.
     * @throws IllegalArgumentException if order of derivative is less than 1.
     */
    public DerivativePolynomialEvaluation(final double x, final double evaluation, final int derivativeOrder) {
        super(evaluation);
        this.x = x;
        setDerivativeOrder(derivativeOrder);
    }

    /**
     * @param x          point where derivative of polynomial has been evaluated.
     * @param evaluation evaluation of nth-derivative of polynomial at point x.
     */
    public DerivativePolynomialEvaluation(final double x, final double evaluation) {
        this(x, evaluation, MIN_DERIVATIVE_ORDER);
    }

    /**
     * Gets point where polynomial derivative has been evaluated.
     *
     * @return point where polynomial derivative has been evaluated.
     */
    public double getX() {
        return x;
    }

    /**
     * Sets point where polynomial derivative has been evaluated.
     *
     * @param x point where polynomial derivative has been evaluated.
     */
    public void setX(final double x) {
        this.x = x;
    }

    /**
     * Gets order of derivative.
     *
     * @return order of derivative.
     */
    public int getDerivativeOrder() {
        return derivativeOrder;
    }

    /**
     * Sets order of derivative.
     *
     * @param derivativeOrder order of derivative.
     * @throws IllegalArgumentException if order of derivative is less than 1.
     */
    public final void setDerivativeOrder(final int derivativeOrder) {
        if (derivativeOrder < MIN_DERIVATIVE_ORDER) {
            throw new IllegalArgumentException("derivative order must be at least 1");
        }

        this.derivativeOrder = derivativeOrder;
    }

    /**
     * Gets type of polynomial evaluation.
     *
     * @return type of polynomial evaluation.
     */
    @Override
    public PolynomialEvaluationType getType() {
        return PolynomialEvaluationType.DERIVATIVE_EVALUATION;
    }
}
