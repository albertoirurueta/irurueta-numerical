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
