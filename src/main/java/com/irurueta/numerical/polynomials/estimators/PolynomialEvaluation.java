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
    public PolynomialEvaluation() { }
    
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
