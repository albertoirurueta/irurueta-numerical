/*
 * Copyright (C) 2012 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical;

/**
 * Class to estimate the gradient of a multidimensional function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the gradient at such point.
 */
public class GradientEstimator {
    
    /**
     * Constant considered as machine precision.
     */
    public static final double EPS = 1e-8;

    /**
     * Listener to evaluate a multidimensional function.
     */
    MultiDimensionFunctionEvaluatorListener listener;
    
    /**
     * Internal array to hold input parameter's values.
     */
    private double[] xh;
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional function.
     */
    public GradientEstimator(MultiDimensionFunctionEvaluatorListener listener) {
        this.listener = listener;
    }
    
    /**
     * Returns the gradient of a multidimensional function at provided point.
     * @param point Input point.
     * @return Gradient.
     * @throws EvaluationException Raised if function cannot be evaluated.
     */
    public double[] gradient(double[] point) throws EvaluationException {
        double[] result = new double[point.length];
        gradient(point, result);
        return result;
    }
    
    /**
     * Sets estimated gradient in provided result array of a multidimensional
     * function at provided point.
     * This method is preferred respect to gradient(double[]) because result
     * array can be reused and hence is more memory efficient.
     * @param point Input point.
     * @param result Output parameter containing estimated array. This parameter
     * must be an array of length equal to point.
     * @throws EvaluationException Raised if function cannot be evaluated.
     * @throws IllegalArgumentException Raised if length of result and point are
     * not equal.
     */
    @SuppressWarnings("Duplicates")
    public void gradient(double[] point, double[] result) 
            throws EvaluationException {
        int length = point.length;
        if (result.length != length) {
            throw new IllegalArgumentException();
        }
        
        if (xh == null || xh.length != length) {
            xh = new double[length];
        }
        System.arraycopy(point, 0, xh, 0, length);
        
        try {
            double temp;
            double h;
            double fh;
            double fold = listener.evaluate(point);
            for (int j = 0; j < length; j++) {
                temp = point[j];
                h = EPS * Math.abs(temp);
                if(h == 0.0) h = EPS;
                xh[j] = temp + h;
                h = xh[j] - temp;
                fh = listener.evaluate(xh);
                xh[j] = temp;
                result[j] = (fh - fold) / h;
            }
        } catch (EvaluationException e) {
            throw e;
        } catch (Exception e) {
            throw new EvaluationException(e);
        }
    }
}
