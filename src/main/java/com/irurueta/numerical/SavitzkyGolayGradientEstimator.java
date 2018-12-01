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

import com.irurueta.algebra.*;

/**
 * Class to estimate the gradient of a multidimensional function.
 * This class evaluates a function at very close locations of a given input
 * point in order to determine the gradient at such point
 * The algorithm used in this implementation is valid for continuous functions
 * only, otherwise inaccurate results might be obtain.
 * This implementation is more robust against small discontinuities than 
 * SymmetricGradientEstimator, but it is also slower to compute.
 * This method interpolates the sampled function values into a polynomial of
 * 2nd degree (parabolic) on each direction, whose derivative is known.
 * Because a linear system of equations has to be solved to determine such
 * polynomials, this method might be less accurate when large values are 
 * involved due to limited machine precision.
 */
@SuppressWarnings("WeakerAccess")
public class SavitzkyGolayGradientEstimator extends GradientEstimator {
    /**
     * Number of required point to evaluate to compute derivative.
     */    
    public static final int N_POINTS = 3;
    
    /**
     * Constructor.
     * @param listener Listener to evaluate a multidimensional functin.
     */    
    public SavitzkyGolayGradientEstimator(
            MultiDimensionFunctionEvaluatorListener listener) {
        super(listener);
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
    @Override
    @SuppressWarnings("Duplicates")
    public void gradient(double[] point, double[] result) 
            throws EvaluationException {
        int n = point.length;
        if (result.length != n) {
            throw new IllegalArgumentException();
        }
        
        double[] xh1 = new double[n];
        double[] xh2 = new double[n];
        
        double f = listener.evaluate(point);
        System.arraycopy(point, 0, xh1, 0, n);
        System.arraycopy(point, 0, xh2, 0, n);
        
        Matrix a;
        try {
            a = new Matrix(N_POINTS, N_POINTS);
        } catch (WrongSizeException e) {
            throw new EvaluationException(e);
        }

        double[] b = new double[N_POINTS];
        
        SingularValueDecomposer decomposer = new SingularValueDecomposer(a);
        
        for (int j = 0; j < n; j++) {
            double temp = point[j];
            double h = EPS * Math.abs(temp);
            if (h == 0.0) {
                h = EPS; //Trick to reduce finite-precision error
            }
            
            double p1 = xh1[j] = temp + h;
            double p2 = xh2[j] = temp - h;
            
            double fh1 = listener.evaluate(xh1);
            double fh2 = listener.evaluate(xh2);

            xh1[j] = temp;
            xh2[j] = temp;

            a.setElementAt(0, 0, temp * temp);
            a.setElementAt(1, 0, p1 * p1);
            a.setElementAt(2, 0, p2 * p2);

            a.setElementAt(0, 1, temp);
            a.setElementAt(1, 1, p1);
            a.setElementAt(2, 1, p2);

            a.setElementAt(0, 2, 1.0);
            a.setElementAt(1, 2, 1.0);
            a.setElementAt(2, 2, 1.0);

            //normalize to increase accuracy
            double normA = Utils.normF(a);
            a.multiplyByScalar(1.0 / normA);

            b[0] = f;
            b[1] = fh1;
            b[2] = fh2;

            //normalize to increase accuracy
            ArrayUtils.multiplyByScalar(b, 1.0 / normA, b);

            double aParam;
            double bParam;
            try {
                decomposer.setInputMatrix(a);
                decomposer.decompose();

                //now solve the system of equations in Least Mean Squared Error
                //because SVD allows the system of equations to be solved using 
                //the pseudo-inverse                
                double[] params = decomposer.solve(b);
                aParam = params[0];
                bParam = params[1];
                //and c = params[2], but we don't need it
                
                //because we have fitted the function in dimension j into a
                //polynomial that has expression: a * x^2 + b * x + c, then its
                //partial derivative on dimension j is:
                //2.0 * a * x + b , therefore:
                result[j] = 2.0 * aParam * temp + bParam;
            } catch (AlgebraException e) {
                result[j] = Double.NaN;
            }                        
        }
    }
}
