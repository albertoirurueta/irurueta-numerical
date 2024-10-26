/*
 * Copyright (C) 2023 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.integration;

import com.irurueta.algebra.WrongSizeException;

/**
 * Computes function integration by using Simpson's rule an Upper Square Root mid-point quadrature
 * when upper integration bound lies at a function singularity.
 */
public class SimpsonUpperSquareRootMidPointQuadratureMatrixIntegrator
        extends SimpsonMatrixIntegrator<UpperSquareRootMidPointMatrixQuadrature> {

   /**
    * Constructor.
    *
    * @param a        Lower limit of integration.
    * @param b        Upper limit of integration.
    * @param listener listener to evaluate a single dimension matrix (multivariate) function at
    *                 required points.
    * @param eps      required accuracy.
    * @throws WrongSizeException if size notified by provided listener is invalid.
    */
   public SimpsonUpperSquareRootMidPointQuadratureMatrixIntegrator(
           final double a, final double b, final MatrixSingleDimensionFunctionEvaluatorListener listener,
           final double eps) throws WrongSizeException {
      super(new UpperSquareRootMidPointMatrixQuadrature(a, b, listener), eps);
   }

   /**
    * Constructor with default accuracy.
    *
    * @param a        Lower limit of integration.
    * @param b        Upper limit of integration.
    * @param listener listener to evaluate a single dimension function at required points.
    * @throws WrongSizeException if size notified by provided listener is invalid.
    */
   public SimpsonUpperSquareRootMidPointQuadratureMatrixIntegrator(
           final double a, final double b, final MatrixSingleDimensionFunctionEvaluatorListener listener)
           throws WrongSizeException {
      this(a, b, listener, EPS);
   }

   /**
    * Gets type of quadrature.
    *
    * @return type of quadrature.
    */
   @Override
   public QuadratureType getQuadratureType() {
      return QuadratureType.UPPER_SQUARE_ROOT_MID_POINT;
   }
}
