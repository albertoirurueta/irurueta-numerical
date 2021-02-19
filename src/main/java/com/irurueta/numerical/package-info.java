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

/**
 * This library contains packages for:
 * <p>
 * -numerical algorithms for function optimization (i.e.finding minima/maxima).
 * Support for unidimensional functions is given using
 * any of the available SingleOptimizer, for multidimensional functions an
 * implementation of MultiOptimizer must be used instead.
 * <p>
 * - classes to find function roots. Any implementation
 * of RootEstimator can be used for that purpose, which can be a
 * SingleRootEstimator for any unidimensional function, or a
 * PolynomialRootsEstimator if the function is polynomial.
 * <p>
 * - robust estimators to discard outliers (RANSAC, LMedS, PROSAC, etc)
 * <p>
 * - classes for model fitting in order to find the parameters of a given model.
 * Both unidimensional and multidimensional implementations exist for both
 * linear and non linear models
 */
package com.irurueta.numerical;
