/**
 * This library contains packages for:
 * 
 * -numerical algorithms for function optimization (i.e.finding minima/maxima). 
 * Support for unidimensional functions is given using 
 * any of the available SingleOptimizer, for multidimensional functions an
 * implementation of MultiOptimizer must be used instead.
 * 
 * - classes to find function roots. Any implementation
 * of RootEstimator can be used for that purpose, which can be a 
 * SingleRootEstimator for any unidimensional function, or a 
 * PolynomialRootsEstimator if the function is polynomial.
 * 
 * - robust estimators to discard outliers (RANSAC, LMedS, PROSAC, etc)
 * 
 * - classes for model fitting in order to find the parameters of a given model.
 * Both unidimensional and multidimensional implementations exist for both 
 * linear and non linear models
 */
package com.irurueta.numerical;
