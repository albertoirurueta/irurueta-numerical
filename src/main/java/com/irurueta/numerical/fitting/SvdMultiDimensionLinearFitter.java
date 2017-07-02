/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fitting.SvdMultiDimensionLinearFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 27, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.SingularValueDecomposer;
import com.irurueta.numerical.NotReadyException;

/**
 * Fits provided data (x,y) to a function made of a linear combination of
 * functions used as a basis (i.e. f(x1, x2, ...) = a * f0(x1, x2, ...) + 
 * b * f1(x1, x2, ...) + ...).
 * Where f0, f1, ... is the function basis which ideally should be formed by
 * orthogonal function.
 * This class is based on the implementation available at Numerical Recipes
 * 3rd Ed, page 795
 */
public class SvdMultiDimensionLinearFitter extends MultiDimensionLinearFitter{
    
    /**
     * Default tolerance
     */
    public static final double TOL = 1e-12;
    
    /**
     * Tolerance to define convergence threshold for SVD
     */
    private double tol;

    /**
     * Constructor
     */
    SvdMultiDimensionLinearFitter(){
        super();
        tol = TOL;        
    }
    
    /**
     * Constructor
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviations of each pair of points (x, y)
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */
    public SvdMultiDimensionLinearFitter(Matrix x, double[] y, double[] sig)
            throws IllegalArgumentException{
        super(x, y, sig);
        tol = TOL;
    }

    /**
     * Constructor
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */
    public SvdMultiDimensionLinearFitter(Matrix x, double[] y, double sig)
            throws IllegalArgumentException{
        super(x, y, sig);
        tol = TOL;
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and
     * obtain the evaluation of function basis at such point
     * @throws FittingException if evaluation fails
     */
    public SvdMultiDimensionLinearFitter(LinearFitterMultiDimensionFunctionEvaluator evaluator)
            throws FittingException{
        super(evaluator);
        tol = TOL;
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and
     * obtain the evaluation of function basis at such point
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviations of each pair of points (x, y)
     * @throws FittingException if evaluation fails
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */
    public SvdMultiDimensionLinearFitter(LinearFitterMultiDimensionFunctionEvaluator evaluator,
            Matrix x, double[] y, double[] sig)
            throws FittingException, IllegalArgumentException{
        super(evaluator, x, y, sig); 
        tol = TOL;
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and
     * obtain the evaluation of function basis at such point
     * @param x input points x where a linear multi dimensional function 
     * f(x1, x2, ...) = a * f0(x1, x2, ...) + b * f1(x1, x2, ...) + ...
     * @param y result of evaluation of linear multi dimensional function 
     * f(x1, x2, ...) at provided x points
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant
     * @throws FittingException if evaluation fails
     * @throws IllegalArgumentException if provided matrix rows and arrays 
     * don't have the same length
     */    
    public SvdMultiDimensionLinearFitter(LinearFitterMultiDimensionFunctionEvaluator evaluator,
            Matrix x, double[] y, double sig)
            throws FittingException, IllegalArgumentException{
        super(evaluator, x, y, sig);
        tol = TOL;
    }     
    
    /**
     * Returns tolerance to define convergence threshold for SVD
     * @return tolerance to define convergence threshold for SVD
     */
    public double getTol(){
        return tol;
    }
    
    /**
     * Sets tolerance to define convergence threshold for SVD
     * @param tol tolerance to define convergence threshold for SVD
     */
    public void setTol(double tol){
        this.tol = tol;
    }    
    
    /**
     * Fits a function to provided data so that parameters associated to that
     * function can be estimated along with their covariance matrix and chi
     * square value
     * @throws FittingException if fitting fails
     * @throws NotReadyException if enough input data has not yet been provided
     */    
    @Override
    public void fit() throws FittingException, NotReadyException {
        if(!isReady()) throw new NotReadyException();
        
        double[] xRow = new double[x.getColumns()];
        int xCols = evaluator.getNumberOfDimensions();
        
        try{
            resultAvailable = false;
            
            int i,j,k;
            double tmp,thresh,sum;
            Matrix aa = new Matrix(ndat, ma);
            double[] b = new double[ndat];
            for(i = 0; i < ndat; i++){
                x.getSubmatrixAsArray(i, 0, i, xCols - 1, xRow);
                evaluator.evaluate(xRow, afunc);
		tmp=1.0/sig[i];
		for(j = 0; j < ma; j++){
                    aa.setElementAt(i, j, afunc[j]*tmp);
                }
		b[i]=y[i]*tmp;
            }
            
            SingularValueDecomposer svd = 
                    new SingularValueDecomposer(aa);
            svd.decompose();
            thresh = (tol > 0. ? tol * svd.getSingularValues()[0] : -1.);
            svd.solve(b, thresh, a);
            chisq = 0.0;
            for(i = 0; i < ndat; i++){
		sum=0.;
		for(j = 0; j < ma; j++){
                    sum += aa.getElementAt(i, j)*a[j];
                }
		chisq += Math.pow(sum-b[i], 2.0);
            }
            for(i = 0; i < ma; i++){
                for(j = 0; j < i + 1; j++){
                    sum=0.0;
                    double[] w = svd.getSingularValues();
                    double tsh = svd.getNegligibleSingularValueThreshold();
                    Matrix v = svd.getV();
                    for (k=0;k<ma;k++){
                        if (w[k] > tsh){
                            sum += v.getElementAt(i, k) * v.getElementAt(j, k) / 
                                    Math.pow(w[k], 2.0);
                        }
                    }
                    covar.setElementAt(j, i, sum);
                    covar.setElementAt(i, j, sum);
		}
            }
            
            resultAvailable = true;
        }catch(FittingException e){
            throw e;
        }catch(Throwable t){
            throw new FittingException(t);
        }        
    }    
    
}
