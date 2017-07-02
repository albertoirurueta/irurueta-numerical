/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.fitting.LevenbergMarquardtSingleDimensionFitter
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date May 28, 2015
 */
package com.irurueta.numerical.fitting;

import com.irurueta.algebra.AlgebraException;
import com.irurueta.algebra.GaussJordanElimination;
import com.irurueta.algebra.Matrix;
import com.irurueta.numerical.NotReadyException;
import java.util.Arrays;

/**
 * Fits provided data (x,y) to a generic non-linear function using 
 * Levenberg-Marquardt iterative algorithm.
 * This class is based on the implementation available at Numerical Recipes 3rd
 * Ed, page 801
 */
public class LevenbergMarquardtSingleDimensionFitter 
        extends SingleDimensionFitter{
    
    /**
     * Default convergence parameter. Number of times that tolerance is assumed
     * to be reached to consider that algorithm has finished iterating
     */
    public static final int NDONE = 4;
    
    /**
     * Default maximum number of iterations
     */
    public static final int ITMAX = 1000;
    
    /**
     * Default tolerance to reach convergence
     */
    public static final double TOL = 1e-3;
    
    /**
     * Convergence parameter
     */
    private int ndone;
    
    /**
     * Maximum number of iterations
     */
    private int itmax;
    
    /**
     * Tolerance to reach convergence
     */
    private double tol;
    
    /**
     * Evaluator of functions
     */
    private LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator;
    
    /**
     * Number of function parameters to be estimated
     */
    private int ma;
    
    /**
     * Determines which parameters can be modified during estimation (if true)
     * and which ones are locked (if false)
     */
    private boolean[] ia;    
    
    /**
     * Curvature matrix
     */
    private Matrix alpha;
    
    /**
     * Number of parameters ot be fitted
     */
    private int mfit = 0;
        
    /**
     * Constructor
     */
    public LevenbergMarquardtSingleDimensionFitter(){
        super();
        ndone = NDONE;
        itmax = ITMAX;
        tol = TOL;
    }
    
    /**
     * Constructor
     * @param x input points x where function f(x) is evaluated
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points
     * @param sig standard deviations of each pair of points (x, y)
     * @throws IllegalArgumentException if provided arrays don't have the same
     * length
     */
    public LevenbergMarquardtSingleDimensionFitter(double[] x, double[] y, 
            double[] sig) throws IllegalArgumentException{
        super(x, y, sig);
        ndone = NDONE;
        itmax = ITMAX;
        tol = TOL;
    }
    
    /**
     * Constructor
     * @param x input points x where function f(x) is evaluated
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant
     * @throws IllegalArgumentException if provided arrays don't have the same 
     * length 
     */
    public LevenbergMarquardtSingleDimensionFitter(double[] x, double[] y, 
            double sig) throws IllegalArgumentException{
        super(x, y, sig);
        ndone = NDONE;
        itmax = ITMAX;
        tol = TOL;
    }

    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and 
     * obtain the evaluation of function basis at such point
     * @throws FittingException if evaluation fails
     */
    public LevenbergMarquardtSingleDimensionFitter(
            LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator)
            throws FittingException{
        this();
        setFunctionEvaluator(evaluator);
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and 
     * obtain the evaluation of function basis at such point
     * @param x input points x where function f(x) is evaluated
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points
     * @param sig standard deviations of each pair of points (x, y)
     * @throws IllegalArgumentException if provided arrays don't have the same
     * length
     * @throws FittingException if evaluation fails
     */
    public LevenbergMarquardtSingleDimensionFitter(
            LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator,
            double[] x, double[] y, double[] sig) 
            throws IllegalArgumentException, FittingException{
        this(x, y, sig);
        setFunctionEvaluator(evaluator);
    }
    
    /**
     * Constructor
     * @param evaluator evaluator to evaluate function at provided point and 
     * obtain the evaluation of function basis at such point
     * @param x input points x where function f(x) is evaluated
     * @param y result of evaluation of linear single dimensional function f(x)
     * at provided x points
     * @param sig standard deviation of all pair of points assuming that 
     * standard deviations are constant
     * @throws IllegalArgumentException if provided arrays don't have the same 
     * length 
     * @throws FittingException if evaluation fails
     */
    public LevenbergMarquardtSingleDimensionFitter(
            LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator, 
            double[] x, double[] y, double sig) throws IllegalArgumentException,
            FittingException{
        this(x, y, sig);
        setFunctionEvaluator(evaluator);
    }    
    
    /**
     * Returns convergence parameter
     * @return convergence parameter
     */
    public int getNdone(){
        return ndone;
    }
    
    /**
     * Sets convergence parameter
     * @param ndone convergence parameter
     * @throws IllegalArgumentException if provided value is less than 1
     */
    public void setNdone(int ndone) throws IllegalArgumentException{
        if(ndone < 1) throw new IllegalArgumentException();
        this.ndone = ndone;
    }

    /**
     * Returns maximum number of iterations
     * @return maximum number of iterations
     */
    public int getItmax(){
        return itmax;
    }
    
    /**
     * Sets maximum number of iterations
     * @param itmax maximum number of iterations
     * @throws IllegalArgumentException if provided value is zero or negative
     */
    public void setItmax(int itmax) throws IllegalArgumentException{
        if(itmax <= 0) throw new IllegalArgumentException();
        this.itmax = itmax;
    }    

    /**
     * Returns tolerance to reach convergence
     * @return tolerance to reach convergence
     */
    public double getTol(){
        return tol;
    }
    
    /**
     * Sets tolerance to reach convergence
     * @param tol tolerance to reach convergence
     * @throws IllegalArgumentException if provided value is zero or negative
     */
    public void setTol(double tol) throws IllegalArgumentException{
        if(tol <= 0.0) throw new IllegalArgumentException();
        this.tol = tol;
    }
    
    /**
     * Returns function evaluator to evaluate function at a given point and
     * obtain function derivatives respect to each provided parameter
     * @return function evaluator
     */
    public LevenbergMarquardtSingleDimensionFunctionEvaluator 
            getFunctionEvaluator(){
        return evaluator;
    }
    
    /**
     * Sets function evaluator to evaluate function at a given point and obtain
     * function derivatives respect to each provided parameter
     * @param evaluator function evaluator
     * @throws FittingException if evaluation fails
     */
    public final void setFunctionEvaluator(
            LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator) 
            throws FittingException{
        internalSetFunctionEvaluator(evaluator);
    }
    
    /**
     * Internal method to set function evaluator to evaluate function at a given 
     * point and obtain function derivatives respect to each provided parameter
     * @param evaluator function evaluator
     * @throws FittingException if evaluation fails
     */    
    private void internalSetFunctionEvaluator(
            LevenbergMarquardtSingleDimensionFunctionEvaluator evaluator) 
            throws FittingException{
        
        try{
            this.evaluator = evaluator;    
        
            if(evaluator != null){
                a = evaluator.createInitialParametersArray();
                ma = a.length;
                covar = new Matrix(ma, ma);
                alpha = new Matrix(ma, ma);
                ia = new boolean[ma];
                Arrays.fill(ia, true);                
            }            
        }catch(AlgebraException e){
            throw new FittingException(e);
        }
    }

    /**
     * Indicates whether provided instance has enough data to start the function
     * fitting.
     * @return true if this instance is ready to start the function fitting, 
     * false otherwise
     */    
    @Override
    public boolean isReady() {
        return evaluator != null && x != null && y != null && 
                x.length == y.length;
    }
    
    /**
     * Returns curvature matrix
     * @return curvature matrix
     */
    public Matrix getAlpha(){
        return alpha;
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

        try{
            resultAvailable = false;        
        
            int j, k, l, iter, done = 0;
            double alamda = .001, ochisq;
            double[] atry = new double[ma];
            double[] beta = new double[ma];
            double[] da = new double[ma];
        
            mfit=0; //number of parameters to be fitted
            for(j = 0; j < ma; j++) if(ia[j]) mfit++;
            Matrix oneda = new Matrix(mfit, 1);
            Matrix temp = new Matrix(mfit, mfit);
            //initialization
            mrqcof(a,alpha,beta);
            for(j = 0 ; j < ma; j++) atry[j]=a[j];
            ochisq=chisq;
            for(iter = 0; iter < itmax; iter++) {
                if(done == ndone) alamda=0.; //last pass. Use zero alamda
                for(j = 0; j < mfit; j++) { 
                    //alter linearized fitting matrix, by augmenting diagonal 
                    //elements
                    for (k = 0; k < mfit; k++){
                        covar.setElementAt(j, k, alpha.getElementAt(j, k));
                    }
                    covar.setElementAt(j, j, alpha.getElementAt(j, j) * (1.0 + alamda));
                    for(k = 0; k < mfit; k++){
                        temp.setElementAt(j, k, covar.getElementAt(j, k));
                    }
                    oneda.setElementAt(j, 0, beta[j]);
                }
                GaussJordanElimination.process(temp, oneda); //matrix solution
                for(j = 0; j < mfit; j++){
                    for(k = 0; k < mfit; k++){
                        covar.setElementAt(j, k, temp.getElementAt(j, k));
                    }
                    da[j] = oneda.getElementAt(j, 0);
                }
                if(done == ndone) { //Converged. Clean up and return
                    covsrt(covar);
                    covsrt(alpha);
                    
                    resultAvailable = true;
                    
                    return;
                }
                for(j = 0, l = 0; l < ma; l++){ //did the trial succeed?
                    if(ia[l]) atry[l] = a[l] + da[j++];
                }
                mrqcof(atry,covar,da);
                if(Math.abs(chisq - ochisq) < Math.max(tol, tol * chisq)) done++;
                if(chisq < ochisq) { //success, accept the new solution
                    alamda *= 0.1;
                    ochisq=chisq;
                    for(j = 0; j < mfit; j++) {
                        for (k = 0; k < mfit; k++){
                            alpha.setElementAt(j, k, covar.getElementAt(j, k));
                        }
                        beta[j] = da[j];
                    }
                    for (l = 0; l < ma; l++) a[l] = atry[l];
                }else{ //failure, increase alamda
                    alamda *= 10.0;
                    chisq=ochisq;
                }
            }

            //too many iterations
            throw new FittingException("too many iterations");
                
        }catch(FittingException e){
            throw e;
        }catch(Throwable t){
            throw new FittingException(t);
        }                
    }
    
    /**
     * Prevents parameter at position i of linear combination of basis functions 
     * to be modified during function fitting
     * @param i position of parameter to be retained
     * @param val value to be set for parameter at position i
     */
    public void hold(int i, double val){
        ia[i] = false;
        a[i] = val;
    }
    
    /**
     * Releases parameter at position i of linear combination of basis functions
     * so it can be modified again if needed
     * @param i position of parameter to be released
     */
    public void free(int i){
        ia[i] = true;
    }
        
    /**
     * Used by fit to evaluate the linearized fitting matrix alpha, and vector 
     * beta to calculate chi square
     * @param a estimated parameters so far
     * @param alpha curvature (i.e. fitting) matrix
     * @param beta array where derivative increments for each parameter are 
     * stored
     * @throws Throwable if evaluation of function fails
     */
    private void mrqcof(double[] a, Matrix alpha, double[] beta) 
            throws Throwable{
	int i,j,k,l,m;
	double ymod,wt,sig2i,dy;
        double[] dyda = new double[ma];
	for(j = 0; j < mfit; j++) {
            for(k = 0; k <= j; k++){
                alpha.setElementAt(j, k, 0.0);
            }
            beta[j] = 0.;
	}
	chisq = 0.;
	for(i = 0; i < ndat; i++) {
            ymod = evaluator.evaluate(i, x[i], a, dyda);
            sig2i = 1.0 / (sig[i] * sig[i]);
            dy = y[i] - ymod;
            for(j = 0, l = 0; l < ma; l++) {
                if(ia[l]){
                    wt = dyda[l] * sig2i;
                    for(k = 0, m = 0; m < l + 1; m++){
                        if(ia[m]){
                            int index = alpha.getIndex(j, k++);
                            alpha.getBuffer()[index] += wt*dyda[m];
                        }                    
                    }
                    beta[j++] += dy*wt;
                }
            }
            chisq += dy*dy*sig2i;
	}
	for (j = 1; j < mfit; j++){
            for(k = 0; k < j; k++){
                alpha.setElementAt(k, j, alpha.getElementAt(j, k));
            }
        }
    }

    /**
     * Expand in storage the covariance matrix covar, so as to take into account
     * parameters that are being held fixed. (For the latter, return zero 
     * covariances)
     * @param covar covariance matrix
     */
    private void covsrt(Matrix covar) {
	int i,j,k;
	for(i = mfit; i < ma; i++){
            for(j = 0; j < i + 1; j++){
                covar.setElementAt(i, j, 0.0);
                covar.setElementAt(j, i, 0.0);
            }
        }
	k = mfit - 1;
	for(j = ma - 1; j >= 0; j--){
            if(ia[j]){
                double[] buffer = covar.getBuffer();
		for(i = 0; i < ma; i++){
                    int pos1 = covar.getIndex(i, k);
                    int pos2 = covar.getIndex(i, j);
                    swap(buffer, buffer, pos1, pos2);
                }
		for(i = 0; i < ma; i++){
                    int pos1 = covar.getIndex(k, i);
                    int pos2 = covar.getIndex(j, i);
                    swap(buffer, buffer, pos1, pos2);
                }
		k--;
            }
	}
    }
    
    /**
     * Swaps values of arrays at provided positions
     * @param array1 1st array
     * @param array2 2nd array
     * @param pos1 1st position
     * @param pos2 2nd position
     */
    private void swap(double[] array1, double[] array2, int pos1, int pos2){
        double value1 = array1[pos1];
        double value2 = array2[pos2];
        array1[pos1] = value2;
        array2[pos2] = value1;
    }    
}
