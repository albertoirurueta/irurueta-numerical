/*
 * Copyright (C) 2015 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.fitting;

import com.irurueta.numerical.NotReadyException;

/**
 * Base class for function fitters used to estimate function parameters along
 * with their covariance matrix and chi square value
 */
public abstract class Fitter {
    
    /**
     * Indicates whether result has been estimated and is available for 
     * retrieval
     */
    protected boolean resultAvailable;
    
    /**
     * Returns boolean indicating whether result has been estimated and is 
     * available for retrieval
     * @return true if result has been estimated and is available for retrieval
     */
    public boolean isResultAvailable(){
        return resultAvailable;
    }
    
    /**
     * Indicates whether this instance is ready because enough input data has 
     * been provided to start the fitting process
     * @return true if this fitter is ready, false otherwise
     */
    public abstract boolean isReady();
    
    /**
     * Fits a function to provided data so that parameters associated to that
     * function can be estimated along with their covariance matrix and chi
     * square value
     * @throws FittingException if fitting fails
     * @throws NotReadyException if enough input data has not yet been provided
     */
    public abstract void fit() throws FittingException, NotReadyException;
}
