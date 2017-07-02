/**
 * @file
 * This file contains unit tests for
 * com.irurueta.numerical.signal.processing.KalmanMeasurementNoiseCovarianceEstimator
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date October 14, 2015
 */
package com.irurueta.numerical.signal.processing;

import com.irurueta.algebra.Matrix;
import com.irurueta.algebra.WrongSizeException;
import com.irurueta.numerical.signal.processing.AccelerationFileLoader.Data;
import java.io.File;
import java.io.IOException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

public class MeasurementNoiseCovarianceEstimatorTest {
    
    public static final String ACCELERATION_NO_MOTION = 
            "./src/test/java/com/irurueta/numerical/signal/processing/acceleration-no-motion.dat";
    
    public static final double ABSOLUTE_ERROR = 1e-6;
    
    public static final double LARGE_ABSOLUTE_ERROR = 1e-3;
    
    public MeasurementNoiseCovarianceEstimatorTest() {}
    
    @BeforeClass
    public static void setUpClass() {}
    
    @AfterClass
    public static void tearDownClass() {}
    
    @Before
    public void setUp() {}
    
    @After
    public void tearDown() {}
    
    @Test
    public void testConstructor() throws SignalProcessingException{
        
        MeasurementNoiseCovarianceEstimator estimator =
                new MeasurementNoiseCovarianceEstimator(3);
        
        //check correctnes
        assertEquals(estimator.getMeasureParams(), 3);
        assertEquals(estimator.getMeasurementNoiseCov().getRows(), 3);
        assertEquals(estimator.getMeasurementNoiseCov().getColumns(), 3);
        assertEquals(estimator.getSampleAverage().length, 3);
        assertEquals(estimator.getSampleCount(), 0);
        
        //Force IllegalArgumentException
        estimator = null;
        try{
            estimator = new MeasurementNoiseCovarianceEstimator(0);
            fail("IllegalArgumentException expected but not thrown");
        }catch(IllegalArgumentException e){}        
        assertNull(estimator);
    }
    
    @Test
    public void testUpdate() throws SignalProcessingException, IOException, 
            WrongSizeException{

        MeasurementNoiseCovarianceEstimator estimator =
                new MeasurementNoiseCovarianceEstimator(3);
        
        File f = new File(ACCELERATION_NO_MOTION);
        Data data = AccelerationFileLoader.load(f);
        
        double[] sample = new double[3];
        
        for(int i = 0; i < data.numSamples; i++){
            sample[0] = data.accelerationX[i];
            sample[1] = data.accelerationY[i];
            sample[2] = data.accelerationZ[i];
            
            estimator.update(sample);
        }
        
        Matrix measurementNoiseCov = estimator.getMeasurementNoiseCov();
        double[] sampleAverage = estimator.getSampleAverage();
        
        assertEquals(data.numSamples, estimator.getSampleCount());
        
        //compute sample average
        double[] sampleAverage2 = new double[3];
        for(int i = 0; i < data.numSamples; i++){
            sampleAverage2[0] += data.accelerationX[i] / data.numSamples;
            sampleAverage2[1] += data.accelerationY[i] / data.numSamples;
            sampleAverage2[2] += data.accelerationZ[i] / data.numSamples;
        }
        
        //check correctness of average
        assertArrayEquals(sampleAverage, sampleAverage2, ABSOLUTE_ERROR);
        
        //compute covariance
        
        //copy all acceleration samples into a matrix (each column contains a 
        //sample)
        Matrix samples = new Matrix(3, data.numSamples);
        for(int i = 0; i < data.numSamples; i++){
            samples.setElementAt(0, i, data.accelerationX[i]);
            samples.setElementAt(1, i, data.accelerationY[i]);
            samples.setElementAt(2, i, data.accelerationZ[i]);
        }
        
        Matrix transSamples = samples.transposeAndReturnNew();
        
        Matrix cov = samples.multiplyAndReturnNew(transSamples);
        cov.multiplyByScalar(1.0 / data.numSamples);
        
        //compare matrices
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                assertEquals(measurementNoiseCov.getElementAt(i, j),
                        cov.getElementAt(i, j), LARGE_ABSOLUTE_ERROR);
            }
        }
    }
}
