/**
 * @file
 * This file contains implementation of
 * com.irurueta.numerical.robust.WeightSelection
 * 
 * @author Alberto Irurueta (alberto@irurueta.com)
 * @date November 12, 2016.
 */
package com.irurueta.numerical.robust;

import com.irurueta.sorting.Sorter;
import com.irurueta.sorting.SortingException;
import java.util.Arrays;

/**
 * Class containing the selection that was made on a weighted algorithm.
 * This is used internally by weighted estimators
 */
public class WeightSelection {
    /**
     * Array indicating which correspondences have been selected (i.e. have
     * a true value), and which ones hasn't (have a false value)
     */
    private boolean[] selected;

    /**
     * Number of correspondences that have been selected
     */
    private int numSelected;

    /**
     * Constructor
     */
    private WeightSelection(){}

    /**
     * Returns array indicating which correspondences have been selected
     * (i.e. have a true value), and which ones hasn't (have a false value)
     * @return array indicating which correspondences have been selected.
     */
    public boolean[] getSelected(){
        return selected;
    }

    /**
     * Sets array indicating which correspondences have been selected (i.e.
     * have a true value), and which ones hasn't (have a false value)
     * @param selected array indicating which correspondences have been 
     * selected
     */
    public void setSelected(boolean[] selected){
        this.selected = selected;
    }

    /**
     * Returns number of correspondences that have been selected
     * @return number of correspondences that have been selected
     */
    public int getNumSelected(){
        return numSelected;
    }

    /**
     * Sets number of correspondences that have been selected
     * @param numSelected number of correspondences that have been selected
     */
    public void setNumSelected(int numSelected){
        this.numSelected = numSelected;
    }

    /**
     * Selects correspondences based on provided weights and creates a
     * weight selection instance
     * @param weights weights. The larger its value the more important a 
     * correspondence is
     * @param sortWeights indicates whether weights must be sorted so that 
     * largest weights are taken into account first
     * @param maxPoints maximum number of correspondences to pick
     * @return instance containing the selection that was made
     * @throws SortingException if weights couldn't be sorted
     */
    public static WeightSelection selectWeights(double[] weights, 
            boolean sortWeights, int maxPoints) throws SortingException{

        int length = weights.length;

        //instantiate selected array with all its values as unselected
        boolean[] selected = new boolean[length];
        int numSelected;

        if(sortWeights){
            //sort weights

            //copy weights because this array will be sorted
            double[] weightsCopy = Arrays.copyOf(weights, length);
            Sorter sorter = Sorter.create();
            //array that will contain original indices in ascending order of
            //weights after sorting
            int[] indices = sorter.sortWithIndices(weightsCopy);

            //traverse indices array from greatest position which corresponds
            //to the greatest weight value after sorting in decreasing order
            //up to maxPoints positions
            int counter = 0;
            for(int i = length - 1; i >= 0; i--){
                selected[indices[i]] = true;
                counter++;
                if(counter >= maxPoints) break;
            }
            numSelected = counter;
        }else{
            //weights aren't sorted
            if(length < maxPoints){
                //we select all points
                Arrays.fill(selected, true);
                numSelected = length;                
            }else{
                //weights aren't sorted so we pick the first maxPoints
                for(int i = 0; i < maxPoints; i++){
                    selected[i] = true;
                }
                numSelected = maxPoints;
            }
        }

        WeightSelection selection = new WeightSelection();
        selection.setSelected(selected);
        selection.setNumSelected(numSelected);
        return selection;
    }       
}
