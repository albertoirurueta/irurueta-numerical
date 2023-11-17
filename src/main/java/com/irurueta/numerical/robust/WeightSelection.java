/*
 * Copyright (C) 2016 Alberto Irurueta Carro (alberto@irurueta.com)
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
package com.irurueta.numerical.robust;

import com.irurueta.sorting.Sorter;
import com.irurueta.sorting.SortingException;

import java.util.Arrays;

/**
 * Class containing the selection that was made on a weighted algorithm.
 * This is used internally by weighted estimators.
 */
public class WeightSelection {

    /**
     * Array indicating which correspondences have been selected (i.e. have
     * a true value), and which ones hasn't (have a false value).
     */
    private boolean[] selected;

    /**
     * Number of correspondences that have been selected.
     */
    private int numSelected;

    /**
     * Constructor.
     */
    private WeightSelection() {
    }

    /**
     * Returns array indicating which correspondences have been selected
     * (i.e. have a true value), and which ones hasn't (have a false value).
     *
     * @return array indicating which correspondences have been selected.
     */
    public boolean[] getSelected() {
        return selected;
    }

    /**
     * Sets array indicating which correspondences have been selected (i.e.
     * have a true value), and which ones hasn't (have a false value).
     *
     * @param selected array indicating which correspondences have been
     *                 selected.
     */
    public void setSelected(final boolean[] selected) {
        this.selected = selected;
    }

    /**
     * Returns number of correspondences that have been selected.
     *
     * @return number of correspondences that have been selected.
     */
    public int getNumSelected() {
        return numSelected;
    }

    /**
     * Sets number of correspondences that have been selected.
     *
     * @param numSelected number of correspondences that have been selected.
     */
    public void setNumSelected(final int numSelected) {
        this.numSelected = numSelected;
    }

    /**
     * Selects correspondences based on provided weights and creates a
     * weight selection instance.
     *
     * @param weights     weights. The larger its value the more important a
     *                    correspondence is.
     * @param sortWeights indicates whether weights must be sorted so that
     *                    largest weights are taken into account first.
     * @param maxPoints   maximum number of correspondences to pick
     * @return instance containing the selection that was made.
     * @throws SortingException if weights couldn't be sorted.
     */
    public static WeightSelection selectWeights(
            final double[] weights, final boolean sortWeights, final int maxPoints) throws SortingException {

        final int length = weights.length;

        // instantiate selected array with all its values as unselected
        boolean[] selected = new boolean[length];
        int numSelected;

        if (sortWeights) {
            // sort weights

            // copy weights because this array will be sorted
            final double[] weightsCopy = Arrays.copyOf(weights, length);
            final Sorter<Double> sorter = Sorter.create();
            // array that will contain original indices in ascending order of
            // weights after sorting
            final int[] indices = sorter.sortWithIndices(weightsCopy);

            // traverse indices array from the greatest position which corresponds
            // to the greatest weight value after sorting in decreasing order
            // up to maxPoints positions
            int counter = 0;
            for (int i = length - 1; i >= 0; i--) {
                selected[indices[i]] = true;
                counter++;
                if (counter >= maxPoints) break;
            }
            numSelected = counter;
        } else {
            // weights aren't sorted
            if (length < maxPoints) {
                // we select all points
                Arrays.fill(selected, true);
                numSelected = length;
            } else {
                // weights aren't sorted so we pick the first maxPoints
                for (int i = 0; i < maxPoints; i++) {
                    selected[i] = true;
                }
                numSelected = maxPoints;
            }
        }

        final WeightSelection selection = new WeightSelection();
        selection.setSelected(selected);
        selection.setNumSelected(numSelected);
        return selection;
    }
}
