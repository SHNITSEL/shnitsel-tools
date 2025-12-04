# Overview
First I will do this for a frames object, then work out whether/how it maps onto a tree...
- Choose which features to analyse (chemical intuition at this stage)
- Choose what limits to set. Following plots/statistics may be helpful here:
    - (1) histograms showing distributions of
        - values over entire ensemble
        - the max of each trajectory
    - (2) cumulative 2D timeplots showing what the 99th 95th 90th etc. percentiles are at each time
        - for this, I need to work out how to calculate percentiles
            - same as for the confidence intervals? I'm not sure the question is the same... but maybe.
                - but in that case, it will be slow...
- Then we have three filtration strategies
    1. Complete exclusion of unwanted trajectories
       For this, it would be good to see which combinations of conditions are violated by which trajectories at what time
        - (3) First plot the populations of still-valid trajectories over time for each condition
        - (4) Then determine the most common sequence of defect types for this ensemble, and plot successive intersections 
        - (8) Bar chart showing how many trajectories would be reincluded by the removal of each criterion alone
        - (9) Heat-maps showing how frequently pairs of violations are solely responsible for exclusion
    2. Individual truncation of trajectories at the first occurance of an unwanted event
        - (6) Histogram of earliest reasons -- how often is each type of unwanted event the earliest unwanted event?
        - (7) Histograms of earliest event times for each reason type
    3. Ensemble level ~~trimming/equalization~~ transection at a common cutoff, with exclusion of trajectories that don't reach that point.
        - (5) Number of frames retained at different transection cutoffs
        - Otherwise, same plots as exclusion
    And apart from making plots, we need to think about how to implement them.

# Workflows
## In general
1. Choose filtranda (manually or by preset) --> filtranda object (unstacked DataArray, or possibly Dataset including cumulative maxima).
2. Use appropriate (strategy-dependent) plots to choose thresholds; these plots require the filtranda object + thresholds (as a coord).
3. Use these criteria (criterion = filtrandum + threshold) to label the data with a `good_upto` coord or similar.
   The rest of the process can be done using the labelled data; the filtranda object can be discarded.
4. At this point, we can analyse the consequences of turning different criteria on and off.
   If we are using the transection strategy, we need to choose a cutoff time at this point.
5. Finally, we apply one of the filtration strategies; the resulting object no longer has labels.
## Omission
1. Choose filtranda
2. `check_thresholds()` or `thesholds_maxima()` Choose thresholds using plot (1): distributions of per-traj maxima
3. `check_criteria()` Consider removing criteria to gain more trajectories using plot (8): bar chart of sole responsability, or (9), heatmap of pairwise responsability
4. `label()` Trajectories are labelled using a coordinate/attribute `keep_trajectory`
    - Or maybe by setting `good_until` to 0 or $end_time as appropriate?
5. `result` The exclusion can be carried out
## Truncation
1. Choose filtranda
2. 
## Transection
1. Choose filtranda
2. `check_thresholds()` or `thresholds_quantiles()` Choose thresholds using plot (2): quantiles of cumulative maxima
3. `valid_populations()` Consider removing criteria to gain more trajectories using plots (3) and (4): Validity populations and intersectional validity/successive intersections
4. Choose a cutoff time using plot (5): number of frames retained at different transection cutoff times
5. Filtration is immediate: we always return the unstacked format, confirming that it has no NaNs

# Objects
filtranda -- this is a DataArray with a `criterion` dimension and coordinate (which need not be an index and may have duplicates)
## states:
- threshold-ready -- if a Dataset contains a 'filtranda' variable
    - for the tree... maybe put this in the root node?
- filter-ready -- if a Dataset contains a `good_upto` or `filtration_mask` variable
    - if both are present, only the `filtration_mask` is used
- normal state -- in the absence of any of the aforementioned variables