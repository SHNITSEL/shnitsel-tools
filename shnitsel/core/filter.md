# The full flow
First I will do this for a frames object, then work out whether/how it maps onto a tree...
- Choose which features to analyse (chemical intuition at this stage)
- Choose what limits to set. Following plots/statistics may be helpful here:
    - cumulative 2D timeplots showing what the 99th 95th 90th etc. percentiles are at each time
        - for this, I need to work out how to calculate percentiles
            - same as for the confidence intervals? I'm not sure the question is the same... but maybe.
                - but in that case, it will be slow...
    - histograms showing distributions of
        - values over entire ensemble
        - the max of each trajectory
- Then we have three filtration strategies
    1. Complete exclusion of unwanted trajectories
       For this, it would be good to see which combinations of conditions are violated by which trajectories at what time
        - First plot the populations of still-valid trajectories over time for each condition
        - Then determine the most common sequence of defect types for this ensemble, and plot successive intersections 
    2. Individual truncation of trajectories at the first occurance of an unwanted event
        - Histogram of earliest reasons -- how often is each type of unwanted event the earliest unwanted event?
        - Histograms of earliest event times for each reason type
    3. Ensemble level trimming/equalization at a common cutoff, with exclusion of trajectories that don't reach that point.
        - Same plots as for (1.)?
    And apart from making plots, we need to think about how to implement them.