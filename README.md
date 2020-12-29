# ProbKMA-FMD

R code implementing **ProbKMA** (probabilistic k-mean with local alignment) for local clustering of functinal data and functional motif discovery, proposed in the paper [Probabilistic K-mean with local alignment for clustering and motif discovery in functional data](https://arxiv.org/abs/1808.04773), by Marzia A. Cremona and Francesca Chiaromonte. 


## Code

#### `ProbKMA-FMD_functions.r`
R functions for ProbKMA, cluster evaluation and functional motif discovery.
- `probKMA`: probabilistic k-mean with local alignment to find candidate motifs
- `probKMA_plot`: plot the results of `probKMA`
- `find_candidate_motifs`: run multiple times `probKMA` function with different K, c and initializations, with the aim to find a set of candidate motifs
- `filter_candidate_motifs`: filter the candidate motifs on the basis of a threshold on the average silhouette index and a threshold on the size of the curves in the motif
- `cluster_candidate_motifs`: determine a global radius, group candidate motifs based on their distance, and determine a group-specific radius
- `cluster_candidate_motifs_plot`: plot the results of `cluster_candidate_motifs`
- `motifs_search`: find occurrences of the candidate motifs in the curves and sort them according to their frequencies and radius
- `motifs_search_plot`: plot the results of motifs_search


## Functional motif discovery example
Functional motif discovery on simulated data: 20 curves embedding 2 functional motifs of length 60, each with 12 occurrences. 
- `len200_sd0.1.RData`: simulated curves
- `len200_sd0.1_simulated_curves_with_motifs.pdf`: plot of curves with true motifs
- `FMD_simulated_data.r`: script to run the example
- `results`: functional motif discovery results

## Probabilistic local clustering examples

### Berkley growth curves
Probabilitstic local clustering of the Berkley Growth Study dataset, provided within the R package `fda` and consisting of the heights of 39 boys and 54 girls recorded from age 1 to 18.
- `growth_smoothed.RData`: smoothed growth curves
- `probKMA_growth.r`: script to run the example

### Italian Covid-19 excess mortality curves
Probabilitic local clustering of Covid-19 excess mortality rate curves (daily difference between 2020 deaths and average deaths in the period 2015-2019) in the 20 regions of Italy. These curves were estimated in the period from February 16, 2020 and April 30, 2020 using the mortality data (due to all causes) from the Italian Institute of Statistics ISTAT. Raw mortality data are available on [ISTAT website](https://www.istat.it/it/files/2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record-4giugno.zip).
- `istat_mortality_rates_smoothed.Rdata`: smoothed excess mortality rate curves
- `probKMA_mortality_regions.r`: script to run the example
