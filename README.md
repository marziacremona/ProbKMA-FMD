# ProbKMA-FMD

R code implementing **ProbKMA** (probabilistic k-mean with local alignment) for local clustering of functinal data and functional motif discovery, proposed in the paper [Probabilistic K-mean with local alignment for clustering and motif discovery in functional data](https://arxiv.org/abs/1808.04773), by Marzia A. Cremona and Francesca Chiaromonte. 


#### Code

## `ProbKMA-FMD_functions.r`
R functions for ProbKMA, cluster evaluation and functional motif discovery.
- `probKMA`: probabilistic k-mean with local alignment to find candidate motifs
- `probKMA_plot`: plot the results of `probKMA`
- `find_candidate_motifs`: run multiple times `probKMA` function with different K, c and initializations, with the aim to find a set of candidate motifs
- `filter_candidate_motifs`: filter the candidate motifs on the basis of a threshold on the average silhouette index and a threshold on the size of the curves in the motif
- `cluster_candidate_motifs`: determine a global radius, group candidate motifs based on their distance, and determine a group-specific radius
- `cluster_candidate_motifs_plot`: plot the results of `cluster_candidate_motifs`
- `motifs_search`: find occurrences of the candidate motifs in the curves and sort them according to their frequencies and radius
- `motifs_search_plot`: plot the results of motifs_search


#### Functional motif discovery example

## `Functional motif discovery on simulated data`


#### Probabilistic local clustering examples

## `Berkley growth curves`

## `Italian Covid-19 excess mortality curves`
