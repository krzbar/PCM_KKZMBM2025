# PCM_KKZMBM2025

Code for simulation study and data analysis for the article "Computational singularities in phylogenetic comparative methods" by K. Bartoszek, B. Brahmantio, J. Muñoz-Durán, J. Fuentes-González, J. Pienaar, P. D. Polly.

This repository contains the following `R` scripts:
* `ML_compare.R`
  Contains code to compare the effects of short tip branch on the maximum likelihood inference quality for trees with $n=4$ and $n=100$ tips (1st scenario).
* `ML_compare2.R`
  Contains code to compare the effects of short tip branch on the maximum likelihood inference quality for trees with $n=4$ and $n=100$ tips (2nd scenario).
* `short_tip_remove.R`
  Contains code to compare the choices of traits replacement in a cherry with short tip branches, tested on trees with $n = 5, 25, 50, 100$.
* `DataPrep.R`
  Script to prepare the Caninae trait measurements (RW.csv-first two shape variables after Procrustes superimposition) resulting in mvData.RData; and phylogenetic tree (TreeNonUltra.tre) resulting in Tree.tre for further analysis by the script `data_analysis.R`.
* `data_analysis.R`
  Contains code to compare different methods to resolve the problem of short branch in the canid phylogeny and morphogeometric data.
* `mvSLOUCHversusmvMORPHtimings.R` (in `TIMINGS` subdirectory)
  Contains code to compare running times of different PCM softwares and methods.
