CLUSTERING PHYLOGENETIC TREES
=============================

This repository contains code for computing clustering of phylogenetic trees
by interpreting the trees as points in the tropical torus $\mathbb{R}^n/\mathbb{R}\mathbf{1}$
and clustering them with respect to the *tropical asymmetric distance*
([Comăneci, Joswig, 2024](https://doi.org/10.1007/s10107-023-01996-8)).
The code builds upon Andrei Comăneci's implementation of this distance within
[Polymake](https://polymake.org) and [OSCAR](https://www.oscar-system.org).

The files in this repository are structured as follows:

* `clustering.jl` contains the essential function `cluster` for computing the clustering of phylogenetic trees.
* `examples.jl` shows how clustering a synthetic dataset of hundred trees on four taxa is performed.
  This case is sufficiently low dimensional to also be visualized.
* `inexact.jl` contains helper functions for working with data with floating point coefficients,
  in particular for working with numerical artifacts.
* `apicomplexa-clustering.jl` contains code for clustering the apicomplexa data set
  reported on in our paper.