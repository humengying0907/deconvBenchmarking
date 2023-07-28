# deconvBenchmarking
**deconvBenchmarking**  is a versatile R package that facilitates comprehensive evaluation of deconvolution methods using bulk data simulated under various strategies. Built with a focus on realistic bulk simulation, this powerful tool enables un-biased benchmarking of deconvolution methods.

## Installation
```````
library("devtools");
install_github("humengying0907/deconvBenchmarking")
```````

## Benchmarking framework
<img src="https://github.com/humengying0907/deconvBenchmarking/assets/54827603/363adbec-7a46-4570-8c2f-aaf7a469acbc" width=60% height=60%>

## Key features
**deconvBenchmarking** provides the following key features that allow for easy and efficient deconvolution benchmarking:

-   **Bulk Simulation Collection**: a curated collection of bulk simulation methods with varying underlying heterogeneity levels to closely mirror real biological variability, as well as simulation methods from previously published benchmarking works ([Avila Cobos et al. 2020](https://doi.org/10.1038/s41467-020-19015-1), [Sturm et al. 2019](https://doi.org/10.1093/bioinformatics/btz363) , and another methods published in the [SCDC](https://github.com/meichendong/SCDC) package. Additionally, the package comes with a set of visualization tools to assess the realism of the bulk simulations in terms of biological variance, pairwise similarity and gene correlations
-   **Built-in Deconvolution Methods**: The package incorporates a diverse selection of pre-implemented deconvolution methods, thoughtfully categorized into four distinct groups: reference-free, regression-based, marker-based, and Bayesian methods. Researchers can easily apply and compare these methods on the simulated bulk data.
-   **Built-in reference construction**: The package provides a collection of reference construction methods to create necessary input for regression-based and marker-based methods, which allows for combined evaluation of both deconvoltion methods and reference construction methods.
-   **User-Friendly Framework**: deconvBenchmarking offers an intuitive and user-friendly benchmarking framework for conducting deconvolution benchmarking experiments. It also offers the flexibility to include user-provided references and deconvolution methods for evaluation, allowing researchers to gain valuable insights into methods being tested.
-   **TCGA Expression Cohorts Integration**: The package also offers options to include TCGA expression cohorts as an additional dataset for deconvolution evaluation. Users can directly compare the pre-calculated purity estimates with the results obtained from deconvolution analyses, allowing for assessment of deconvolution methods in the context of real-world data from large-scale cancer studies.

## Reference
Heterogeneous pseudobulk simulation enables realistic benchmarking of cell-type deconvolution methods \
Mengying Hu, Maria Chikina \
bioRxiv 2023.01.05.522919; doi: https://doi.org/10.1101/2023.01.05.522919
