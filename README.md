# deconvBenchmarking
Built with a focus on realistic bulk simulation, **deconvBenchmarking** offers a benchmarking framework for comprehensive evaluations of the deconvolution methods.

## Installation
```````
library("devtools");
install_github("humengying0907/deconvBenchmarking")
```````

## Benchmarking framework
The entire benchmarking framework can be operated through the following three core functions:
* `benchmarking_init()`: creates simulated bulk expression and construct references using the training cells
* `benchmarking_deconv()`: performs deconvolution on the benchmarking_obj with user-selected deconvolution methods
* `benchmarking_evalu()`: evaluates the performance of each method in terms of per-cell type correlation and RMSE (root mean square error)
<img src="https://github.com/humengying0907/deconvBenchmarking/assets/54827603/363adbec-7a46-4570-8c2f-aaf7a469acbc" width=60% height=60%>

## Key features
**deconvBenchmarking** provides the following key features that allow for easy and efficient deconvolution benchmarking:

-   **Bulk Simulation Collection**: a curated collection of bulk simulation methods with varying underlying heterogeneity levels to closely mirror real biological variability, as well as simulation methods from previously published benchmarking works ([Avila Cobos et al. 2020](https://doi.org/10.1038/s41467-020-19015-1), [Sturm et al. 2019](https://doi.org/10.1093/bioinformatics/btz363)), and another methods published in the [SCDC](https://github.com/meichendong/SCDC) package. Additionally, the package comes with a set of visualization tools to assess the realism of the bulk simulations in terms of biological variance, pairwise similarity and gene correlations
-   **Built-in Deconvolution Methods**: The package incorporates a diverse selection of pre-implemented deconvolution methods, thoughtfully categorized into four distinct groups: reference-free, regression-based, marker-based, and Bayesian methods. Researchers can easily apply and compare these methods on the simulated bulk data.
-   **Built-in reference construction**: The package provides a collection of reference construction methods to create necessary input for regression-based and marker-based methods, which allows for combined evaluation of both deconvoltion methods and reference construction methods.
-   **User-Friendly Framework**: deconvBenchmarking offers an intuitive and user-friendly benchmarking framework for conducting deconvolution benchmarking experiments. It also offers the flexibility to include user-provided references and deconvolution methods for evaluation, allowing researchers to gain valuable insights into methods being tested.
-   **TCGA Expression Cohorts Integration**: The package also offers options to include TCGA expression cohorts as an additional dataset for deconvolution evaluation. Users can directly compare the pre-calculated purity estimates with the results obtained from deconvolution analyses, allowing for assessment of deconvolution methods in the context of real-world data from large-scale cancer studies.

## Bulk simulation
We provided a useful `bulkSimulator()` function with 7 built-in simulation strategies to select from:
-   **homo**: Homogeneous bulk simulation, single cells are aggregated randomly within each cell type, regardless of where they are from. 
-   **semi**: Semi-heterogeneous bulk simulation, upon single cell aggregation, cells from the `heter_cell_type` are constrained to originate from the same biological sample, where for other cell-types, single cells are aggregated randomly without any consideration of their origin.
-   **heter**: Heterogeneous bulk simulation, each cell-type component in the simulated bulk sample is constrained to originating from the same biological sample.
-   **heter_sampleIDfree**: sampleID-independent heterogeneous bulk simulation, each cell-type component in the simulated bulk sample is constrained to originating from the same sub-cluster. 

And three other methods sourced from previously published benchmarking work or resources that offer bulk simulation capabilities.

-   **favilaco**: [source code](https://github.com/favilaco/deconv_benchmark/blob/master/helper_functions.R#L27C1-L27C1)
-   **immunedeconv**: [source code](https://www.rdocumentation.org/packages/immunedeconv/versions/1.1.0/topics/make_bulk_eset)
-   **SCDC**: [source code](https://rdrr.io/github/meichendong/SCDC/man/generateBulk_norep.html)

## Evaluation of simulated bulk expression
We provided a set of visualization tools to assess the realism of the bulk simulations in terms of biological variance, pairwise similarity and gene correlations. 
The example plots shown below are derived from practical function applications featured in the [tutorial](https://humengying0907.github.io/deconvBenchmarking_tutorial.html).
-  `plot_variance_comparison()`: Comparison of variance: simulated bulk vs real bulk
![p1](https://github.com/humengying0907/deconvBenchmarking/assets/54827603/d83cf638-a53b-41db-85e6-ec0854e75581)

-  `plot_PairwiseSampleCorrelation()`: Comparison of pairwise similarities between simulated samples
<img src="https://github.com/humengying0907/deconvBenchmarking/assets/54827603/4e828e41-0551-4297-9d48-8e6b4b1ec7d0" width=60% height=60%>

-  `plot_PairwiseGeneCorrelation()`: Plot pairwise gene correlation in simulated bulk and real bulk expression
![p3](https://github.com/humengying0907/deconvBenchmarking/assets/54827603/51d1ecbe-a0bd-4a6a-a187-f5b56954f9f1)

## Tutorial
We provided a [tutorial](https://humengying0907.github.io/deconvBenchmarking_tutorial.html) for detailed implementation of the benchmarking framework and bulk simulation using an example dataset.

## Reference
Heterogeneous pseudobulk simulation enables realistic benchmarking of cell-type deconvolution methods \
Mengying Hu, Maria Chikina \
bioRxiv 2023.01.05.522919; doi: https://doi.org/10.1101/2023.01.05.522919
