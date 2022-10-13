# Overview
This directory contains source code to reproduce the results in the manuscript Hu et al., "Heterogeneous bulk simulation enables realistic benchmarking of cell-type".

# Benchmarking framework in this study

<img src="https://github.com/humengying0907/deconvBenchmarking/blob/main/images/framework.png" width=60% height=60%>

# Hetergeneous bulk simulation pipeline in this study

<img src="https://github.com/humengying0907/deconvBenchmarking/blob/main/images/heter_pipeline.png" width=80% height=80%>

Datasets
========
All the datasets used in the manuscript can be downloaded from their respective sources:

R examples
========
## load library dependencies and source code
```r
source('libraries.R')
source('funcs.R')
```
## Input data
We provided an example single-cell profile under data folder that can be used to run our pipeline. It contains an scExpression matrix and a pheno data with cellular annotations. 
```r
# load expression data
scExpr <- as.matrix(fread("data/example_scExpression.csv.gz"),rownames=1)
scExpr[1:5,1:4]
```

| | cell_1 | cell_2 | cell_3 | cell_4
------- | ------- | ------- | ------- | -------
|RPS11 | 445.743085 | 336.21511 | 302.9024 | 394.77306
|ELMO2 | 2.324252 | 50.20218 | 0 | 0
|PNMA1 | 26.842141 | 0 | 0 | 27.67098
|MMP2 | 0 | 0 | 0.075 | 0
|TMEM216 | 0 | 0 | 0 | 52.89860

```r
# load phenotype data
scMeta <- read.delim('data/example_pheno.csv',sep = ',',row.names = 1)
head(scMeta)

```
| | sampleID | cell_type 
------- | ------- | ------- 
|cell_1 | patient1 | malignant 
|cell_2 | patient1 | malignant 
|cell_3 | patient1 | Tcell
|cell_4 | patient1 | malignant
|cell_5 | patient1 | NK 
|cell_6 | patient1 | Tcell 

## Create simulated cell fractions
We provide a function to simulate cell proportions based on the cell fraction distributions of the single-cell profiles. Calling simulated_frac() will automatically return paired histograms comparing distribution of cell fraction distributions for the selected cell type. 
```r
simulated_frac=simulate_frac(scMeta,100,'sampleID','cell_type','malignant')
```
<img src="https://github.com/humengying0907/deconvBenchmarking/blob/main/images/simulated_frac.png" width=50% height=50%>

## Create simulated bulk samples using single cell profiles as input
We provide three bulk simulation strategies that aggregate single cells in different ways.
#### 1. homogeneous simulation
This method aggregates random single cells with pre-defined fractions and adds up the expression values on linear scale.
```r
homo=create_homoSimulation(scExpr,scMeta,colnames_of_cellType = 'cell_type',simulated_frac,ncells_perSample = 500)
```
#### 2. semi-heterogeneous simulation

With semi simulation, we restricted that the malignant parts of each synthetic bulk sample come from the same patient, while the non-malignant parts are randomly selected regardless of where they are from. chunk_size_threshold_for_fixed_cell controls the minimum number of malignant cells required for a synthetic bulk sample.
```r
semi=create_semiheterSimulation(simulated_frac,scExpr,scMeta,
                                colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                                fixed_cell_type='malignant',ncells_perSample = 500,chunk_size_threshold_for_fixed_cell = 10)
```
#### 3. heterogeneouse simulation

We restricted that both malignant and non-malignant parts of each synthetic bulk sample come from the same patient. chunk_size_threshold controls the minimum number of cells to aggregate for each cellular components.
```r
heter=create_heterSimulation(simulated_frac,scExpr,scMeta,
                             colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                             chunk_size_threshold = 3)
```
## Run benchmarking pipeline







