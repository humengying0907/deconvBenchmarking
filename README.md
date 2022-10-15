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
We provide a function to simulate cell proportions based on the cell fraction distributions of the single-cell profiles. Calling ```simulated_frac()``` will automatically return paired histograms comparing distribution of cell fraction distributions for the selected cell type. 
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

With semi simulation, we restricted that the malignant parts of each synthetic bulk sample come from the same patient, while the non-malignant parts are randomly selected regardless of where they are from. ```chunk_size_threshold_for_fixed_cell``` controls the minimum number of malignant cells required for a synthetic bulk sample.
```r
semi=create_semiheterSimulation(simulated_frac,scExpr,scMeta,
                                colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                                fixed_cell_type='malignant',ncells_perSample = 500,chunk_size_threshold_for_fixed_cell = 10)
```
#### 3. heterogeneouse simulation

We restricted that both malignant and non-malignant parts of each synthetic bulk sample come from the same patient. ```chunk_size_threshold``` controls the minimum number of cells to aggregate for each cellular components.
```r
heter=create_heterSimulation(simulated_frac,scExpr,scMeta,
                             colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                             chunk_size_threshold = 3)
```
## Run benchmarking pipeline
We provide a benchmarking pipeline that evaluates the performance of different deconvolution methods. 
#### 1. create simulation object

Given a single cell profile, we split the data into training (50%) and testing (50%) sets. Training cells are used to build reference matrices and cell-type specific markers, and testing cells are used to generate synthetic bulk samples. The above steps are repeated 10 times for a comprehensive evaluation.
```r
simulation_repeats=list()
for(i in 1:10){
  simulation_repeats[[i]]=create_simulation_fold(scExpr,scMeta,colnames_of_sample = 'sampleID',colnames_of_cellType = 'cell_type',
                                                 simulated_frac = simulated_frac,
                                                 max.spec_cutoff_for_DE=0.3, # parameters for informative gene filtering
                                                 chunk_size_threshold=3,scale_to_million=T, # heter-simulation parameters
                                                 ncells_perSample=500, # homo-simulation parameter
                                                 fixed_cell_type='malignant',chunk_size_threshold_for_fixed_cell=10, #semiheter-simulation parameters
                                                 log2FC=2,minimum_n=15,maximum_n=50,# marker gene list parameters
                                                 create_autogeneS_input=F,autogeneS_input_file_name=NULL, max.spec_cutoff_for_autogeneS=0.5, # create input for autogeneS
                                                 create_cibersortx_input=F,cibersort_input_file_name=NULL,cibersort_downsample=F,cibersort_downsample_scale=0.2)
  message(paste('>>>>>>>>>>>>>>>>>>>>>>>>> finish',i,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))
}
```

#### 2. run deconvolution on the simulation object
```r
deconv_results=run_deconv(simulation_repeats,run_linseed = T,CellTypeNumber = ncol(simulated_frac))
```

#### 3. evaluate the performance of each deconvolution method
```r
# prepare a mapping table and a mapping list to evaluate the performance
cell_map=get_cell_map(deconv_results)
mapping=get_mapping_list(cell_map,simulated_frac)
deconv_performance_results=run_evalu(deconv_results,simulated_frac,cell_map,mapping)
# deconvolution results (pearson correlation or RMSE value) can be summarized in a long table
performance_table=quick_gather(deconv_performance_results,'maxcor') 
head(performance_table)
```
| cell_type | method | value | group | rep 
------- | ------- | ------- | ------- | -------
|Bcell | MarkerBased_gsva_limmaMarkers | 0.79 | homo | 1
|CAF | MarkerBased_gsva_limmaMarkers | 0.81 | homo | 1
|Endo | MarkerBased_gsva_limmaMarkers | 0.85 | homo | 1
|Macro | MarkerBased_gsva_limmaMarkers | 0.74 | homo | 1
|malignant | MarkerBased_gsva_limmaMarkers | 0.80 | homo | 1
|NK | MarkerBased_gsva_limmaMarkers | 0.60 | homo | 1

#### 4. visualize deconvolution results
```r
deconv_evalu_plot(deconv_performance_results[[1]]$homo_evalu$R$linseed)
```
<img src="https://github.com/humengying0907/deconvBenchmarking/blob/main/images/performance_visualization.png" width=70% height=70%>


