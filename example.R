# load library dependencies and source code
source('libraries.R')
source('funcs.R')

# input data
scExpr <- as.matrix(fread("data/example_scExpression.csv.gz"),rownames=1)
scMeta <- read.delim('data/example_pheno.csv',sep = ',',row.names = 1)
scExpr[1:5,1:5]
head(scMeta)

# Create simulated cell fractions
simulated_frac=simulate_frac(scMeta,100,'sampleID','cell_type','malignant')

# Create simulated bulk samples using single cell profiles as input
heter=create_heterSimulation(simulated_frac,scExpr,scMeta,
                             colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                             chunk_size_threshold = 3)

semi=create_semiheterSimulation(simulated_frac,scExpr,scMeta,
                                colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                                fixed_cell_type='malignant',ncells_perSample = 500,chunk_size_threshold_for_fixed_cell = 10)

homo=create_homoSimulation(scExpr,scMeta,
                           colnames_of_cellType = 'cell_type',simulated_frac,ncells_perSample = 500)

# Run benchmarking pipeline
# 1. create simulation object
simulation_repeats=list()
for(i in 1:10){
  simulation_repeats[[i]]=create_simulation_fold(scExpr,scMeta,colnames_of_sample = 'sampleID',colnames_of_cellType = 'cell_type',
                                                 simulated_frac = simulated_frac,
                                                 max.spec_cutoff_for_DE=0.3, # parameters for hv genes 
                                                 chunk_size_threshold=3,use_chunk = 'all',scale_to_million=T, # heter-simulation parameters
                                                 ncells_perSample=500, # homo-simulation parameter
                                                 fixed_cell_type='malignant',chunk_size_threshold_for_fixed_cell=10, # semiheter-simulation parameters
                                                 log2FC=2,minimum_n=15,maximum_n=50,# marker gene list parameters (maximum_n is set higher in case of downsampling)
                                                 create_autogeneS_input=F,autogeneS_input_file_name=NULL, max.spec_cutoff_for_autogeneS=0.5, # create input for autogeneS
                                                 create_cibersortx_input=F,cibersort_input_file_name=NULL,cibersort_downsample=F,cibersort_downsample_scale=0.2)
  message(paste('>>>>>>>>>>>>>>>>>>>>>>>>> finish',i,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))
}

# 2. run deconvolution on the simulation object
deconv_results=run_deconv(simulation_repeats,run_linseed = F,CellTypeNumber = ncol(simulated_frac))

# 3. evaluate the performance of each deconvolution method
# prepare a mapping table and a mapping list to evaluate the performance
cell_map=get_cell_map(deconv_results)
mapping=get_mapping_list(cell_map,simulated_frac)

deconv_performance_results=run_evalu(deconv_results,simulated_frac,cell_map,mapping)
performance_table=quick_gather(deconv_performance_results,'maxcor')
head(performance_table)

# 4. visualize deconvolution results
deconv_evalu_plot(deconv_performance_results[[1]]$homo_evalu$R$MarkerBased_CAMTHC_scranMarkers)


