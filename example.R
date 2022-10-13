
source('libraries.R')
source('funcs.R')

scExpr <- as.matrix(fread("data/example_scExpression.csv.gz"),rownames=1)
scMeta <- read.delim('data/example_pheno.csv',sep = ',',row.names = 1)

scExpr[1:5,1:5]
head(scMeta)

simulated_frac=simulate_frac(scMeta,100,'sampleID','cell_type','malignant')

heter=create_heterSimulation(simulated_frac,scExpr,scMeta,
                             colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                             chunk_size_threshold = 3)

semi=create_semiheterSimulation(simulated_frac,scExpr,scMeta,
                                colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
                                fixed_cell_type='malignant',ncells_perSample = 500,chunk_size_threshold_for_fixed_cell = 10)

homo=create_homoSimulation(scExpr,scMeta,
                           colnames_of_cellType = 'cell_type',simulated_frac,ncells_perSample = 500)
