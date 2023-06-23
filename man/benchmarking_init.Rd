% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmarking.R
\name{benchmarking_init}
\alias{benchmarking_init}
\title{Create an object for deconvolution benchmarking}
\usage{
benchmarking_init(
  scExpr,
  scMeta,
  fixed_cell_type = NA,
  min.frac = 0.01,
  showFractionPlot = T,
  training_ratio = 0.5,
  bulkSimulator_methods = NULL,
  colnames_of_cellType = NA,
  colnames_of_cellState = NA,
  colnames_of_sample = NA,
  simulated_frac = NULL,
  heter_cell_type = NA,
  ncells_perSample = 500,
  min_chunkSize = 5,
  use_chunk = "all",
  dirichlet_cs_par = 0.1,
  min.percentage = 1,
  max.percentage = 99,
  nbulk = 100,
  seed = 24,
  use_simulated_frac_as_prop_mat = FALSE,
  disease = NULL,
  ct.sub = NULL,
  samplewithRep = TRUE,
  refMarkers_methods = c("limma", "scran"),
  hv_genes = NULL,
  log2FC = 2,
  log2FC_flexible = 1,
  minimum_n = 15,
  maximum_n = 50,
  max.spec_cutoff_for_DE = 0.3,
  sigMatrixList = NULL,
  refMatrix_methods = c("raw", "limma", "scran"),
  markerList = NULL,
  include_tcga = F,
  tcga_abbreviation = NA,
  purity_method = "CPE",
  create_autogeneS_input = F,
  max.spec_cutoff_for_autogeneS = 0.5,
  autogeneS_input_file_name = NULL,
  display_autogeneS_command = T,
  ngen = 5000,
  seed_autogeneS = 0,
  nfeatures = 400,
  mode = "fixed",
  create_cibersortx_input = F,
  downsample = F,
  downsample_ratio = 0.2,
  cibersortx_input_file_name = NULL,
  n.core = 1
)
}
\arguments{
\item{scExpr}{single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns}

\item{scMeta}{a dataframe that stores annotation info of each cells}

\item{fixed_cell_type}{argument for fracSimulator_Beta() function: a character denotes the target cell type for which we strive to faithfully preserve its distribution.
It is recommended to set this parameter to the name of the malignant cell types. If left undefined, the function will automatically
select the most abundant cell type as 'fixed_cell_type'.}

\item{min.frac}{fracSimulator_Beta() argument: minimum fraction in the simulated fraction, values below this threshold will be set to zero. Default = 0.01}

\item{showFractionPlot}{fracSimulator_Beta() argument: a logical variable determining whether to display simulated fraction distribution for the fixed_cell_type}

\item{training_ratio}{ratio of training cells. Default = 0.5}

\item{bulkSimulator_methods}{bulkSimulator() argument: a character vector indicating which bulk simulation methods to use. Use list_bulkSimulator() to check for available method names.}

\item{colnames_of_cellType}{column name that corresponds to cellType in scMeta}

\item{colnames_of_cellState}{column name that corresponds to cellState in scMeta, where cellState contains sub-clustering information for each cellType.  This is an argument required for 'heter_sampleIDfree' method only.}

\item{colnames_of_sample}{column name that corresponds to sampleID in scMeta. Required for 'semi', 'heter', and 'SCDC' methods.}

\item{simulated_frac}{a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns. This argument if required for 'homo', 'semi', 'heter', 'heter_sampleIDfree', 'immunedeconv' methods}

\item{heter_cell_type}{name of the cell_type to maintain the highest level of heterogeneity. It is recommended to set this parameter to the name of the malignant cell types.
This argument is required for 'semi' and 'heter_sampleIDfree' methods}

\item{ncells_perSample}{number of cells to aggregate for each simulated bulk sample. This is an argument required for 'homo', 'semi', 'favilaco', 'immunedeconv' and 'SCDC' methods}

\item{min_chunkSize}{minimum number of cells to aggregate to construct a given cell-type component in the simulated bulk. This is an argument required for 'semi' and 'heter' methods}

\item{use_chunk}{a character indicating which cells to pool together for the a given cell_type. Default='all' other options include 'random'.
When use_chunk = 'all', use all the cells belonging to the same patient for a given cell type to generate the certain cell type component in the simulated bulk;
when use_chunk = 'random', randomly select 50-100\% of the cells belonging to the same patient for a given cell type. This is an argument required for 'semi' and 'heter' methods}

\item{dirichlet_cs_par}{a numeric value determine the dispersion level of the simulated fractions. With lower value indicating higher dispersion level. Default = 0.1. This is an argument required for 'heter_sampleIDfree' method.}

\item{min.percentage}{minimum percentage of cellType fraction to generate in fraction simulation. Default = 1. This argument is only required for 'favilaco'}

\item{max.percentage}{maximum percentage of cellType fraction to generate in fraction simulation. Default = 99. This argument is only required for 'favilaco'}

\item{nbulk}{number of simulated bulk samples. This argument is required for 'favilaco' and 'SCDC' methods.}

\item{seed}{a seed value for 'favilaco' method. Default = 24}

\item{use_simulated_frac_as_prop_mat}{a logical variable to determine whether to use 'simulated_frac' as 'prop_mat' as input for 'SCDC' method}

\item{disease}{indicate the health condition of subjects. This argument is only required for 'SCDC' method.}

\item{ct.sub}{a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used. This argument is only required for 'SCDC' method.}

\item{samplewithRep}{logical, randomly sample single cells with replacement. Default is T. This argument is only required for 'SCDC' method.}

\item{refMarkers_methods}{a character vector specifying the desired methods for generating cell-type specific markers. Use list_refMarkers() to check for available method names.}

\item{hv_genes}{a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes}

\item{log2FC}{log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.}

\item{log2FC_flexible}{a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
all the genes that pass the threshold will be considered as marker genes.}

\item{minimum_n}{minimum number of marker genes for a cell-type}

\item{maximum_n}{maximum number of marker genes of a cell-type}

\item{max.spec_cutoff_for_DE}{specificity score threshold to select for hv_genes. Default = 0.3}

\item{sigMatrixList}{a list of signature matices to derive markers from. This argument is required for 'sigMatrixList' method}

\item{refMatrix_methods}{a character vector specifying the desired methods for generating signature matrices. Use list_refMarix() to check for available method names.}

\item{markerList}{a list of pre-calculated cell-type marker list, where the first level of this list represents different methods used to derive cell-type specific markers, and the second level
comprises the actual cell-type specific genes identified by each method. This argument is required for 'markerList' method}

\item{include_tcga}{a logicial variable determining whether to include tcga in the output object}

\item{tcga_abbreviation}{a character indicating tcga abbreviation for the tcga cohort to include, for example 'SKCM'}

\item{purity_method}{tumor purity estimation method that is utilized as a means of estimating the malignant proportion within the exported object for TCGA expression data.
Available methods include 'ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC' and 'CPE', Default = 'CPE'}

\item{create_autogeneS_input}{a logical variable determine whether to create input data for autogeneS, which is a python based approach to construct signature matrix}

\item{max.spec_cutoff_for_autogeneS}{specificity score threshold to select for hv_genes. Default = 0.5}

\item{autogeneS_input_file_name}{desired file name to save the processed file}

\item{display_autogeneS_command}{a logical variable indicating whether to display the command lines to run autogeneS. By pasting the generated code into the command line, an external Python file will be executed,
which will return the generated signature matrices in the 'autogeneS_output' folder.}

\item{ngen}{a numeric variable indicating number of generations used in autogeneS}

\item{seed_autogeneS}{autogeneS seed argument. Default = 0}

\item{nfeatures}{autogeneS nfeatures argument. Default = 400}

\item{mode}{autogeneS mode argument. Default = 'fixed'}

\item{create_cibersortx_input}{a logical variable determine whether to create input data for cibersortx, which is a web-server to construct signature matrix}

\item{downsample}{a logical variable indicating whether to downsample the processed file of not. This argument is useful when the processed file exceeds the storage limit for the cibersortx server}

\item{downsample_ratio}{a numeric value indicating downsampling ratio. This argument determines the fraction of the original size that will be retained in the downsampling process}

\item{cibersortx_input_file_name}{desired file name to save the processed file}

\item{n.core}{number of cores to use for parallel programming. Default = 1}
}
\value{
a list containing the following elements: 1) a list of training/testing cells; 2) a list of simulated bulk object and/or tcga expression;
3) a list of cell-type specific markers; 4) a list of signature matrices
}
\description{
This function takes scRNA profile as input and generate an object intended for future deconvolution benchmarking.
The function performs the following steps: (1) It divides the cells into training and testing cells;
(2) the training cells are utilized to generate reference profiles, such as markers and signature matrices;
(3) the testing cells are used to generate simulated bulk expression, which is then employed for deconvolution purposes;
(4) additionally, this function offers the flexibility to include a TCGA cohort as part of the object for future deconvolution benchmarking
}
\examples{
\dontrun{
# a standard benchmarking pipeline
benchmarking_init(scExpr = scExpr,
                  scMeta = scMeta,
                  fixed_cell_type = 'malignant',
                  bulkSimulator_methods = c('homo', 'semi','heter','favilaco','immunedeconv','SCDC'),
                  colnames_of_cellType = 'cell_type',
                  colnames_of_sample = 'sampleID',

                  # argument for semi bulk simulation method
                  heter_cell_type = 'malignant',

                  # argument for marker constructions
                  refMarkers_methods = c('limma','scran'),

                  # arguments for signature matrics construction
                  refMatrix_methods = c('raw','limma','scran'),

                  # arguments to include tcga
                  include_tcga = T,
                  tcga_abbreviation = 'SKCM',
                  purity_method = 'CPE',

                  # export files for autogeneS and cibersortx
                  create_autogeneS_input = T,
                  create_cibersortx_input = T

                  n.core = 4
                  )

# generate a benchmarking object containing only TCGA cohort
# and use all the single cells to generate reference markers and signature matrices
benchmarking_init(scExpr = scExpr,
                  scMeta = scMeta,
                  training_ratio = 1,
                  bulkSimulator_methods = NULL,

                  # argument for marker constructions
                  refMarkers_methods = c('limma','scran'),

                  # arguments for signature matrics construction
                  refMatrix_methods = c('raw','limma','scran'),

                  # arguments to include tcga
                  include_tcga = T,
                  tcga_abbreviation = 'SKCM',
                  purity_method = 'CPE'
                  )
}
}