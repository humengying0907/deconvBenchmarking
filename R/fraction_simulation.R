##################### fraction simulation ##################
#' Simulate realistic cell-type fractions from beta distribution
#' @description Simulate cell-type fractions by utilizing the Cell-Type Proportion Distribution derived from scRNA.
#'    This approach is recommended when the scMeta contains a sufficient number of distinct samples to accurately capture the distribution
#'    and generate realistic proportions (n >= 10 is recommended).
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param n number of simulated bulk samples
#' @param colnames_of_sample column name that corresponds to cellType in scMeta
#' @param colnames_of_cellType column name that corresponds to sampleID in scMeta
#' @param fixed_cell_type a character denotes the target cell type for which we strive to faithfully preserve its distribution.
#'    It is recommended to set this parameter to the name of the malignant cell types. If left undefined, the function will automatically
#'    select the most abundant cell type as 'fixed_cell_type'.
#' @param min.frac minimum fraction in the simulated fraction, values below this threshold will be set to zero. Default = 0.01
#' @param showFractionPlot a logical variable determining whether to display simulated fraction distribution for the fixed_cell_type
#'
#' @return a matrix of simulated fraction, with samples in rows and cell_types in columns
#' @export
#'
#' @examples
#' \dontrun{
#' fracSimulator_Beta(scMeta, n = 100, colnames_of_sample = 'sampleID',
#'                    colnames_of_cellType = 'cell_type', fixed_cell_type = 'malignant')
#' }
fracSimulator_Beta<-function(scMeta,n,
                             colnames_of_sample = NA,
                             colnames_of_cellType = NA,
                             fixed_cell_type = NA,
                             min.frac = 0.01,
                             showFractionPlot = T){

  if(is.na(colnames_of_sample)|is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to sampleID and cell_type in scMeta')
  }

  ct_table=scMeta %>% group_by_(colnames_of_sample, colnames_of_cellType) %>%
    summarise(n=n()) %>%
    mutate(freq=n/sum(n))

  if(is.na(fixed_cell_type)){
    fixed_cell_type=names(table(scMeta$cell_type))[which.max(table(scMeta$cell_type))]
  }

  # for fixed cell type, simulate cell fraction based on beta distribution
  x=ct_table[,'freq'][ct_table[,colnames_of_cellType]==fixed_cell_type]
  fit_x <- fitdistrplus::fitdist(x, "beta")
  fixed_frac=rbeta(n,fit_x$estimate[1],fit_x$estimate[2])

  quick_hist=function(){
    par(mfrow=c(1,2))
    hist(x,main = 'original frac for fixed cell type')
    hist(fixed_frac,main = 'simulated frac for fixed cell type')
  }

  if(showFractionPlot==T){
    quick_hist()
    par(mfrow=c(1,1))
  }

  # for other cell types, simulate cell fraction based on beta distribution, and scaled to the remaining frac
  quick_rbeta<-function(cell_type){
    x=ct_table[,'freq'][ct_table[,colnames_of_cellType]==cell_type]
    if(1 %in% x){
      x[x==1]=x[x==1]-0.01
    }
    fit_x <- fitdistrplus::fitdist(x, "beta")
    rbeta(n,fit_x$estimate[1],fit_x$estimate[2])
  }

  all_cell_types=unique(ct_table[,colnames_of_cellType])
  other_cell_types=all_cell_types[all_cell_types!=fixed_cell_type]
  flexible_frac=do.call(cbind,lapply(other_cell_types,quick_rbeta))

  remaining_frac=1-fixed_frac
  scaling_factor=remaining_frac/rowSums(flexible_frac)
  flexible_frac=sweep(flexible_frac,1,scaling_factor,'*')

  simulated_frac=cbind(fixed_frac,flexible_frac)
  colnames(simulated_frac)=c(fixed_cell_type,other_cell_types)

  simulated_frac[simulated_frac<min.frac] = 0
  simulated_frac = sweep(simulated_frac,1,rowSums(simulated_frac),'/')

  return(simulated_frac)
}


#' Simulate cell-type fractions from dirichlet distribution
#'
#' @param frac_prior A named numeric vector representing the fraction prior. The simulated fraction will fluctuate around the values specified in this vector.
#' @param n number of random fractions to generate
#' @param dispersion_par a numeric value determine the dispersion level of the simulated fractions. With lower value indicating higher dispersion level.
#' @param min.frac minimum fraction in the simulated fraction, values below this threshold will be set to zero
#' @param showFractionPlot a logical variable determining whether to display the boxplot of simulated fractions. It is recommended to set this parameter to TRUE to
#'    visualize the effect dispersion_par and choose the optimal value of it
#'
#' @return a matrix of simulated fraction, with samples in rows and cell_types in columns
#' @export
#'
#' @examples
#' frac_prior = c(100,200,300,400)
#' names(frac_prior) = paste0('cellType',seq(1,4))
#' fracSimulator_Dirichlet(frac_prior, n = 10, dispersion_par = 0.01)

fracSimulator_Dirichlet<-function(frac_prior,n, dispersion_par = 0.1, min.frac = 0, showFractionPlot=T){

  simulated_frac = gtools::rdirichlet(n,frac_prior*dispersion_par)
  colnames(simulated_frac) = names(frac_prior)

  simulated_frac[simulated_frac<min.frac] = 0
  simulated_frac = sweep(simulated_frac,1,rowSums(simulated_frac),'/')

  if(showFractionPlot){
    boxplot(simulated_frac)
  }

  return(simulated_frac)
}


#' Simulate cell-type fractions using Generator() function from Favilaco et al.
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param nbulk number of simulated bulk samples
#' @param pool.size number of cells to aggregate for each simulated bulk sample
#' @param min.percentage minimum percentage of cellType fraction to generate in fraction simulation. Default = 1
#' @param max.percentage maximum percentage of cellType fraction to generate in fraction simulation. Default = 99
#' @param seed a single value. Default = 24
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
fracSimulator_favilaco <- function(scExpr, scMeta, colnames_of_cellType = NA, nbulk = 100, pool.size = 100,
                                   min.percentage = 1, max.percentage = 99, seed = 24){

  # this function incorporated the Generator() function from https://github.com/favilaco/deconv_benchmark/blob/master/helper_functions.R
  # the Generator() function expects the 'phenoData' argument to contain two columns named 'cellID' and 'cellType'
  # Furthermore, this function does not require a pre-defined proportion matrix, it will generate and provide the simulated fraction matrix autmatically

  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  phenoData <- scMeta %>% mutate_('cellType'=colnames_of_cellType)
  phenoData$cellID = rownames(phenoData)

  v = Generator(scExpr,phenoData,Num.mixtures = nbulk, pool.size , min.percentage , max.percentage , seed )
  simulated_frac = as.matrix(t(v$P))
  rownames(simulated_frac) = NULL

  return(simulated_frac)
}


#' Simulate cell-type fractions using generateBulk_norep() function from SCDC package
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param disease indicate the health condition of subjects
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.
#' @param nbulk number of simulated bulk samples
#' @param samplewithRep logical, randomly sample single cells with replacement. Default is T
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
fracSimulator_SCDC = function(scExpr,scMeta,colnames_of_cellType = NA, colnames_of_sample = NA, disease = NULL, ct.sub = NULL,
                              nbulk = 10, samplewithRep = T, prop_mat = NULL){
  require(SCDC,quietly = T)
  eset = Biobase::ExpressionSet(assayData = scExpr,phenoData = new("AnnotatedDataFrame", data = scMeta))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.na(colnames_of_sample)){
    stop('please provide column name that corresponds to sampleID in scMeta')
  }

  if(is.null(ct.sub)){
    ct.sub = unique(scMeta[,colnames_of_cellType])
  }
  v = SCDC::generateBulk_norep(eset,ct.varname=colnames_of_cellType,sample = colnames_of_sample,disease, ct.sub, prop_mat, nbulk, samplewithRep)

  frac = v$true_p
  rownames(frac) = NULL

  return(frac)
}

