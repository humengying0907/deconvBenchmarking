
# This module contains functions from external resources.

# This function implements the bulk simulation pipeline described in the article Benchmarking of cell type deconvolution pipelines for transcriptomics data. (Nature Communications; https://doi.org/10.1038/s41467-020-19015-1)
# The original function can be found here: https://github.com/favilaco/deconv_benchmark/blob/master/helper_functions.R
Generator <- function(sce, phenoData, Num.mixtures = 1000, pool.size = 100, min.percentage = 1, max.percentage = 99, seed = 24){

  CT = unique(phenoData$cellType)
  ?stopifnot(length(CT) >= 2)

  set.seed(seed)
  require(dplyr)
  require(gtools)

  cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE)
  colnames(cell.distribution) = c("CT","max.n")

  Tissues = list()
  Proportions = list()

  for(y in 1:Num.mixtures){

    #Only allow feasible mixtures based on cell distribution
    while(!exists("P")){

      num.CT.mixture = sample(x = 2:length(CT),1)
      selected.CT = sample(CT, num.CT.mixture, replace = FALSE)

      P = runif(num.CT.mixture, min.percentage, max.percentage)
      P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
      P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)

      missing.CT = CT[!CT %in% selected.CT]
      missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)

      P = rbind.data.frame(P, missing.CT)
      potential.mix = merge(P, cell.distribution)
      potential.mix$size = potential.mix$expected * pool.size

      if( !all(potential.mix$max.n >= potential.mix$size) | sum(P$expected) != 1){
        rm(list="P")
      }

    }

    # Using info in P to build T simultaneously
    chosen_cells <- sapply(which(P$expected != 0), function(x){

      n.cells = P$expected[x] * pool.size
      chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]],
                      n.cells)

      chosen
    }) %>% unlist()


    T <- Matrix::rowSums(sce[,colnames(sce) %in% chosen_cells]) %>% as.data.frame()
    colnames(T) = paste("mix",y,sep="")

    P = P[,c("CT","expected")]
    P$mix = paste("mix",y,sep="")

    Tissues[[y]] <- T
    Proportions[[y]] <- P

    rm(list=c("T","P","chosen_cells","missing.CT"))

  }

  P = do.call(rbind.data.frame, Proportions)
  T = do.call(cbind.data.frame, Tissues)

  P = data.table::dcast(P, CT ~ mix,
                        value.var = "expected",
                        fun.aggregate = sum) %>% data.frame(.,row.names = 1)

  P = P[,gtools::mixedsort(colnames(P))]

  return(list(T = T, P = P))

}



# This function extract a table of the top-ranked genes from limma model fit, and it is part of DE genes identification pipeline from the article Benchmarking of cell type deconvolution pipelines for transcriptomics data. (Nature Communications; https://doi.org/10.1038/s41467-020-19015-1)
# The original function can be found here: https://github.com/favilaco/deconv_benchmark/blob/master/helper_functions.R
marker.fc <- function(fit2, cont.matrix, log2.threshold = 1, output_name = "markers"){

  topTable_RESULTS = limma::topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2.threshold)
  AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
  topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]

  # since limma::topTable applies make.names() automatically to column names, reverse the names back to original format
  colnames(topTable_RESULTS) = colnames(cont.matrix)

  if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }

  markers <- apply(topTable_RESULTS,1,function(x){
    temp = sort(x)
    ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)

  })

  topTable_RESULTS = topTable_RESULTS[markers,]

  markers <- cbind.data.frame(rownames(topTable_RESULTS),
                              t(apply(topTable_RESULTS, 1, function(x){
                                temp = max(x)
                                if(temp < log2.threshold){
                                  temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
                                } else {
                                  temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
                                }
                                temp
                              })))

  colnames(markers) <- c("gene","log2FC","CT")
  markers$log2FC = as.numeric(as.character(markers$log2FC))
  markers <- markers %>% dplyr::arrange(CT,desc(log2FC))

  markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
  markers$gene <- as.character(markers$gene)
  markers$CT <- as.character(markers$CT)

  return(markers)

}

# The following three functions implement specificity score calculations as described in BayesPrism package
# The original functions can be found at: https://github.com/Danko-Lab/BayesPrism/tree/main/BayesPrism/R
compute.specificity <- function(input.matrix,
                                pseudo.min = 1E-8){

  ref.ct <- norm.to.one (ref = input.matrix, pseudo.min = pseudo.min)

  #compute gene expression specificity score
  exp.spec <- t(ref.ct) /  colSums(ref.ct)
  max.spec <- apply(exp.spec,1,max)

  return(max.spec)
}

norm.to.one <- function(ref,
                        pseudo.min){

  G <- ncol(ref)
  phi <- ref/rowSums(ref) * (1-pseudo.min*G) + pseudo.min

  #if the minimum value is greater than zero. simply normalize by total depth
  min.value <- apply(ref,1,min)
  which.row <- min.value>0
  if(any(which.row)){
    #cat("One or more cell types have all genes with non-zero expression. pseudo.min is not applied to these cell types. \n")
    phi[which.row,] <- ref[which.row,,drop=F]/rowSums(ref[which.row,,drop=F])
  }
  return(phi)
}

collapse <- function(ref, labels){

  stopifnot(nrow(ref) == length(labels))

  #remove NA in labels
  non.na.idx <- !is.na(labels)
  if(sum(!non.na.idx)>0) print("Warning: NA found in the cell type/state labels. These cells will be excluded!")
  labels <- labels[non.na.idx]
  ref <- ref[non.na.idx,]

  labels.uniq <- unique(labels)

  ref.collapsed <- do.call(rbind,
                           lapply(labels.uniq,
                                  function(label.i)
                                    colSums(ref[labels==label.i,,drop=F])
                           )
  )

  rownames(ref.collapsed) <- labels.uniq

  return(ref.collapsed)
}
