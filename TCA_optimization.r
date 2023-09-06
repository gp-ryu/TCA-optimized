source('~/.scripts/zz_header.r')
library(microbenchmark)
library(TCA)
library(Rcpp)
library(RcppEigen)
library(RcppProgress)
library(RcppClock)
library(matrixcalc)
library(pracma)
library(pbmcapply)
library(futile.logger)
library(nloptr)
library(matrixStats)
library(Matrix)
# library(bigmemory)

source('anRpackage/R/misc.R')
source('anRpackage/R/tca.R')
source('anRpackage/R/tensor.R')
sourceCpp('anRpackage/src/estimate_Z_j.cpp')
sourceCpp('anRpackage/src/matrixMultiplications.cpp')
sourceCpp('anRpackage/src/estimate_Z_j_optimize.cpp')


args <- commandArgs(trailingOnly = T)
cores = 16


# run  ----
## input call----
# args[[1]] <- '../LEVEL3_QCDB/trans/COPD_CHIP_beta_inner.txt'
# args[[2]] <- '../LEVEL3_QCDB/decompose/COPD_CHIP_beta_inner_cellFraction.txt'
#input_X <- fread(args[[1]], nThread = 4, verbose = T)
input_X <- vroom::vroom(args[[1]], delim = ',', altrep = T, num_threads = 2, progress = T ) %>% 
  column_to_rownames('probe') %>% 
  as.matrix

# if (args[[3]])
# tmp_out <- epidish(tmp, ref.m = centDHSbloodDMC.m, method = 'RPC')

epi_W <- vroom::vroom(args[[2]]) %>% 
  column_to_rownames('probe') %>% 
  as.matrix

gc()


## run tca ----
# tca.mdl <- tca_pbmc(input_X, epi_W, constrain_mu = T, debug = T, max_iters = 20)
write_rds(tca.mdl, file = paste0(args[[2]],'.rds'))
gc()
# tca.mdl <- readRDS(paste0(args[[2]],'.rds'))
## run tensor ####
tensor_cpp(input_X, tca.mdl, debug = T)
# tensor_res_cpp <- tensor_cpp(input_X, tca.mdl, debug = T)
# tensor_res_cpp_tmp <- tensor_tmp(input_X, tca.mdl, debug = T)
# names(tensor_res) <- colnames(tca.mdl$W)
## write tensor res ####
# for(cells in colnames(tca.mdl$W)){
#   print (cells)
#   tensor_res[[cells]] %>% as.data.table(keep.rownames = T) %>% setnames('rn','probe') %>% fwrite(., paste0(args[[3]], '_', cells,'.txt'), sep = '\t')
# }













