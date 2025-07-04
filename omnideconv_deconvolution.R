library(omnideconv)
library(tidyverse)
library(data.table)
library(ggpubr)

setwd('/nfs/proj/omnideconv_benchmarking/omnideconv/methyldeconv_paper/')
source('utils.R')

#### Read processed data ####

tpm_matrix <- readRDS('data/tpm.rds')

meta <- readRDS('data/meta.rds')

#### Use Hao-sampled as scRNA-seq reference ####

norm_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/matrix_norm_counts.rds')
counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/matrix_counts.rds')
sc_annotation <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/celltype_annotations.rds')
sc_batch <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/batch.rds')


res_bisque <- omnideconv::deconvolute(bulk_gene_expression = tpm_matrix, 
                               single_cell_object = norm_counts, 
                               cell_type_annotations = sc_annotation, 
                               batch_ids = sc_batch, model = NULL, 
                               method = 'bisque')
saveRDS(res_bisque, 'deconv_results/omnideconv/bisque.rds')

sig_dwls <- omnideconv::build_model(single_cell_object = counts, cell_type_annotations = sc_annotation, method = 'dwls', batch_ids = sc_batch, bulk_gene_expression = tpm_matrix, verbose = T)
colnames(sig_dwls) <- gsub(' ','_',colnames(sig_dwls))

res_dwls <- omnideconv::deconvolute(bulk_gene_expression = tpm_matrix, model = sig_dwls, method = 'dwls', verbose=T)
colnames(res_dwls) <- gsub('_',' ',colnames(res_dwls))
saveRDS(res_dwls, 'deconv_results/omnideconv/dwls.rds')


sig_scaden <- build_model(single_cell_object = counts, cell_type_annotations = sc_annotation, method = 'scaden', batch_ids = sc_batch, bulk_gene_expression = tpm_matrix, verbose = T)
saveRDS(sig_scaden, 'deconv_results/omnideconv/scaden_sig.rds')

res_scaden <- omnideconv::deconvolute(bulk_gene_expression = tpm_matrix, model = sig_scaden, method = 'scaden', verbose=T)
saveRDS(res_scaden, 'deconv_results/omnideconv/scaden.rds')

plot_deconvolution(list('dwls'=res_dwls,'bisque'=res_bisque, 'scaden'=res_scaden), plot_method = 'jitter')
