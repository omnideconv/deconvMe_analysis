library(UpSetR)
library(ComplexHeatmap)

source('utils.R')

#### load signatures ####

sig_epidish <- EpiDISH::centDHSbloodDMC.m
sig_houseman <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable
sig_methylcc <- methylCC:::.find_dmrs()
sig_methylcc <- sig_methylcc$regions_all
sig_methylresolver <- as.matrix(MethylResolver::MethylSig)
sig_methatlas <- unique(data.table::fread('signatures/reference_atlas.csv'))
cpgs <- sig_methatlas$CpGs
sig_methatlas$CpGs <- NULL
sig_methatlas <- as.matrix(sig_methatlas)
rownames(sig_methatlas) <- cpgs


sig_quantiseq <- data.table::fread(system.file("extdata", "quantiseq", "TIL10_signature.txt", package = "immunedeconv", mustWork = TRUE))
genes <- sig_quantiseq$ID
sig_quantiseq <- as.matrix(sig_quantiseq[,-1])
rownames(sig_quantiseq) <- genes

sig_epic <- EPIC::BRef$refProfiles[EPIC::BRef$sigGenes,]

sig_cibersort <- data.table::fread('signatures/LM22.txt')
genes <- sig_cibersort$`Gene symbol`
sig_cibersort <- as.matrix(sig_cibersort[,-1])
rownames(sig_cibersort) <- genes

# unify cell type names in columns
colnames(sig_epidish) <- c('B cells','NK cells','T cells CD4', 'T cells CD8','Monocytes','Neutrophils','Eosinophils')
colnames(sig_epic) <- c('B cells','T cells CD4', 'T cells CD8','Monocytes','Neutrophils','NK cells')
colnames(sig_quantiseq) <- c('B cells','Macrophages M1','Macrophages M2','Monocytes','Neutrophils','NK cells','T cells CD4','T cells CD8','Tregs','Dendritic cells')
colnames(sig_methylresolver) <- c('Monocytes','Dendritic cells','Macrophages','Neutrophils','Eosinophils','Tregs','T cells CD4','T memory cells','T cells CD8','NK cells','B cells')
colnames(sig_houseman) <- c('T cells CD8','T cells CD4','NK cells','B cells','Monocytes','Neutrophils')
colnames(sig_cibersort) <- c('B cells','B cells memory','Plasma cells','T cells CD8','T cells CD4','T memory cells','T memory cells activated','T cells helper','Tregs','T cells gd','NK cells','NK cells activated','Monocytes','Macrophages','Macrophages M1','Macrophages M2','Dendritic cells','Dendritic cells activated','Mast resting','Mast activated','Eosinophils','Neutrophils')
sig_methylcc$cellType <- recode(sig_methylcc$cellType, 'Gran' = 'Granulocytes','CD4T'='T cells CD4','Bcell'='B cells','Mono'='Monocytes','NK'='NK cells','CD8T'='T cells CD8')
colnames(sig_methatlas) <- recode(colnames(sig_methatlas), 'Monocytes_EPIC' = 'Monocytes','B-cells_EPIC'='B cells','CD4T-cells_EPIC'='T cells CD4','NK-cells_EPIC'='NK cells','CD8T-cells_EPIC'='T cells CD8','Neutrophils_EPIC'='Neutrophils')

signatures <- list('cibersort'=sig_cibersort, 
                   'epic'=sig_epic,
                   'epidish'=sig_epidish,
                   'houseman'=sig_houseman,
                   'methatlas'=sig_methatlas,
                   'methylcc'=sig_methylcc,
                   'methylresolver'=sig_methylresolver,
                   'quantiseq'=sig_quantiseq)


#### Count overlaps between signatures ####

lt <- list('EpiDISH'=get_sig_annotation(sig_epidish, '450k'), 
           'Houseman'=get_sig_annotation(sig_houseman, 'epic'), 
           'MethylResolver'=get_sig_annotation(sig_methylresolver, '450k'),
           'MethylCC'=sig_methylcc,
           'MethAtlas'=get_sig_annotation(sig_methatlas,'epic'),
           'quanTIseq'=get_sig_annotation_genes(sig_quantiseq, include_promoters = T),
           'EPIC'=get_sig_annotation_genes(sig_epic, include_promoters = T),
           'CIBERSORT'=get_sig_annotation_genes(sig_cibersort, include_promoters = T))
mt <- suppressWarnings(make_comb_mat(lt))

p_upset <- UpSet(mt, 
                 top_annotation = upset_top_annotation(mt, add_numbers=T), 
                 right_annotation = upset_right_annotation(mt, add_numbers = T))

pdf('plots/final_versions/upset_basepairs.pdf', width = 1200, height = 400)
p_upset
dev.off()
