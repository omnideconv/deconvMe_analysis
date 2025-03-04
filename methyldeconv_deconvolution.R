library(immunedeconv)
library(methyldeconv)
library(tidyverse)
library(data.table)
library(ggpubr)

source('utils.R')

#### Read processed data ####

tpm_matrix <- readRDS('/nfs/data/NetfLID/methyldeconv_paper/processed_data/tpm.rds')


meta <- readRDS('/nfs/data/NetfLID/methyldeconv_paper/processed_data/meta.rds')


facs <- readRDS('/nfs/data/NetfLID/methyldeconv_paper/processed_data/facs.rds')


meth_matrix <- readRDS('/nfs/data/NetfLID/methyldeconv_paper/processed_data/meth.rds')
unmeth_matrix <- readRDS('/nfs/data/NetfLID/methyldeconv_paper/processed_data/unmeth.rds')

meth_matrix_complete <- meth_matrix[complete.cases(cbind(meth_matrix, unmeth_matrix)),]
unmeth_matrix_complete <- unmeth_matrix[complete.cases(cbind(meth_matrix, unmeth_matrix)),]

methyl_set <- minfi::MethylSet(Meth = meth_matrix_complete, 
                               Unmeth = unmeth_matrix_complete, 
                               preprocessMethod = 'PreprocessRaw')
methyl_set$Sex <- meta$Gender



#### Run methyldeconv ####

methylresolver <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'methylresolver', doPar = T, numCores = 8, alpha = 1) |> data.frame(check.names = F)

epidish <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'epidish', reference = 'blood', maxit= 100) |> data.frame(check.names = F)

methylcc <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'methylcc', array = 'EPIC') |> data.frame(check.names = F)

houseman <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'houseman', array = 'EPIC', compositeCellType = 'Blood', referencePlatform = 'IlluminaHumanMethylationEPIC') |> data.frame(check.names = F)

methatlas <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'methatlas', use_epic_reference=T) |> data.frame(check.names = F)


#### Run immunedeconv ####

quantiseq <- immunedeconv::deconvolute(gene_expression = tpm_matrix, method = 'quantiseq', tumor = F, arrays = F, scale_mrna = T)

epic <- immunedeconv::deconvolute(gene_expression = tpm_matrix, method = 'epic', tumor = F, arrays = F, scale_mrna = T)

immunedeconv::set_cibersort_binary('~/tools_index/cibersort/CIBERSORT.R')
immunedeconv::set_cibersort_mat('~/tools_index/cibersort/LM22.txt')
cibersort <- immunedeconv::deconvolute(gene_expression = tpm_matrix, method = 'cibersort', tumor = F, arrays = F, scale_mrna = T)


#### Prepare results for plotting ####

methylresolver_clean <- methylresolver %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'MethylResolver')

epidish_clean <- epidish %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'EpiDISH')

methylcc_clean <- methylcc %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'MethylCC')

houseman_clean <- houseman %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'Houseman')

methatlas_clean <- methatlas %>%
  rownames_to_column(var = 'sample') %>% 
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'MethAtlas')

methyldeconv_df <- bind_rows(methylresolver_clean, epidish_clean, methylcc_clean, houseman_clean, methatlas_clean)

# clean up cell type labels
methyldeconv_df$celltype_clean <- recode(methyldeconv_df$celltype,
                                         .default = 'other',
                                         
                                         "Mon"="Monocytes",
                                         "Dendritic"="Dendritic cells",
                                         "Macro"="Macrophages",
                                         "Neu"="Neutrophils",
                                         "Eos"="Eosinophils",
                                         "Treg"="Tregs",
                                         "Tnaive"="T cells CD4",
                                         "Tmem"="T cells CD4",     
                                         "CD8"="T cells CD8",
                                         "NK"="NK cells",
                                         "Bcell"="B cells",
                                         "B"="B cells",
                                         "CD4T"="T cells CD4",
                                         "CD8T"="T cells CD8",
                                         "Eosino"="Eosinophils",
                                         "Mono"="Monocytes",
                                         "Neutro"="Neutrophils",
                                         "Gran"="Granulocytes",
                                         "B-cells_EPIC" = "B cells",
                                         "CD4T-cells_EPIC" = "T cells CD4",
                                         "CD8T-cells_EPIC" = "T cells CD8",
                                         "Monocytes_EPIC" = "Monocytes",
                                         "Neutrophils_EPIC" = "Neutrophils",
                                         "NK-cells_EPIC" = "NK cells")

methyldeconv_df <- methyldeconv_df |> 
  group_by(sample, celltype_clean) |> 
  dplyr::summarize(mean_value = mean(value), .groups = 'drop') |>
  group_by(sample) |> 
  do(mutate(., value = mean_value / sum(mean_value))) |>
  ungroup() |>
  mutate(method = 'combined') |> 
  bind_rows(methyldeconv_df) 

methyldeconv_df$celltype_rough <- recode(methyldeconv_df$celltype_clean,
                                         .default = 'other',  
                                         
                                         "B cells" = "B cells",
                                         "Monocytes" = "Monocytes",
                                         "NK cells" = "NK cells",
                                         "T cells CD4" = "T cells CD4",
                                         "T cells CD8" = "T cells CD8")

samples <- colnames(quantiseq)[-1]
quantiseq_clean <- quantiseq %>% 
  data.table::transpose(make.names = 1) %>% 
  mutate(sample = samples) %>%
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'quanTIseq')

samples <- colnames(epic)[-1]
epic_clean <- epic %>% 
  data.table::transpose(make.names = 1) %>% 
  mutate(sample = samples) %>%
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'EPIC')

samples <- colnames(cibersort)[-1]
cibersort_clean <- cibersort %>% 
  data.table::transpose(make.names = 1) %>% 
  mutate(sample = samples) %>%
  pivot_longer(cols = -sample, names_to = 'celltype') %>% 
  mutate(method = 'CIBERSORT')

immunedeconv_df <- bind_rows(quantiseq_clean, epic_clean, cibersort_clean)

# clean up cell type labels
immunedeconv_df$celltype_clean <- recode(immunedeconv_df$celltype,
                                         `B cell` = "B cells",
                                         `B cell naive` = "B cells",
                                         `B cell memory` = "B cells",
                                         `B cell plasma` = "B cells",
                                         
                                         `Macrophage M1` = "Macrophages",
                                         `Macrophage M2` = "Macrophages",
                                         `Macrophage M0` = "Macrophages",
                                         `Macrophage` = "Macrophages",
                                         
                                         `Monocyte` = "Monocytes",
                                         
                                         `Neutrophil` = "Neutrophils",
                                         
                                         `NK cell` = "NK cells",
                                         `NK cell resting` = "NK cells",
                                         `NK cell activated` = "NK cells",
                                         
                                         `T cell CD4+ (non-regulatory)` = "T cells CD4",
                                         `T cell CD4+` = "T cells CD4",
                                         `T cell CD4+ naive` = "T cells CD4",
                                         `T cell CD4+ memory resting` = "T cells CD4",
                                         `T cell CD4+ memory activated` = "T cells CD4",
                                         `T cell regulatory (Tregs)` = "Tregs",
                                         `T cell follicular helper` = "T cells CD4",
                                         
                                         `T cell CD8+` = "T cells CD8",
                                         
                                         `Myeloid dendritic cell` = "Dendritic cells",
                                         `Myeloid dendritic cell resting` = "Dendritic cells",
                                         `Myeloid dendritic cell activated` = "Dendritic cells",
                                         
                                         `Endothelial cell` = "Endothelial cell",
                                         
                                         `Eosinophil` = "Eosinophils",
                                         
                                         .default = "other"
)

immunedeconv_df$celltype_rough <- recode(immunedeconv_df$celltype,
                                         `B cell` = "B cells",
                                         `B cell naive` = "B cells",
                                         `B cell memory` = "B cells",
                                         `B cell plasma` = "B cells",
                                         
                                         `Monocyte` = "Monocytes",
                                         
                                         `NK cell` = "NK cells",
                                         `NK cell resting` = "NK cells",
                                         `NK cell activated` = "NK cells",
                                         
                                         `T cell CD4+ (non-regulatory)` = "T cells CD4",
                                         `T cell CD4+` = "T cells CD4",
                                         `T cell CD4+ naive` = "T cells CD4",
                                         `T cell CD4+ memory resting` = "T cells CD4",
                                         `T cell CD4+ memory activated` = "T cells CD4",
                                         `T cell follicular helper` = "T cells CD4",
                                         
                                         `T cell CD8+` = "T cells CD8",
                                         
                                         .default = "other"
)

#### Load FACS data ####

facs_df <- facs %>% 
  dplyr::select(-`OLINK ID`, -`Metabolom ID`, -`Sample_Name_transcriptomics (bulk)`, -`ID_FACS`,-`non single cells`) %>%
  dplyr::rename('sample'=Netflid_ID) %>%
  pivot_longer(cols = c(-sample), names_to = 'celltype') %>% 
  mutate(method = 'FACS', value=value/100) |> 
  subset(sample %in% measured_samples)

facs_df$celltype_clean <- recode(facs_df$celltype,
                                 `CD19+ (B-Zellen)` = "B cells",
                                 `CD14+ (Monozyten)` = "Monocytes",
                                 `MAIT Zellen` = "MAIT cells",
                                 `gd T-Zellen` = "gd T cells",
                                 `CD8+, CD4+` = "other",
                                 `CD8-, CD4-` = "other",
                                 `T- Helferzellen` = "T cells CD4",
                                 `zytotoxische T-Zellen` = "T cells CD8",
                                 `DCs` = "Dendritic cells",
                                 `non DCs` = "other",
                                 `NK bright` = "NK cells bright",
                                 `NK dim` = "NK cells dim",
                                 `non NK` = "other",
                                 `non LD-` = "other",
                                 
                                 .default = "other"
)

facs_df$celltype_rough <- recode(facs_df$celltype,
                                 `CD19+ (B-Zellen)` = "B cells",
                                 `CD14+ (Monozyten)` = "Monocytes",
                                 `T- Helferzellen` = "T cells CD4",
                                 `zytotoxische T-Zellen` = "T cells CD8",
                                 `NK bright` = "NK cells",
                                 `NK dim` = "NK cells",
                                 
                                 .default = "other"
)

facs_df$datatype <- 'FACS'
immunedeconv_df$datatype <- 'gene expression'
methyldeconv_df$datatype <- 'DNAm'
results_df <- bind_rows(immunedeconv_df, methyldeconv_df, facs_df)

results_df$method <- factor(results_df$method, levels = c('CIBERSORT','EPIC','quanTIseq','EpiDISH','Houseman','MethylCC','MethylResolver','MethAtlas','combined','FACS'))


#### Scatterplot ####

plot_df <- results_df |> dplyr::select(celltype_rough, sample, datatype, method, value) |> 
  subset(celltype_rough != 'other') |>
  group_by(sample, celltype_rough, method) |>
  dplyr::summarize(value_sum = sum(value)) |>
  pivot_wider(names_from = method, values_from = value_sum) |>
  pivot_longer(cols = c(CIBERSORT,EPIC,quanTIseq,EpiDISH,Houseman,MethylCC,MethylResolver,MethAtlas,combined), names_to = 'method', values_to = 'estimated fraction') |>
  dplyr::rename('true fraction' = 'FACS')

plot_df$method <- factor(plot_df$method, levels = c('EpiDISH','Houseman','MethAtlas','MethylCC','MethylResolver','combined','CIBERSORT','EPIC','quanTIseq'))

p_scatter <- ggplot(plot_df, aes(x=`true fraction`, y=`estimated fraction`))+
  geom_point(size = .8, aes(color=celltype_rough))+
  facet_grid(celltype_rough ~ method, scales = 'free_y')+
  theme_bw()+
  geom_abline(linetype = 'dashed')+
  stat_cor(size=2.5)+
  theme(panel.grid = element_blank(), 
        legend.position = 'top')+
  scale_x_continuous(breaks = c(0, 0.2, 0.4))+
  scale_color_manual('',values = celltype_rough_palette)+
  guides(color = guide_legend(override.aes = list(size = 3.5)))

ggsave(plot = p_scatter, device = 'pdf',filename = 'plots/final_versions/scatter.pdf', width = 3800, height = 2000, units = 'px')


#### Heatmaps ####

tmp <- results_df |> dplyr::select(celltype_rough, sample, datatype, method, value) |> 
  subset(celltype_rough != 'other') |>
  group_by(sample, celltype_rough, method) |>
  dplyr::summarize(value_sum = sum(value)) |>
  pivot_wider(names_from = method, values_from = value_sum) |>
  pivot_longer(cols = c(CIBERSORT,EPIC,quanTIseq,EpiDISH,Houseman,MethylCC,MethylResolver,MethAtlas,combined), names_to = 'method', values_to = 'estimated fraction') |>
  dplyr::rename('true fraction' = 'FACS') |>
  drop_na( `estimated fraction`) 


celltype_metrics <- tmp |> 
  group_by(celltype_rough, method) |>
  dplyr::summarize(correlation = cor.test(`true fraction`, `estimated fraction`)$estimate, 
                   rmse = compute_rmse(`true fraction`, `estimated fraction`)) 

global_metrics <- tmp |> 
  group_by(method) |>
  dplyr::summarize(correlation = cor.test(`true fraction`, `estimated fraction`)$estimate, 
                   rmse = compute_rmse(`true fraction`, `estimated fraction`)) |>
  dplyr::mutate(celltype_rough = 'global')

metrics_df <- bind_rows(celltype_metrics, global_metrics) |> 
  mutate(datatype = ifelse(method %in% c('quanTIseq','EPIC','CIBERSORT'),'gene expression based','DNAm based'))
metrics_df$method <- factor(metrics_df$method, levels = c('EpiDISH','Houseman','MethAtlas','MethylCC','MethylResolver','combined','CIBERSORT','EPIC','quanTIseq'))
metrics_df$celltype_rough <- factor(metrics_df$celltype_rough, levels = c("global", "T cells CD8", "T cells CD4", "NK cells",    "Monocytes",   "B cells"))

metrics_df$y_position <- as.numeric(factor(metrics_df$celltype_rough))
metrics_df$y_position[metrics_df$celltype_rough == "global"] <- metrics_df$y_position[metrics_df$celltype_rough == "global"] - 0.3

metrics_df$x_position <- as.numeric(factor(metrics_df$method))
metrics_df$x_position[metrics_df$method == "combined"] <- metrics_df$x_position[metrics_df$method == "combined"] + 0.3

p_r <- ggplot(metrics_df, aes(x=x_position, y=y_position, fill=correlation, label=round(correlation, 2)))+
  geom_tile(aes(height = 0.95, width = 0.95))+
  facet_grid(~datatype, scales='free')+
  xlab('')+ylab('')+
  scale_fill_gradient2('Pearsons R with \nFACS ground truth', low='#b80000', high = '#1a6c9a', mid = 'white')+
  scale_y_continuous(
    breaks = unique(metrics_df$y_position),
    labels = unique(metrics_df$celltype_rough)
  ) +
  scale_x_continuous(
    breaks = unique(metrics_df$x_position),
    labels = unique(metrics_df$method)
  ) +
  theme_minimal()+
  theme(strip.text = element_text(size=14, face = 'bold'),
        axis.text = element_text(size=12), axis.ticks = element_line(color='black'))+
  geom_text(color=ifelse(metrics_df$correlation > 0.7, 'white','black'))

ggsave(plot = p_r, device = 'pdf',filename = 'plots/final_versions/heatmap_R.pdf', width = 4800, height = 2000, units = 'px')


p_rmse <- ggplot(metrics_df, aes(x=x_position, y=y_position, fill=rmse, label=round(rmse, 2)))+
  geom_tile(aes(height = 0.95, width = 0.95))+
  facet_grid(~datatype, scales='free')+
  xlab('')+ylab('')+
  scale_fill_gradient('RMSE with \nFACS ground truth', low='#1a6c9a', high = 'white', )+
  theme_minimal()+
  scale_y_continuous(
    breaks = unique(metrics_df$y_position),
    labels = unique(metrics_df$celltype_rough)
  ) +
  scale_x_continuous(
    breaks = unique(metrics_df$x_position),
    labels = unique(metrics_df$method)
  ) +
  theme(strip.text = element_text(size=14, face = 'bold'),
        axis.text = element_text(size=12), axis.ticks = element_line(color='black'))+
  geom_text(color=ifelse(metrics_df$rmse < 0.15, 'white','black'))

ggsave(plot = p_rmse, device = 'pdf',filename = 'plots/final_versions/heatmap_RMSE.pdf', width = 4800, height = 2000, units = 'px')

