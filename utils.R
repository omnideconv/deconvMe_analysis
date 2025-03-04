celltype_rough_palette <- c(
  "Monocytes" = "#AA0DB4",
  "other" = "#000000",
  "T cells CD4" = "#2271B2",
  "T cells CD8" = "#3DB7E9",
  "NK cells" = "#e69f00",
  "B cells" = "#359B73"
)

compute_rmse <- function(x, y, zscored = FALSE, weighted = FALSE){
  if(zscored == TRUE){
    sd_x = sd(x)
    sd_y = sd(y)
    if(sd_x == 0){sd_x = 1}
    if(sd_y == 0){sd_y = 1}
    x <- (x - mean(x))/(sd_x)
    y <- (y - mean(y))/(sd_y)
  }
  
  rmse <- sqrt(mean((x - y)^2)) 
  if(weighted == TRUE){rmse <- rmse/max(y)}
  return(rmse)
}

get_sig_annotation <- function(signature, array_type, celltype=NA, liftover=TRUE){
  
  if(array_type == '450k'){
    anno <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations)[rownames(signature),]
  }else{
    anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations)[rownames(signature),]
  }
  
  if(any(is.na(anno$chr))){
    message('Not every CpG could be found in the annotation object. Removing CpGs with no match in annotation.')
    signature <- signature[!is.na(anno$chr),]
    anno <- anno[!is.na(anno$chr),]
  }
  
  if(is.na(celltype)){
    cpg_df <- cbind(anno, signature)
  }else{
    cpg_df <- cbind(anno, 'methylation'=signature[,which(colnames(signature) == celltype)])
  }
  
  cpg_df$start <- as.character(cpg_df$pos)
  cpg_df$end <- as.character(cpg_df$pos)
  dnam_gr <- GenomicRanges::makeGRangesFromDataFrame(df = cpg_df, 
                                                     keep.extra.columns = T, 
                                                     start.field = 'start',
                                                     end.field = 'end', 
                                                     seqnames.field = 'chr', 
                                                     strand.field = 'strand')
  
  chain <- rtracklayer::import.chain("/nfs/data/references/liftover/hg19ToHg38.over.chain") # file downloaded from UCSC
  dnam_gr_hg38 <- unlist(rtracklayer::liftOver(dnam_gr, chain))
  dnam_gr_hg38$pos_hg38 <- ranges(dnam_gr_hg38)@start
  
  
  return(dnam_gr_hg38)
  
}

get_sig_annotation_genes <- function(signature, include_promoters = FALSE){
  
  if(include_promoters){
    genome_annotation_genes <- annotatr::build_annotations(genome = 'hg38', 
                                                           annotations=c('hg38_genes_exons','hg38_genes_introns'))
  }else{
    genome_annotation_genes <- annotatr::build_annotations(genome = 'hg38', 
                                                           annotations=c('hg38_genes_exons','hg38_genes_introns','hg38_genes_promoters'))
  }
  
  sig_genes <- rownames(signature)
  sig_genes <- remap_genes(sig_genes)
  
  sig_gene_coordinates <- parallel::mclapply(sig_genes, function(g){
    anno_g <- genome_annotation_genes |> subset(symbol == g)
    
    if(length(anno_g) == 0 | !any(unique(seqnames(anno_g)) %in% c(paste0('chr', 1:22),'chrX','chrY'))){
      warning(paste0('No annotation found for ', g))
      #return(data.frame(chr=NA, start=NA, end=NA, strand=NA, gene=g))
      return(GRanges())
    }
    return(anno_g)
    
  }, mc.cores=10)
  
  
  
  return(do.call(c, GenomicRanges::reduce(GRangesList(sig_gene_coordinates))))
  
}

# change gene symbol manually in case it does not exist in genome annotation
remap_genes <- function(genes){
  genes <- recode(genes, 
                  "CECR1" = "ADA2",
                  "HIST1H2BC" = "H2BC4",
                  "HIST1H3D" = "H3C4",
                  "FAM46C" = "TENT5C",
                  "FAM160B1" = "FHIP2A",
                  "FAIM3" = "FCMR",
                  'FAM65B' = 'RIPOR2',
                  'KIAA0226L' = 'RUBCNL',
                  'LRMP' = 'IRAG2',
                  'ATHL1' = 'PGGHG',
                  'EMR1' = 'ADGRE1',
                  'EMR2' = 'ADGRE2',
                  'EMR3' = 'ADGRE3',
                  'GPR97' = 'ADGRG3',
                  'VNN3' = 'VNN3P')
  
  return(genes)
}

