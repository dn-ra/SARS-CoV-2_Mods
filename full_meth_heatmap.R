library(ComplexHeatmap)
library(viridisLite)
library(circlize)
library(RColorBrewer)

source("read_tombo_wig.R")

wig_files <- list.files(pattern='wig')
null_files <- list.files(pattern='null_samples')
ref_file <- list.files(pattern='ref_coords.txt')
null_wig <- list.files(pattern = "virion")

i <- which(wig_files == null_wig)
if (length(i) > 0) wig_files <- wig_files[-i]

splice <- data.frame(read.csv('orf_alignments.csv', header=T))
rownames(splice) <- splice[,'orf']
splice <- splice[,-1]

meth_matrix <- methpipe(wig_files, null_files, splice, target = 'Fraction')
#meth_matrix_coded <- methpipe(wig_files, null_files, splice, target = 'coded')

#return_diff only needed here to remove transcripts with missing data
base_diff <- lapply(meth_matrix, return_diff, diff = FALSE)
#max_diff <- lapply(meth_matrix_coded, return_diff, max=TRUE)

order_transcripts <- function(matrix_list) {
  tx_order <- c("leader,ORF10", "leader,N" ,    "leader,ORF8" , "leader,ORF7a" ,"leader,ORF6" , "leader,M"  ,   "leader,E", "leader,ORF3a", "leader,S")
  matrix_list <- lapply(matrix_list, function(x) {
    cnam <- colnames(x)[2:(length(x)-1)]
    idx <- order(sort(cnam), cnam)
    return(x[c(1,length(x),idx+1)])
  }
  )
  matrix_list <- matrix_list[tx_order]
  if (any(names(matrix_list) %in% NA)){
    idx <- which(is.na(names(matrix_list)))
    matrix_list <- matrix_list[-idx]
  }
  return(matrix_list)
}

#meth_matrix <- order_transcripts(meth_matrix)

#to apply to basediff list with lapply
full_heatmap <- function(base_diff_element, genome_length = 29893, subtract_base = FALSE) {
  genome_vector <- seq(1,genome_length)
  if (subtract_base == TRUE) {
    baseline_idx <- which(colnames(base_diff_element) == 'baseline_meth') 
    hm_matrix <- matrix(NA,nrow=genome_length, ncol=length(base_diff_element)-2)
    element_vector <- match(base_diff_element$Base, genome_vector)
    #get matrix with meth fractions as difference from baseline (virion)
    hm_matrix[element_vector,] <- as.matrix(apply(base_diff_element[,-c(1, baseline_idx)], 2, FUN = 
                                                  function(x) as.numeric(x) - as.numeric(base_diff_element[,baseline_idx]) ) )
    #attr(hm_matrix, 'exps') <- colnames(base_diff_element)[-c(1, baseline_idx)]
  } else {
  
  hm_matrix <- matrix(NA,nrow=genome_length, ncol=length(base_diff_element)-1)
  element_vector <- match(base_diff_element$Base, genome_vector)
  hm_matrix[element_vector,] <- as.matrix(base_diff_element[,-1])
  #base_diff_element[is.na(base_diff_element)] <- 0
  #levels(base_diff_element$baseline_meth) <- c(levels(base_diff_element$baseline_meth))
  #attr(hm_matrix, 'exps') <- colnames(base_diff_element)[-1]
  }
  return(apply(hm_matrix, MARGIN = 2, FUN = as.numeric)) 
}

base_diff <- order_transcripts(base_diff)
all_heatmap <- lapply(base_diff, full_heatmap, subtract_base = TRUE)

#get rid of transcrripts that don't have all the experiments
expected_experiments = length(wig_files)
idx <- which(lapply(all_heatmap, ncol) == expected_experiments) 
all_heatmap <- all_heatmap[idx]

#experiment names
rem_names <- which(colnames(meth_matrix[[1]]) %in% c("baseline_meth", "Base"))
a <- colnames(meth_matrix[[1]])[-rem_names]

#tail of matrices
b <- do.call(cbind, all_heatmap)

ma <- nrow(all_heatmap[[1]])
mi <- ma-10000
s <- seq.int(mi, ma)
lab <- seq(20000,30000,250)
lab_c <- which(s %in% lab == TRUE)
lab <- s[lab_c]

annot_cols <- as.character(c( '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', 'black', 'grey'))
names(annot_cols) <- names(meth_matrix)
transcript_annot <- rowAnnotation('transcript' = rep(names(all_heatmap), each = expected_experiments),col = list(transcript = annot_cols) )
experiment_annot <- rowAnnotation('experiment' = rep(as.character(a) , times=length(all_heatmap)))
coord_annot <- HeatmapAnnotation('Genome Position' = anno_mark(at = lab_c, labels= lab, which = 'column', side = 'bottom'))


cols <- brewer.pal(n = 11, name = "RdYlBu")
#cols <- viridisLite::viridis(21, option = "A")
col_fun <- colorRamp2(c(-1,0,1), c('blue','antiquewhite','red'))
#col_fun <- colorRamp2(c(-1,0,1), rev(cols[c(2,6,10)]))
b_tail <- t(tail(b, n=10000)) 
b_square <- b_tail^2
b_sign <- sign(b_tail)
  
hm_complete <- Heatmap(b_square * b_sign, na_col = 'white', cluster_rows=FALSE,cluster_columns = FALSE, row_labels = a,  show_column_names = FALSE, show_row_names = FALSE,
row_split = rep(as.character(a) , times=length(all_heatmap)) , col=col_fun, use_raster = F,
row_gap = unit(1, 'mm'), bottom_annotation = coord_annot, width = unit(30, 'cm'), height = unit(20, 'cm'))  +  #right_annotation = transcript_annot,
#name = 'Methylated fraction', row_title_rot = 0,
#row_title_gp = gpar(cex=0.9)) +
transcript_annot + experiment_annot

draw(hm_complete)

pdf('heatmap_mod_denovo_8transcripts_subtract.pdf', width = 15, height = 15)
draw(hm_complete)
dev.off()

writexl::write_xlsx(as.data.frame(b_square * b_sign), path = 'raw_heatmap_data_8txn.xlsx')
