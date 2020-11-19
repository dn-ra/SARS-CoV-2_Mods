#Author: Daniel Rawlinson
#add this into R with 'source(read_tombo_wig.R)' then run the function with 'read_wig(your file)'
#TO plot for a specific transcript, run 'plot_meth(name of data, your transcript of interest)'

#setwd('C:\\Users\\DRAWLINSON\\OneDrive - The University of Melbourne\\Documents\\Manuscripts\\early_evolution_and_genomics')

library(ggplot2)

#read in wig file from tombo. each file is from one experiment, and each experiment can contain multiple transcripts
read_wig <- function(wig_file) {
  wig_data <- list()
  meth_list <- list()
  #open for counting lines
  open_wig <- file(wig_file, 'r')
  i=0
  f_size <- length(readLines(open_wig))
  close(open_wig)
  pb <- txtProgressBar(min=i, max = f_size, initial = 0, title = "Reading methylation Data...")
  #reopen for actual reading
  open_wig <- file(wig_file, 'r')
  while ( TRUE ) {
    line = readLines(open_wig, n=1)
    if ( length(line) == 0 ) {
      #save last list
      current <- as.data.frame(do.call(rbind, wig_data))
      colnames(current) <- c('Base','Fraction')
      current$Base <- as.integer(as.character(current$Base))
      current$Fraction <- as.numeric(as.character(current$Fraction))
      meth_list[[unlist(strsplit(transcript,';')) [2]]] <- current	
      message('Finished!')
      break
    }
    else if ( startsWith(line, 'variableStep') ) {
      if (exists('transcript')) {
        current <- as.data.frame(do.call(rbind, wig_data))
        colnames(current) <- c('Base','Fraction')
        current$Base <- as.integer(as.character(current$Base))
        current$Fraction <- as.numeric(as.character(current$Fraction))
        meth_list[[unlist(strsplit(transcript,';')) [2]]] <- current
        wig_data <- list() #remove all read data
      }
      transcript <- gsub('chrom=', '',strsplit(line, ' ')[[1]][2])
      getTxtProgressBar(pb)
    }
    else if (startsWith(line, 'track')) {
      next
    }
    else { 
      i=i+1
      wig_data[[i]] <- c(strsplit(line, ' ')[[1]])
      setTxtProgressBar(pb, i)
      
      
    }
  }
  
  close(open_wig)
  close(pb)
  message('Imported methylation data for ', length(meth_list), ' transcripts')
  return(meth_list)
}


read_null <- function(null_file) {
  null_data <- read.csv(null_file, sep='\t', header=F)
  name <- gsub('_null_samples.txt','',null_file)
  null_data$experiment <- name
  return(null_data)
}
#add missing bases, set to 0 methylation. applies to object: transcrript within experiment
fill_missing <- function(meth_list_element) {
  all <- seq.int(1:max(meth_list_element$Base))
  #all <- seq.int(1:29893)
  w_missing <- meth_list_element$Fraction[match(all,meth_list_element$Base)]
  w_missing[is.na(w_missing)] <- 0
  missing_added <- as.data.frame(cbind(all, w_missing))
  colnames(missing_added) <- colnames(meth_list_element)
  return(missing_added)
}
#set as 0,1,2 for degree of methyation at each base. applies to object: transcript within experiment
meth_encode <- function(meth_list_element, low_for_0 = 0.5, high_for_2 = 0.8) {
  coded_meth <- cut(meth_list_element$Fraction, c(0, low_for_0, high_for_2, 1.1), labels = c(0,1,2), right=FALSE)
  return(coded_meth)
}
#add experiment names to list. applies to object: all meth data in a list
name_by_experiment <- function(metameth_list, filenames) {
  name_split <- strsplit(filenames, '\\.')
  exps <- sapply(name_split, function(x) x[1])
  names(metameth_list) <- exps
  return(metameth_list)
}

#not in use. Refs wrong.
read_refcoord <- function(ref_file) {
  refs <- read.csv(ref_file, header=F, sep=" ")
  names <- sapply(refs$V1, function(x) strsplit(as.character(x), split = ";"))
  names <- sapply(names, '[', 2)
  coords <- sapply(refs$V3, function(x) strsplit(as.character(x), split = ","))
  coords <- t(sapply(coords, '[', c(1,2,3)))
  diff <- as.numeric(coords[,1])-1
  before_splice <- cbind(as.numeric(coords[,2]) - diff, diff)
  splice_refs <- cbind(before_splice, as.numeric(coords[,3]) - before_splice[,1] -1)
  rownames(splice_refs) <- names
  colnames(splice_refs) <- c('splice_loc', 'before_add', 'after_add')
  return(splice_refs)
}
#add in missing bases and codify. applies to object: experiment
adjust_experiment <- function(meth_list) {
  meth_list  <- lapply(meth_list, fill_missing)
  coded_meth <- lapply(meth_list, meth_encode)
  for (i in 1:length(meth_list)) {
    meth_list[[i]]$coded <- coded_meth[[i]]
  }
  return(meth_list)
}
#where 'target' is either 'coded' or 'Fraction', depending on if you want to extract raw methylation or classified methylation.
get_coded_columns <- function(metameth_list, idxs, transcript, target) {
  transcript_dfs <- lapply(metameth_list[idxs], `[`, transcript)
  codes <- sapply(unlist(transcript_dfs, recursive = FALSE), `[[`, target)
  lens <- sapply(codes, length)
  #check if lengths are different between experiments
  if (!(all(lens == mean(lens)))) codes <- equalize_code_lengths(codes)
  return(as.data.frame(codes))
}
equalize_code_lengths <- function(codes) {
  m <- max(sapply(codes, length))
  for (i in 1:length(codes)) {
    if (length(codes[[i]]) < m) {
      codes[[i]] <- as.character(c(as.vector(codes[[i]]), rep(0, m-length(codes[[i]]))))
    }
  }
  return(codes)
}
#Change the base for each transcript to reference co-ordinates. Applies to object: matrix of coded meth data, output of read_refcooord
convert_to_ref_coords <- function(meth_matrix, ref_coord_data) {
  for (i in names(meth_matrix)) {
    s <- ref_coord_data[i,]
    b <- meth_matrix[[i]]$Base
    #meth_matrix[[i]]$Base <- ifelse(b >= s[1], b + s[2], NA)
    meth_matrix[[i]]$Base <- ifelse(b >= as.integer(s[1]), b + as.integer(s[3]), b)
  }
  return(meth_matrix)
}
#target is either Fraction or coded. Extracts that column to apply to meth_matrix
apply_baseline_meth <- function(meth_matrix_element, base_meth_data, target ) {
  to_take <- meth_matrix_element$Base
  idx <- match(to_take, base_meth_data[[1]]$Base)
  codes <- unlist(sapply(base_meth_data, '[', target), use.names = FALSE)[idx]
  #codes <- base_meth_data[[1]]$coded[idx]
  meth_matrix_element$baseline_meth <- codes
  return(meth_matrix_element)
}
#merge coded methylation data. applies to object: all meth data in list
matrixify <- function(metameth_list, target) {
  meth_matrix <- list()
  tn <- lapply(metameth_list, names)
  tn_u <- unique(unlist(tn))
  tn_link <- lapply(tn_u, grep, tn)
  names(tn_link) <- tn_u
  for (t in tn_u) {
    to_take <- tn_link[[t]]
    codes <- get_coded_columns(metameth_list, to_take, t, target)
    names(codes) <- names(tn)[to_take]
    meth_matrix[[t]] <- codes
    Base <- 1:nrow(codes)
    meth_matrix[[t]] <- cbind(Base, codes)
  }
  return(meth_matrix)
}
#barcode plot for a specific transcript. applies to object: experiment + transcript name
plot_meth <- function(meth_df, transcript_name) {
  meth_50 <- subset(meth_df[[transcript_name]], Fraction > 0.50)
  meth_90 <- subset(meth_df[[transcript_name]], Fraction > 0.90)
  ggplot() + geom_vline(xintercept = meth_50$Base, alpha = meth_50$Fraction, color = '#94bdc9') +
    geom_vline(xintercept = meth_90$Base, color = '#2e66a9') +
    labs(x = 'Genome Position (bp)', subtitle = transcript_name) +
    theme_classic()
}
#get differential methyation sites from coded_meth_matrix. applies to object: transcript target from meth_matrix
return_diff <- function(meth_matrix_element, ignore = TRUE, max=FALSE, diff=TRUE) {
  if (isTRUE(ignore)) {
    ignore_exp <- attr(meth_matrix_element, 'null_exp')
    if (length(ignore_exp) != 0) {
      ignore_idx <- which(colnames(meth_matrix_element) %in% ignore_exp)
      meth_matrix_element <- meth_matrix_element[,-ignore_idx]
    }
    
    if (diff==FALSE) { return(meth_matrix_element)} #exit here if diff==FALSE 
  }
  
  diff_idx <- apply(meth_matrix_element[,-1], MARGIN = 1, FUN=check_diff_row)
  
  if (max==TRUE) {
    
    max_idx <- apply(meth_matrix_element[,-1][diff_idx,], MARGIN=1, FUN=check_max_diff)
    return(apply(meth_matrix_element[diff_idx,][max_idx,]), 2, as.numeric)
    
  } else {
    return(apply(meth_matrix_element[diff_idx,]), 2, as.numeric)
  }
}
#check if methylation is different at base for each experiment. applies to object: row vector of coded methylation data.
check_diff_row <- function(row) {
  if (any(is.na(row))) {
    return(FALSE) } else {
      u <- unique(row)
      return(length(u) != 1)
    }
}
check_max_diff <- function(row) {
  return((is.element(0,row) && is.element(2, row)))
}
#splice has been added to provide a work around for reference coordinates, while the ref_coords from npTranscript are false.
methpipe <-function(wig_files, null_files, splice, target) {
  meth <- lapply(wig_files, read_wig)
  meth <- name_by_experiment(meth, wig_files)
  meth <- lapply(meth, adjust_experiment)
  meth_matrix <- matrixify(meth, target = 'Fraction')
  null_meth <- lapply(null_files, read_null)
  null_meth <- do.call(rbind, null_meth)
  for (i in names(meth_matrix)) {
    ignore_exp <- null_meth[which(null_meth$V1 == i), 'experiment']
    attr(meth_matrix[[i]], 'null_exp') <- ignore_exp
  }
  #correct bases to reference space
  #ref_coord_data <- read_refcoord(ref_file)
  ref_coord_data <- splice
  meth_matrix <- convert_to_ref_coords(meth_matrix, ref_coord_data)
  base <- read_baseline_wig(null_wig, 19999)
  meth_matrix <- lapply(meth_matrix, apply_baseline_meth, base, target)
  return(meth_matrix)
}
#read in null meth (baseline) data, and correct for offset to reference. applies to object: path to null experiment wig file
read_baseline_wig <- function(null_wig_file, ref_start) {
  null_meth <- read_wig(null_wig_file)
  null_meth <- lapply(null_meth, fill_missing)
  for (i in 1:length(null_meth)) {
    null_meth[[i]]$Base <- null_meth[[i]]$Base + ref_start
    null_meth[[i]]$coded <- meth_encode(null_meth[[i]])
  }
  return(null_meth)
}
