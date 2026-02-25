# helper functions for calling alleles from OSIRIS tab-delimited output

#'@param f OSIRIS tab-delimited file
#'@return logical value whether size and intensity vectors are equal lengths
check_data_vec_len <- function(f) {
  f2 <- f[,c("Allele", "RFU")] |> na.omit()  
  
  stopifnot(length(f2[,"Allele"]) == length(f2[,"RFU"]))
  
  frag_size <- lapply(f2[,"Allele"], chr_to_num)
  frag_rfu <- lapply(f2[,"RFU"], chr_to_num)
  
  is_len_same <- lapply(seq_along(frag_size), 
                        function(i) {
                          length(frag_size[[i]]) == length(frag_rfu[[i]])
                          })
  are_all_true <- all(unlist(is_len_same))
  return(are_all_true)
}

#'@param x a numeric vector
#'@return integer index of non-maximum values
which_not_max <- function(x) {
  which(!seq_along(x) %in% which.max(x))
}

#'@param x a vector
#'@return NA if length of x is zero or x is NULL
length_zero_to_na <- function(x) {
  ifelse(length(x) == 0 || is.null(x), return(NA), return(x))
}

#'@param x a vector
#'@param select_index an integer index to select elements at these positions
#'@return x if select_index is NA, else return subset of x
select_subset <- function(x, select_index) {
  if (all(is.na(select_index))) {
    return(x)
  } else {
    return(x[select_index])
  }
}

#'@param x a vector
#'@param remove_index an integer index to remove elements at these positions
#'@return x if remove_index is NA, else return subset of x
remove_subset <- function(x, remove_index) {
  if (all(is.na(remove_index))) {
    return(x)
  } else {
    return(x[-remove_index])
  }
}

#'@param x a character string of numbers separated by comma
#'@retrun a numeric vector
chr_to_num <- function(x) {
  if (is.na(x)) {
    out <- NA
  } else if (grepl(paste(c(LETTERS, letters), collapse = "|"), x)) {
    out <- x
  } else {
    out <- as.numeric(unlist(strsplit(x, ",")))
  }
  return(out)
}

#'@param x a character vector of numbers separated by comma
#'@return a numeric vector
chr_to_num_vec <- function(x) {
  out <- ( sapply(x, FUN = chr_to_num, USE.NAMES = FALSE) 
           |> unlist() 
           |> unique() 
           )
  return(out)
}

#'@param x a numeric vector
#'@return a character string separated by comma
num_to_chr <- function(x) {
  if (unique(is.na(x)) || length(x) == 0) {
    out <- NA
  } else {
    out <- paste(x, collapse = ",")
  }
  return(out)
}

#'@param x a numeric vector
#'@param val_to_replace values to be replaced by mean
#'@param group group numbers of values to be replaced
#'@return a numeric vector x with grouped values replaced by their means
replace_with_mean <- function(x, val_to_replace, group) {
  index <- which(x %in% val_to_replace)
  index2 <- aggregate(index ~ group, FUN = max)$index
  
  x[index] <- NA
  
  val_means <- aggregate(val_to_replace ~ group, FUN = mean)$val_to_replace
  
  x[index2] <- val_means
  
  x <- x[!is.na(x)]
  return(x)
}

#'@param x a vector
#'@param y a vector
#'@param new_chr a character string
#'@return the character string if y is not NA
replace_num_with_chr <- function(x, y, new_chr) {
  if (all(is.na(y))) {
    out <- x
  } else {
    out <- new_chr
  }
  return(out)
}

#'@param x numeric vector
#'@param max_diff the maximum difference between two adjacent values
#'@return integer index of adjacent values
find_adjacent_values <- function(x, max_diff) {
  index <- which(diff(x) <= max_diff)
  
  if (length(index) == 0) {
    adj_index <- NA
  } else {
    adj_index <- unique(sort(c(index,index+1)))
  }
  return(adj_index)
} 

#'@param x a numeric vector or missing value (NA)
#'@param max_diff the maximum difference between two adjacent values
#'@return a numeric vector of adjacent values
get_adjacent_values <- function(x, max_diff) {
  index <- find_adjacent_values(x, max_diff)

  if (length(index) == 0) {
    adj_val <- NA
  } else {
    adj_val <- x[index]
  }
  return(adj_val)
}

#'@param x a numeric vector
#'@param max_diff the maximum difference between two adjacent values
#'@return a list of values, integer index, and group number of adjacent values
group_adjacent_values <- function(x, max_diff) {
  if (all(is.na(x)) || length(x) == 0 || length(x) == 1) {
    return(NA)
  } 
  
  group <- c() 
  id = 1
  adj_val <- get_adjacent_values(x, max_diff)
  adj_index <- find_adjacent_values(x, max_diff)
  
  if (all(is.na(adj_val))) {
    return(NA)
  }
  
  for (i in seq_along(adj_val)) {
    if (i == 1 ) { # first value is always put in group 1
      group <- append(group, 1)
    }
    if (i > 1 && i < length(adj_val)) { # middle values
      if (adj_val[i] - adj_val[i-1] <= max_diff) {
        group <- append(group, id)
      } else {
        id = id + 1
        group <- append(group, id) 
      }
    }
    if (i == length(adj_val)) { # last value
      if (adj_val[i] - adj_val[i-1] <= max_diff) {
        group <- append(group, id)
      } else {
        id = id + 1
        group <- append(group, id) 
      }
    }
  }
  out <- list(value = adj_val, index = adj_index, group = group)
  return(out)
}

#'@param x a numeric vector 
#'@param y a numeric vector with the same length as x
#'@param max_diff the maximum difference between two adjacent values
#'@return integer index of adjacent values in x which x with non-maximum y removed
get_adjacent_x_with_non_max_y <- function(x, y, max_diff) {
  adj_x <- group_adjacent_values(x, max_diff)
  
  if (all(is.na(adj_x))) {
    return(NA)
  }
  
  adj_y <- y[adj_x$index]
  
  no_max_y <- aggregate(adj_y ~ adj_x$group, FUN = which_not_max) 
  no_max_y <- no_max_y$adj_y
  
  adj_index <- split(adj_x$index, factor(adj_x$group))
  
  adj_x_no_max_y <- ( lapply(seq_along(adj_index)
                             , function(i) {
                               select_subset(adj_index[[i]], no_max_y[[i]])
                             })
                      |> unlist())
  return(adj_x_no_max_y)
} 

#'@param f path to OSIRIS tab-delimited file
#'@return a data frame
load_osiris_tab <- function(f) {
  read.table(f, sep = "\t", header = TRUE, na.strings = "")
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@return logical vector indicating which rows correspond to control samples
find_ctrl <- function(f, ctrl) {
  ctrl_spl <- grepl(pattern = paste(ctrl, collapse='|'), f$File.Name)
  return(ctrl_spl)
} 

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param A_dist maximum base pair distance to find non-template addition
#'@return a list of fragment size and intensity 
get_high_A <- function (f, r, A_dist) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  frag_to_remove <- get_adjacent_x_with_non_max_y(frag_size, frag_rfu, A_dist)
  
  frag_size2 <- remove_subset(frag_size, frag_to_remove) |> num_to_chr()
  frag_rfu2 <- remove_subset(frag_rfu, frag_to_remove) |> num_to_chr()
  
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param A_dist maximum base pair distance to find non-template addition
#'@return a data frame
get_high_A_all <- function(f, ctrl, A_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  
  f2 <- f
  
  row <- row.names(sample) |> as.numeric()
  
  for (i in row) {
    f2$Allele[i] <- get_high_A(f, i, A_dist)[["size"]] 
    f2$RFU[i] <- get_high_A(f, i, A_dist)[["intensity"]]
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@return a list of sample fragment size and intensity
get_sample_dat <- function(f, s) {
  osiris_subset <- subset(f, File.Name == s)
  dat <- split(osiris_subset, seq_len(nrow(osiris_subset)))
  
  frag_size <- lapply(seq_along(dat)
                      , function (i) {chr_to_num(dat[[i]]$Allele)}
  )
  frag_rfu <- lapply(seq_along(dat)
                     , function(i) {chr_to_num(dat[[i]]$RFU)}
  )
  spl_dat <- list(frag_size, frag_rfu)
  names(spl_dat) <- c("size", "intensity")
  return(spl_dat)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@return a list of off-scale fragment size and intensity
get_off_scale <- function(f, s, min_off_rfu) {
  spl_dat <- get_sample_dat(f, s)
  size <- spl_dat[["size"]]
  rfu <- spl_dat[["intensity"]]
  
  off_index <- lapply(rfu
                      , function(x) {which(x >= min_off_rfu)}) 
  
  off_size <- lapply(seq_along(size)
                     , function(x) {size[[x]][off_index[[x]]]})
  
  off_rfu <- lapply(seq_along(rfu)
                     , function(x) {rfu[[x]][off_index[[x]]]})
  
  off_size <- lapply(off_size, length_zero_to_na)
  off_rfu <- lapply(off_rfu, length_zero_to_na)
  
  off <- list(off_size, off_rfu)
  names(off) <- c("size", "intensity")
  return(off)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param ru_len repeat unit length, base pair distance to look for stutters
#'@return a list of off-scale fragment sizes with pull-up removed
get_off_scale2 <- function(f, s, min_off_rfu, ru_len) {
  spl_dat <- get_sample_dat(f, s)
  frag_size <- spl_dat[["size"]]
  frag_rfu <- spl_dat[["intensity"]]
  
  off <- get_off_scale(f, s, min_off_rfu)
  off_size <- off$size
  
  adj_size <- lapply(frag_size, get_adjacent_values, ru_len)
  
  # pull-up peaks have no stutter
  not_pull_up <- lapply(seq_along(off_size)
                        , function(i) {
                          off_size[[i]] %in% adj_size[[i]]
                        })
  off_size2 <- lapply(seq_along(off_size)
                      , function(i) {
                        off_size[[i]][not_pull_up[[i]]]
                      })
  off_size2 <- lapply(off_size2, length_zero_to_na)
  return(off_size2)
}


#'@param spl_dat a list of sample fragment size and intensity
#'@param off_size a list of off-scale fragment sizes
#'@param off_dist base pair distance around off-scale peak to find pull-up
#'@return a list of integer index of pull-up fragment sizes
find_pull_up <- function(f, s, min_off_rfu, ru_len, off_dist) {
  spl_dat <- get_sample_dat(f, s)
  frag_size <- spl_dat$size
  num_colours <- length(frag_size)
  
  off_size <- get_off_scale2(f, s, min_off_rfu, ru_len)
  
  pull_up_index <- vector(mode = "list", length = num_colours)
  
  for (i in seq_along(off_size)) { # off-scale channel, i = 1 to 4
    for (j in seq_along(off_size[[i]])) { # off-scale size in a channel
      for (k in seq_along(off_size)[-i]) { # not off-scale channel
        
        if (length(off_size[[i]]) == 0 || 
            length(frag_size[[k]]) == 0) {
          next
        }
        
        left_side = off_size[[i]][j] - off_dist
        right_side = off_size[[i]][j] + off_dist
        
        pos <- which(frag_size[[k]] >= left_side & 
                       frag_size[[k]] <= right_side)
        
        pull_up_index[[k]] <- ( pull_up_index[[k]] 
                                |> append(pos) 
                                |> unique() 
                                |> sort() 
                                )
        
        pull_up_index[[k]] <- length_zero_to_na(pull_up_index[[k]])
      }
    }
  }
  return(pull_up_index)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param ru_len repeat unit length, base pair distance to look for stutters
#'@param off_dist base pair distance around off-scale peak to find pull-up
#'@return fragment size and intensity with pull-up removed
remove_pull_up <- function(f, s, min_off_rfu, ru_len, off_dist) {
  spl_dat <- get_sample_dat(f, s)
  frag_size <- spl_dat$size
  frag_rfu <- spl_dat$intensity
  
  pull_up_index <- find_pull_up(f, s, min_off_rfu, ru_len, off_dist)
  
  frag_size2 <- lapply(seq_along(frag_size)
                       , function(i) {
                         remove_subset(frag_size[[i]], pull_up_index[[i]])
                       })
  
  frag_rfu2 <- lapply(seq_along(frag_rfu)
                      , function(i) {
                        remove_subset(frag_rfu[[i]], pull_up_index[[i]])
                      })
  
  frag_size2 <- lapply(frag_size2, num_to_chr)
  frag_rfu2 <- lapply(frag_rfu2, num_to_chr)
  
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param ru_len repeat unit length, base pair distance to look for stutters
#'@param off_dist base pair distance around off-scale peak to find pull-up
#'@return a data frame
remove_pull_up_all <- function(f, ctrl, min_off_rfu, ru_len, off_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  sample_name <- unique(sample$File.Name)
  
  f2 <- f[,c("File.Name", "Locus", "Allele", "RFU")]
  
  for (i in seq_along(sample_name)) {
    current_sample <- sample_name[i]
    
    out <- remove_pull_up(f, current_sample, min_off_rfu, ru_len, off_dist)
    
    current_rows <- which(f2$File.Name %in% sample_name[i])
    
    for (j in seq_along(current_rows)) {
      r <- current_rows[j]
      
      f2[r,]$Allele <- out[["size"]][[j]]
      f2[r,]$RFU <- out[["intensity"]][[j]]
    }
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param split_dist maximum base pair distance to find off-scale split peaks
#'@return fragment size and intensity with split off-scale values replaced by mean
get_mean_split <- function(f, r, min_off_rfu, split_dist) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  off_index <- which(frag_rfu >= min_off_rfu)
  off_size <- frag_size[off_index]
  off_rfu <- frag_rfu[off_index]
  
  split <- group_adjacent_values(off_size, split_dist)

  if (all(is.na(split))) { # no change if no off-scale peaks
    frag_size2 <- num_to_chr(frag_size)
    frag_rfu2 <- num_to_chr(frag_rfu)
  } else {
    split_index <- split$index
    split_group <- split$group
    
    split_size <- off_size[split_index]
    split_rfu <- off_rfu[split_index]
    
    frag_size2 <- ( replace_with_mean(frag_size, split_size, split_group)
                    |> round(digits = 1)
                    |> num_to_chr())
    frag_rfu2 <- ( replace_with_mean(frag_rfu, split_rfu, split_group)
                   |> round(digits = 0)
                   |> num_to_chr())
  }
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param split_dist maximum base pair distance to find off-scale split peaks
#'@return a data frame
get_mean_split_all <- function(f, ctrl, min_off_rfu, split_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  
  f2 <- f[,c("File.Name", "Locus", "Allele", "RFU")]
  
  row <- row.names(sample) |> as.numeric()
  
  for (i in row) {
    f2$Allele[i] <- get_mean_split(f, i, min_off_rfu, split_dist)[["size"]] 
    f2$RFU[i] <- get_mean_split(f, i, min_off_rfu, split_dist)[["intensity"]]
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param ploidy ploidy of organism, number of allelic peaks to be selected
#'@param peak_ratio minimum ratio of light intensity of heterozygous peaks
#'@param stutter_ratio maximum ratio between stutter and allelic peak
#'@param noise_level no extra peaks allowed above shortest allelic peak minus this number
#'@return a list of characters of allelic peaks
allele_caller <- function(f, r, ploidy, peak_ratio, stutter_ratio, noise_level) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  # no change if NA or only 1 peak
  if (unique(is.na(frag_size)) || length(frag_rfu) == 1) {
    frag_size2 <- num_to_chr(frag_size) 
    frag_rfu2 <- num_to_chr(frag_rfu)
  } 
  
  rfu_dec <- sort(frag_rfu, decreasing = TRUE, index.return = TRUE)
  size_dec <- frag_size[rfu_dec$ix]
    
  noise_max_rfu <- ( rfu_dec$x[ploidy] - noise_level )
  peaks_above_noise <- frag_rfu[frag_rfu > noise_max_rfu]
  # 
  # allele_out <- function(rfu_dec, frag_size, ploidy) {
  #   frag_rfu2 <- rfu_dec$x[1:ploidy] |> num_to_chr()
  #   frag_rfu2_index <- rfu_dec$ix[1:ploidy]
  #   frag_size2 <- frag_size[frag_rfu2_index] |> num_to_chr()
  #   
  #   out <- list(frag_size2, frag_rfu2)
  #   names(out) <- c("size", "intensity")
  #   return(out)
  # }
  
  # select the number of alleles equal as expected from ploidy
  if (length(peaks_above_noise) == ploidy) {
    if ( rfu_dec$x[ploidy]/rfu_dec$x[1] >= peak_ratio &&
         size_dec[1] - size_dec[ploidy] < 0) { # heterozygous
      
      frag_rfu2 <- rfu_dec$x[1:ploidy] |> num_to_chr()
      frag_rfu2_index <- rfu_dec$ix[1:ploidy]
      frag_size2 <- frag_size[frag_rfu2_index] |> num_to_chr()
      
    } else if (rfu_dec$x[ploidy]/rfu_dec$x[1] >= peak_ratio &&
               size_dec[1] - size_dec[ploidy] > 0 && 
               size_dec[1] - size_dec[ploidy] > ru_len) { # heterozygous
      
      frag_rfu2 <- rfu_dec$x[1:ploidy] |> num_to_chr()
      frag_rfu2_index <- rfu_dec$ix[1:ploidy]
      frag_size2 <- frag_size[frag_rfu2_index] |> num_to_chr()
      
    } else if (rfu_dec$x[ploidy]/rfu_dec$x[1] >= peak_ratio && 
               size_dec[1] - size_dec[ploidy] > 0 &&
               size_dec[1] - size_dec[ploidy] <= ru_len &&
               rfu_dec$x[ploidy]/rfu_dec$x[1] >= stutter_ratio) { # heterozygous
      
      frag_rfu2 <- rfu_dec$x[1:ploidy] |> num_to_chr()
      frag_rfu2_index <- rfu_dec$ix[1:ploidy]
      frag_size2 <- frag_size[frag_rfu2_index] |> num_to_chr()
      
    } else { # homozygous
      frag_rfu2 <- rfu_dec$x[1] |> num_to_chr()
      frag_rfu2_index <- rfu_dec$ix[1]
      frag_size2 <- frag_size[frag_rfu2_index] |> num_to_chr()
    }
  } 
  
  # too many peaks than expected from ploidy
  if (length(peaks_above_noise) > ploidy) {
    frag_rfu2 <- "too many peaks"
    frag_size2 <- "too many peaks"
  }
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
} 

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param ploidy ploidy of organism, number of allelic peaks to be selected
#'@param peaks_ratio minimum ratio of light intensity of heterozygous peaks
#'@param stutter_ratio maximum ratio between stutter and allelic peak
#'@param noise_level no extra peaks above shortest allelic peak minus this number
#'@return a data frame
allele_caller_all <- function(f, ctrl, ploidy, peak_ratio, stutter_ratio, noise_level) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  
  f2 <- f[,c("File.Name", "Locus", "Allele", "RFU")]
  
  row <- row.names(sample) |> as.numeric()
  
  for (i in row) {
    f2$Allele[i] <- allele_caller(f
                                  , i
                                  , ploidy
                                  , peak_ratio
                                  , stutter_ratio
                                  , noise_level)[["size"]] 
    f2$RFU[i] <- allele_caller(f
                               , i
                               , ploidy
                               , peak_ratio
                               , stutter_ratio
                               , noise_level)[["intensity"]]
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@return a list with fragment sizes in control
get_contamination <- function(f, ctrl){
  ctrl_sample <- find_ctrl(f, ctrl)
  ctrl_dat <- f[ctrl_sample,]
  cont <- split(ctrl_dat$Allele, factor(ctrl_dat$Locus))
  cont <- lapply(cont, chr_to_num_vec)
  return(cont)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param cont_dist base pair distance around peaks in control to be also treated as contamination
#'@return a list of integer index of where is contamination in sample
find_contamination <- function(f, s, ctrl, cont_dist) {
  spl_dat <- get_sample_dat(f, s)
  frag_size <- spl_dat$size
  num_colours <- length(frag_size)
  
  cont <- get_contamination(f, ctrl)
  
  cont_index <- vector(mode = "list", length = num_colours)
  
  for (i in seq_along(cont)) { # contamination colour channels
    for (j in seq_along(cont[[i]])) {
      current_cont <- cont[[i]][j] # each contamination fragment in a channel
      
      if (is.na(current_cont) || !is.numeric(frag_size[[i]])) { next }
      
      left_side <- ( current_cont - cont_dist )
      right_side <- ( current_cont + cont_dist )
      
      pos <- which(frag_size[[i]] >= left_side & frag_size[[i]] <= right_side)
      
      cont_index[[i]] <- ( cont_index[[i]] 
                           |> append(pos) 
                           |> unique() 
                           |> sort() )
    }
    cont_index[[i]] <- length_zero_to_na(cont_index[[i]])
  }
  return(cont_index)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param cont_dist base pair distance around peaks in control to be also treated as contamination
#'@return character string of fragment size or word "contaminated"
remove_contamination <- function(f, s, ctrl, cont_dist) {
  spl_dat <- get_sample_dat(f, s)
  frag_size <- spl_dat$size
  frag_rfu <- spl_dat$intensity
  
  cont_index <- find_contamination(f, s, ctrl, cont_dist)
  
  frag_size2 <- lapply(seq_along(frag_size)
                       , function(i) {
                         replace_num_with_chr(frag_size[[i]]
                                              , cont_index[[i]]
                                              , "contaminated")
                       })
  
  frag_rfu2 <- lapply(seq_along(frag_rfu)
                      , function(i) {
                        replace_num_with_chr(frag_rfu[[i]]
                                             , cont_index[[i]]
                                             , "contaminated")
                      })
  
  frag_size2 <- lapply(frag_size2, num_to_chr)
  frag_rfu2 <- lapply(frag_rfu2, num_to_chr)
  
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param cont_dist base pair distance around peaks in control to be also treated as contamination
#'@return a data frame
remove_contamination_all <- function(f, ctrl, cont_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  sample_name <- unique(sample$File.Name)
  
  f2 <- f[,c("File.Name", "Locus", "Allele", "RFU")]
  
  for (i in seq_along(sample_name)) {
    current_sample <- sample_name[i]
    
    out <- remove_contamination(f, current_sample, ctrl, cont_dist)
    
    current_rows <- which(f2$File.Name %in% sample_name[i])
    
    for (j in seq_along(current_rows)) {
      r <- current_rows[j]
      
      f2[r,]$Allele <- out[["size"]][[j]]
      f2[r,]$RFU <- out[["intensity"]][[j]]
    }
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@return a data frame in wide format
long_to_wide_data <- function(f, ctrl) {
  ctrl_sample <- find_ctrl(f, ctrl)
  f <- f[!ctrl_sample,]
  
  f$Locus <- gsub("Channel1-1", "BLUE",
                  gsub("Channel2-2", "GREEN",
                       gsub("Channel3-3", "YELLOW",
                            gsub("Channel4-4", "RED", f$Locus))))
  
  f2 <- ( f
          |> reshape(drop = c("RFU")
                     , idvar = "File.Name"
                     , timevar = c("Locus")
                     , direction = "wide"))
  return(f2)
} 