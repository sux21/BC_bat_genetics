# helper functions for calling alleles from OSIRIS tab-delimited output

#'@param x a vector
#'@return NA if length of x is zero
length_zero_to_na <- function(x) {
  ifelse(length(x) == 0, return(NA), return(x))
}

#'@param x a vector
#'@return TRUE if x is all NA, else FALSE
is_all_na <- function(x) {
  ! FALSE %in% is.na(x)
}

#'@param f path to osiris tab-delimited output
#'@return a data frame of osiris tab-delimited output
load_osiris_tab <- function(f) {
  read.table(f, sep = "\t", header = TRUE, na.strings = "")
}

#'@param x a character string separated by comma
#'@retrun a numeric vector
chr_to_num <- function(x) {
  if (is.na(x)) {
    out <- NA
  } else {
    out <- as.numeric(unlist(strsplit(x, ",")))
  }
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
#'@return a numeric vector x with some values replaced by their means
replace_with_mean <- function(x, val_to_replace, group) {
  index <- which(x %in% val_to_replace)
  index2 <- aggregate(index ~ group, FUN = max)$index
  
  x[index] <- NA
  
  val_means <- ( aggregate(val_to_replace ~ group, FUN = mean)$val_to_replace
                # |> round(digits = 1)
  )
  
  x[index2] <- val_means
  
  x <- x[!is.na(x)]
  
  return(x)
}

#'@param x a numeric vector
#'@param y a numeric vector
#'@return a numeric vector with y removed from x
remove_y_from_x <- function(x, y) {
  y_rm <- which(! x %in% y)
  return(x[y_rm])
}

#'@param x a numeric vector
#'@param y a numeric vector contained in x
#'@return integer index
find_y_in_x <- function(x, y) {
  stopifnot(!is.na(x))
  
  index <- which(x %in% y)
  return(index)
}

#'@param x numeric vector
#'@param max_diff the maximum difference between two adjacent values
#'@return integer index
find_adjacent_values <- function(x, max_diff) {
  index <- which(diff(x) <= max_diff)
  
  if (length(index) == 0) {
    index2 <- NA
  } else {
    index2 <- unique(sort(c(index,index+1)))
  }
  return(index2)
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
#'@return a numeric vector of adjacent value's group numbers
group_adjacent_values <- function(x, max_diff) {
  if (unique(is.na(x)) || length(x) == 0 || length(x) == 1) {
    return(NA)
  } 
  
  group <- c() 
  id = 1
  adj_val <- get_adjacent_values(x, max_diff)
  
  if (unique(is.na(adj_val))) {
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
  return(group)
}

#'@param x a numeric vector 
#'@param y a numeric vector with the same length as x
#'@param max_diff the maximum difference between two adjacent values
#'@return a list of numeric vectors
get_adjacent_x_with_max_y <- function(x, y, max_diff) {
  stopifnot(!is.na(x))
  
  adj_x <- get_adjacent_values(x, max_diff)
  adj_y <- y[find_adjacent_values(x, max_diff)]
  
  adj_group <- group_adjacent_values(x, max_diff)
  
  max_y <- aggregate(adj_y ~ adj_group, FUN = max)
  names(max_y) <- c("group", "max_y")
  
  out <- list(adj_x = adj_x
              , adj_y = adj_y
              , max_y = max_y$max_y)
  return(out)
} 

#'@param x a numeric vector
#'@param y a numeric vector with the same length as x
#'@param max_diff the maximum difference between two adjacent values
#'@return a list of two numeric vectors 
remove_adjacent_x_with_max_y <- function(x, y, max_diff) {
  stopifnot(!is.na(x))
  
  max_y <- get_adjacent_x_with_max_y(x, y, max_diff)
  max_index <- find_y_in_x(max_y[["adj_y"]], max_y[["max_y"]])
  
  adj_x <- max_y[["adj_x"]][-max_index]
  adj_y <- max_y[["adj_y"]][-max_index]
  
  out <- list(adj_x = adj_x
              , adj_y = adj_y)
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param A_dist maximum base pair distance to find non-template addition
#'@return a list of fragment size and intensity 
remove_A_addition <- function (f, r, A_dist) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  adj_frag_size <- get_adjacent_values(frag_size, A_dist) 
  adj_index <- find_adjacent_values(frag_size, A_dist)
  adj_frag_rfu <- frag_rfu[adj_index]
  
  if (unique(is.na(adj_frag_size)) ) { # no change if NA
    frag_size2 <- num_to_chr(frag_size) 
    frag_rfu2 <- num_to_chr(frag_rfu)
  } else {
    adj_frag2 <- remove_adjacent_x_with_max_y(frag_size, frag_rfu, A_dist) 
    
    # remove adjacent fragments with non-maximum intensity
    frag_size2_index <- which(! frag_size %in% adj_frag2$adj_x)
    frag_size2 <- frag_size[frag_size2_index] |> num_to_chr() 
    
    # remove corresponding intensity of these fragments
    frag_rfu2_index <- which(! frag_rfu %in% adj_frag2$adj_y)
    frag_rfu2 <- frag_rfu[frag_rfu2_index] |> num_to_chr() 
  } 
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@return a list of numeric vectors
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

#'@param spl_dat a list of sample fragment size and intensity
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@return a list of off-scale fragment size and intensity
get_off_scale <- function(spl_dat, min_off_rfu) {
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

#'@param spl_dat a list of sample fragment size and intensity
#'@param off a list of off-scale fragment size and intensity
#'@param ru_len repeat unit length, base pair distance to look for stutters
#'@return a list of off-scale fragment sizes with pull-up removed
remove_high_pull_up <- function(spl_dat, off, ru_len) {
  frag_size <- spl_dat[["size"]]
  frag_rfu <- spl_dat[["intensity"]]
  
  off_size <- off[["size"]]
  
  adj_size <- lapply(frag_size, get_adjacent_values, ru_len)
  
  # pull-up peaks have no adjacent peaks
  not_pull_up <- lapply(seq_along(off_size)
                          , function(i) {
                            off_size[[i]] %in% adj_size[[i]]
                          })
  off_size <- lapply(seq_along(off_size)
                            , function(i) {
                              off_size[[i]][not_pull_up[[i]]]
                            })
  off_size <- lapply(off_size, length_zero_to_na)
  return(off_size)
}

#'@param spl_dat a list of sample fragment size and intensity
#'@param off_size a list of off-scale fragment sizes
#'@param off_dist base pair distance around off-scale peak to find pull-up
#'@return a list of integer index of pull-up fragment sizes
find_pull_up <- function(spl_dat, off_size, off_dist) {
  num_colours <- length(spl_dat[["size"]])
  frag_size <- spl_dat[["size"]]
  
  pull_up_index <- vector(mode = "list", length = num_colours)
  
  for (i in seq_along(off_size)) {
    for (j in seq_along(off_size[[i]])) {
      for (k in seq_along(off_size)[-i]) {
        
        if (length(off_size[[i]]) == 0 || length(frag_size[[k]]) == 0) {
          next
        }
        
        left_side = off_size[[i]][j] - off_dist
        right_side = off_size[[i]][j] + off_dist
        
        pos <- which(frag_size[[k]] >= left_side & frag_size[[k]] <= right_side)
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

#'@param spl_dat a list of sample fragment size and intensity
#'@param pull_up_index a list of integer index of pull-up fragment sizes
#'@return a list of pull-up fragment sizes 
get_pull_up <- function(spl_dat, pull_up_index) {
  num_colours <- length(spl_dat[["size"]])
  frag_size <- spl_dat[["size"]]
  frag_rfu <- spl_dat[["intensity"]]
  
  pull_up_size <- lapply(seq_along(frag_size)
                         , function(x) {frag_size[[x]][pull_up_index[[x]]]}
  )
  
  pull_up_rfu <- lapply(seq_along(frag_rfu)
                        , function(x) {frag_rfu[[x]][pull_up_index[[x]]]}
  )
  
  pull_up <- list(pull_up_size, pull_up_rfu)
  names(pull_up) <- c("size", "intensity")
  return(pull_up)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param ru_len repeat unit length, base pair distance to look for stutters
#'@param off_dist base pair distance around off-scale peak to find pull-up
#'@return a character string with pull-up fragments removed
remove_pull_up <- function(f, s, min_off_rfu, ru_len, off_dist) {
  spl_dat <- get_sample_dat(f, s)
  num_colours <- length(spl_dat[["size"]])
  frag_size <- spl_dat[["size"]]
  frag_rfu <- spl_dat[["intensity"]]
  
  off <- get_off_scale(spl_dat, min_off_rfu)
  off_size2 <- remove_high_pull_up(spl_dat, off, ru_len)
  
  if (is_all_na(unlist(off_size2, use.names = FALSE))) { # no change if no off-scale
    frag_size2 <- lapply(frag_size, num_to_chr)
    frag_rfu2 <- lapply(frag_rfu, num_to_chr)
  } else {
    pull_up_index <- find_pull_up(spl_dat, off_size2, off_dist)
    # bug: if pull_up_index is NA, it will return NA but should return original fragments
    frag_size2 <- lapply(seq_along(frag_size)
                         , function(x) {frag_size[[x]][-pull_up_index[[x]]]})
    
    frag_rfu2 <- lapply(seq_along(frag_rfu)
                        , function(x) {frag_rfu[[x]][-pull_up_index[[x]]]})
    
    frag_size2 <- lapply(frag_size2, num_to_chr)
    frag_rfu2 <- lapply(frag_rfu2, num_to_chr)
  }
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param split_dist maximum base pair distance to find off-scale split peaks
#'@return OSIRIS tab-delimited file with split peaks replaced by mean
replace_split_by_mean <- function(f, r, min_off_rfu, split_dist) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  off_index <- which(frag_rfu >= min_off_rfu)
  off_size <- frag_size[off_index]
  off_rfu <- frag_rfu[off_index]
  
  split_index <- find_adjacent_values(off_size, split_dist)
  split_size <- off_size[split_index]
  split_rfu <- off_rfu[split_index]
  split_group <- group_adjacent_values(off_size, split_dist)
  
  if (unique(is.na(split_group))) {
    frag_size2 <- num_to_chr(frag_size)
    frag_rfu2 <- num_to_chr(frag_rfu)
  } else {
    frag_size2 <- ( replace_with_mean(frag_size, split_size, split_group)
                    |> round(digits = 1)
                    |> num_to_chr()
    )
    frag_rfu2 <- ( replace_with_mean(frag_rfu, split_rfu, split_group)
                   |> round(digits = 0)
                   |> num_to_chr())
  }
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param ploidy ploidy of organism, number of allelic peaks to be selected
#'@param peaks_ratio minimum ratio of light intensity of heterozygous peaks
#'@param noise_level no extra peaks above shortest allelic peak - this number
#'@return a list of characters of allelic peaks
allele_caller <- function(f, r, ploidy, peaks_ratio, noise_level) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  if (unique(is.na(frag_size)) ) { # no change if NA
    frag_size2 <- frag_size 
    frag_rfu2 <- frag_rfu
  } 
  
  rfu_dec <- sort(frag_rfu, decreasing = TRUE)
  
  if (length(rfu_dec) == 1) { # no change if only 1 peak
    frag_size2 <- num_to_chr(frag_size) 
    frag_rfu2 <- num_to_chr(frag_rfu)
  }
  
  noise_max_rfu <- ( rfu_dec[ploidy] - noise_level )
  peaks_above_noise <- frag_rfu[frag_rfu > noise_max_rfu]
  
  if (length(peaks_above_noise) == ploidy) {
    if (rfu_dec[ploidy]/rfu_dec[1] >= peaks_ratio) { # heterozygous
      frag_rfu2 <- rfu_dec[1:ploidy]
    } else { # homozygous
      frag_rfu2 <- rfu_dec[1]
    }
    frag_rfu2_index <- which(frag_rfu %in% frag_rfu2)
    frag_size2 <- frag_size[frag_rfu2_index]
  } else if (length(peaks_above_noise) > ploidy) {
    frag_rfu2 <- "too many peaks"
    frag_size2 <- "too many peaks"
  }
  frag_rfu2 <- num_to_chr(frag_rfu2)
  frag_size2 <- num_to_chr(frag_size2)
  
  out <- list(frag_size2, frag_rfu2)
  names(out) <- c("size", "intensity")
  return(out)
} 

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@return logical vector indicating which rows correspond to control samples
find_ctrl <- function(f, ctrl) {
  ctrl_sample <- grepl(pattern = paste(ctrl, collapse='|'), f$File.Name)
  return(ctrl_sample)
} 

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@return a numeric vector of fragment sizes in control samples
get_cont_size <- function(f, ctrl) {
  ctrl_sample <- find_ctrl(f, ctrl)
  
  cont_size <- f[ctrl_sample,]$Allele
  
  cont_size <- ( sapply(cont_size, chr_to_num) 
                      |> unlist(use.names = FALSE) 
  )
  cont_size <- cont_size[!is.na(cont_size)]
  return(cont_size)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
remove_cont <- function(f, ctrl){
  cont <- get_cont_size(f, ctrl)
  
  ctrl_sample <- find_ctrl(f, ctrl)
  
  sample <- f[!ctrl_sample,]$Allele |> sapply(chr_to_num, USE.NAMES = FALSE)
  
  sample_cont_rm <- lapply(sample, rm_val, cont)
  return(sample_cont_rm)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param ru_len repeat unit length, base pair distance to look for stutters
#'@param off_dist base pair distance around off-scale peak to find pull-up
#'@return OSIRIS tab-delimited file with pull up removed
remove_pull_up_all <- function(f, ctrl, min_off_rfu, ru_len, off_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  sample_name <- unique(sample$File.Name)
  
  f2 <- f
  
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
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param split_dist maximum base pair distance to find off-scale split peaks
#'@return OSIRIS tab-delimited file with split peaks replaced by its mean
replace_split_by_mean_all <- function(f, ctrl, min_off_rfu, split_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  
  f2 <- f
  
  row <- row.names(sample) |> as.numeric()
  
  for (i in row) {
    f2$Allele[i] <- replace_split_by_mean(f, i, min_off_rfu, split_dist)[["size"]] 
    f2$RFU[i] <- replace_split_by_mean(f, i, min_off_rfu, split_dist)[["intensity"]]
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param A_dist maximum base pair distance to find non-template addition
#'@return OSIRIS tab-delimited file with highest non-template addition peaks selected
remove_A_addition_all <- function(f, ctrl, A_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  
  f2 <- f
  
  row <- row.names(sample) |> as.numeric()
  
  for (i in row) {
    f2$Allele[i] <- remove_A_addition(f, i, A_dist)[["size"]] 
    f2$RFU[i] <- remove_A_addition(f, i, A_dist)[["intensity"]]
  }
  return(f2)
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, c("negative", "NEG", "Ladder")
#'@param ploidy ploidy of organism, number of allelic peaks to be selected
#'@param peaks_ratio minimum ratio of light intensity of heterozygous peaks
#'@param noise_level no extra peaks above shortest allelic peak - this number
#'@return OSIRIS tab-delimited file with allelic peak(s) selected
allele_caller_all <- function(f, ctrl, ploidy, peaks_ratio, noise_level) {
  ctrl_sample <- find_ctrl(f, ctrl)
  sample <- f[!ctrl_sample,]
  
  f2 <- f
  
  row <- row.names(sample) |> as.numeric()
  
  for (i in row) {
    f2$Allele[i] <- allele_caller(f
                                  , i
                                  , ploidy
                                  , peaks_ratio
                                  , noise_level)[["size"]] 
    f2$RFU[i] <- allele_caller(f
                               , i
                               , ploidy
                               , peaks_ratio
                               , noise_level)[["intensity"]]
  }
  print("done")
  #return(f2)
}