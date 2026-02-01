# helper functions for calling alleles from OSIRIS tab-delimited output

#'@param x a vector
#'@return NA if length of x is zero
length_zero_to_na <- function(x) {
  ifelse(length(x) == 0, return(NA), return(x))
}

#'@param f path to osiris tab-delimited output
#'@return a data frame of osiris tab-delimited output
load_osiris_tab <- function(f) {
  read.table(f, sep = "\t", header = TRUE, na.strings = "")
}

#'@param x a character string separated by comma or missing value (NA)
#'@retrun a numeric vector
chr_to_num <- function(x) {
  if (is.na(x)) {
    out <- NA
  } else {
    out <- as.numeric(unlist(strsplit(x, ",")))
  }
  return(out)
}

#'@param x a numeric vector or missing value (NA)
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
#'@param group group numbers of values to be replaced be means
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

#'@param x a numeric vector or missing value (NA)
#'@param max_diff the maximum difference between two adjacent values
#'@return integer index of adjacent values if there are adjacent values. Return NA if no adjacent values. 
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
#'@param y a numeric vector contained in x
#'@return integer index of elements of y in x
find_y_in_x <- function(x, y) {
  stopifnot(!is.na(x))
  
  index <- which(x %in% y)
  return(index)
}

#'@param x a numeric vector 
#'@param y a numeric vector with the same number of elements as x
#'@param max_diff the maximum difference between two adjacent values
#'@return a list of three numeric vectors of adjacent x and corresponding y and maximum y of each adjacent group
get_max_y_for_adjacent_x <- function(x, y, max_diff) {
  stopifnot(!is.na(x))
  
  adj_x <- get_adjacent_values(x, max_diff)
  
  adj_group <- group_adjacent_values(x, max_diff)
  
  adj_y <- y[find_adjacent_values(x, max_diff)]
    
  max_y_by_group <- aggregate(adj_y ~ adj_group, FUN = max)
  names(max_y_by_group) <- c("group", "max_y")
  
  out <- list(adj_x = adj_x
              , adj_y = adj_y
              , max_y = max_y_by_group$max_y)
  return(out)
} 

#'@param x a numeric vector
#'@param y a numeric vector with the same number of elements as x
#'@param max_diff the maximum difference between two adjacent values
#'@return a list of two numeric vectors of adjacent x and y which x with maximum y removed
remove_x_with_max_y <- function(x, y, max_diff) {
  stopifnot(!is.na(x))
  
  max_y <- get_max_y_for_adjacent_x(x, y, max_diff)
  
  max_index <- find_y_in_x(max_y[["adj_y"]], max_y[["max_y"]])
  
  adj_x <- max_y[["adj_x"]][-max_index]
  
  adj_y <- max_y[["adj_y"]][-max_index]
  
  out <- list(adj_x = adj_x
              , adj_y = adj_y)
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param max_diff the maximum difference between two adjacent values
#'@return a list of fragment size and intensity with non-template addition removed
remove_A_addition <- function (f, r, max_diff) {
  frag_size <- f$Allele[r] |> chr_to_num()
  frag_rfu <- f$RFU[r] |> chr_to_num()
  
  adj_frag_size <- get_adjacent_values(frag_size, max_diff) 
  adj_frag_rfu <- frag_rfu[find_adjacent_values(frag_size, max_diff)]
  
  if (unique(is.na(adj_frag_size))) { # no change if NA
    new_frag_size <- num_to_chr(frag_size) 
    new_frag_rfu <- num_to_chr(frag_rfu)
    
    out <- list(new_frag_size, new_frag_rfu)
    return(out)
  } else {
    adj_frag2 <- remove_x_with_max_y(frag_size, frag_rfu, max_diff) 
    # adjacent fragments with non-maximum intensity
    
    frag_size2 <- ( frag_size[! frag_size %in% adj_frag2$adj_x] 
                       |> num_to_chr() 
                       # remove adjacent fragments with non-maximum intensity
    )
    frag_rfu2 <- ( frag_rfu[! frag_inte %in% adj_frag2$adj_y] 
                       |> num_to_chr() 
                       # remove corresponding intensity of these fragments
    )
    out <- list(frag_size2, frag_rfu2)
    return(out)
  }  
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@return a list of numeric vectors of fragment size and intensity in each colour channel of a sample
get_sample_dat <- function(f, s) {
  osiris_subset <- subset(f, File.Name == s)
  
  dat <- split(osiris_subset, seq_len(nrow(osiris_subset)))
  
  sample_frag_size <- lapply(seq_along(dat)
                             , function (i) {chr_to_num(dat[[i]]$Allele)}
  )
  sample_frag_rfu <- lapply(seq_along(dat)
                            , function(i) {chr_to_num(dat[[i]]$RFU)}
  )
  
  sample_dat <- list(sample_frag_size, sample_frag_rfu)
  names(sample_dat) <- c("size", "intensity")
  
  return(sample_dat)
}

#'@param spl_dat a list of sample fragment size and intensity of each colour
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@return off-scale fragment size and intensity in each colour channel
get_off_scale <- function(spl_dat, min_off_rfu) {
  size <- spl_dat[["size"]]
  rfu <- spl_dat[["intensity"]]
  
  off_index <- lapply(rfu
                      , function(x) {which(x >= min_off_rfu)}) 
  
  off_size <- lapply(seq_along(size)
                     , function(x) {size[[x]][off_index[[x]]]})
  
  off_size <- lapply(off_size, length_zero_to_na)
  
  off_rfu <- lapply(seq_along(spl_dat[["intensity"]])
                     , function(x) {
                       spl_dat[["intensity"]][[x]][off_index[[x]]]
                     })
  
  off <- list(off_size, off_rfu)
  names(off) <- c("size", "intensity")
  
  return(off)
}

#'@param spl_dat a list of sample fragment size and intensity of each colour
#'@param off a list of off-scale fragment size and intensity
#'@param ru_len maximum distance to look for stutter around off-scale peak
#'@return a list of off-scale fragment sizes with pull-up removed in each colour channel
remove_high_pull_up <- function(spl_dat, off, ru_len) {
  sample_frag_size <- spl_dat[["size"]]
  sample_frag_rfu <- spl_dat[["intensity"]]
  
  off_size <- off[["size"]]
  
  adj_size <- lapply(sample_frag_size, get_adjacent_values, ru_len)
  
  have_adj_frag <- lapply(seq_along(off_size)
                          , function(i) {
                            off_size[[i]] %in% adj_size[[i]]
                          })
  off_scale_size2 <- lapply(seq_along(off_size)
                            , function(i) {
                              off_size[[i]][have_adj_frag[[i]]]
                            })
  return(off_scale_size2)
}

#'@param sample_dat a list of fragment size and intensity per sample by split_data_per_sample function
#'@param off_size a list of off-scale fragment sizes
#'@param off_dist maximum distance to the left and right the off-scale peaks to look for pull-up
#'@return a list of integer index of pull-up fragment sizes in each colour channel
find_pull_up <- function(sample_dat, off_size, off_dist) {
  frag_size <- sample_dat[["size"]]

  pull_up_index <- vector(mode = "list", length = 4)
  
  for (i in seq_along(off_size)) {
    for (j in seq_along(off_size[[i]])) {
      for (k in seq_along(off_size)[-i]) {
        if (length(off_size[[i]]) == 0 || length(frag_size[[k]]) == 0) {
          next
        }
        
        lower_bound = off_size[[i]][j] - off_dist
        upper_bound = off_size[[i]][j] + off_dist
        
        pos <- which(frag_size[[k]] >= lower_bound & frag_size[[k]] <= upper_bound)
        pull_up_index[[k]] <- ( pull_up_index[[k]] |> append(pos) |> unique() |> sort() )
        
        pull_up_index[[k]] <- length_zero_to_na(pull_up_index[[k]])
        # if (length(pull_up_index[[k]]) == 0) {
        #   pull_up_index[[k]] <- NA
        # }
      }
    }
  }
  return(pull_up_index)
}

#'@param per_sample_dat a list of fragment size and intensity per sample by split_data_per_sample function
#'@param pull_up_index a list of integer index of pull-up fragment sizes made by find_pull_up function
#'@return a list of pull-up fragment sizes in each colour channel
get_pull_up <- function(per_sample_dat, pull_up_index) {
  size_4channel <- per_sample_dat[["size"]]
  intensity_4channel <- per_sample_dat[["intensity"]]
    
  pull_up_size <- vector(mode = "list", length = 4)
  pull_up_intensity <- vector(mode = "list", length = 4)
  
  for (i in seq_along(size_4channel)) {
    pull_up_size[[i]] <- size_4channel[[i]][pull_up_index[[i]]]
    pull_up_intensity[[i]] <- intensity_4channel[[i]][pull_up_index[[i]]]
  }
  
  pull_up <- list(pull_up_size, pull_up_intensity)
  names(pull_up) <- c("size", "intensity")
  return(pull_up)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param dist_to_stutter maximum distance to look for adjacent stutter products from off-scale peak to tell it from strong pull-up
#'@param off_dist maximum distance to the left and right the off-scale peaks to look for pull-up
#'@return a character string with pull-up fragments removed
remove_pull_up <- function(f, s, min_off_rfu, dist_to_stutter, off_dist) {
  sample_dat <- get_sample_dat(f, s)
  
  off_size <- get_off_scale(sample_dat, min_off_rfu)
  off_size2 <- remove_high_pull_up(sample_dat, off_size, dist_to_stutter)
  
  pull_up_index <- find_pull_up(sample_dat, off_size2, off_dist)
  pull_up <- get_pull_up(sample_dat, pull_up_index)
  
  new_frag_size <- vector(mode = "list", length = 4)
  new_frag_intensity <- vector(mode = "list", length = 4)
  
  for (i in seq_along(sample_dat[["size"]])) {
    new_frag_size[[i]] <- remove_y_from_x(sample_dat[["size"]][[i]], pull_up[["size"]][[i]]) |> num_to_chr()
    new_frag_intensity[[i]] <- remove_y_from_x(sample_dat[["intensity"]][[i]], pull_up[["intensity"]][[i]]) |> num_to_chr()
  }
  
  out <- list(new_frag_size, new_frag_intensity)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param split_dist maximum base pair distance to find off-scale split peaks
average_split <- function(f, r, min_off_rfu, split_dist) {
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
#'@param p patterns of control sample name, e.g. p=c("negative", "NEG", "Ladder")
#'@return logical vector indicating which rows correspond to control samples
find_ctrl <- function(f, p) {
  ctrl_sample <- grepl(pattern = paste(p, collapse='|'), f$File.Name)
  return(ctrl_sample)
} 

#'@param f OSIRIS tab-delimited file
#'@param p patterns of control sample name, e.g. p=c("negative", "NEG", "Ladder")
#'@return a numeric vector of fragment sizes in control samples
get_cont_size <- function(f, p) {
  ctrl_sample <- find_ctrl(f, p)
  
  cont_size <- f[ctrl_sample,]$Allele
  
  cont_size <- ( sapply(cont_size, chr_to_num) 
                      |> unlist(use.names = FALSE) 
  )
  cont_size <- cont_size[!is.na(cont_size)]
  return(cont_size)
}

#'@param f OSIRIS tab-delimited file
#'@param p patterns of control sample name, e.g. p=c("negative", "NEG", "Ladder")
remove_cont <- function(f, p){
  cont <- get_cont_size(f, p)
  
  ctrl_sample <- find_ctrl(f, p)
  
  sample <- f[!ctrl_sample,]$Allele |> sapply(chr_to_num, USE.NAMES = FALSE)
  
  sample_cont_rm <- lapply(sample, rm_val, cont)
  return(sample_cont_rm)
}

select_alleles <- function(file, row_num) {
  # convert character string to numeric vector
  peaks <- as.numeric(unlist(strsplit(file$Allele[row_num], ",")))
  
  rfu <- as.numeric(unlist(strsplit(file$RFU[row_num], ",")))
  
  if (unique(is.na(peaks))) {
    # no change to the vector if the channel is missing or there is only one value
    peaks2 <- peaks
    
    rfu2 <- rfu
    
    output <- list(peaks2, rfu2)
    
  } else {
    # no change to the vector if there is only one value
    if (length(peaks) == 1) {
      peaks2 <- as.character(peaks)
      
      rfu2 <- as.character(rfu)
      
      output <- list(peaks2, rfu2)
      
    } else {
      # sort RFU in decreasing order, so the first and second values are the two highest peaks
      rfu_dec <- sort(rfu, decreasing = TRUE)
      
      if (rfu_dec[2]/rfu_dec[1] >= 0.5) {
        # select both peaks if RFU of the second peak is >= 50% of RFU of the highest peak
        rfu2 <- rfu_dec[1:2] 
        
        rfu2_index <- which(rfu %in% rfu2) 
        
        peaks2 <- peaks[rfu2_index]
        
        # convert back to character string
        peaks2 <- paste(peaks2, collapse = ",")
        
        rfu2 <- paste(rfu2, collapse = ",")
        
        # store output in a list
        output <- list(peaks2, rfu2)
        
      } else {
        # select the highest peak if not
        rfu2 <- rfu_dec[1]
        
        rfu2_index <- which(rfu %in% rfu2)
        
        peaks2 <- peaks[rfu2_index]
        
        # convert back to character string
        peaks2 <- paste(peaks2, collapse = ",")
        
        rfu2 <- paste(rfu2, collapse = ",")
        
        # store output in a list
        output <- list(peaks2, rfu2)
        
      }   
    }
  }
}

#'@param f OSIRIS tab-delimited file
#'@param ctrl patterns of control sample name, e.g. p=c("negative", "NEG", "Ladder")
#'@param min_off_rfu minimum intensity to find off-scale peaks
#'@param ru_len maximum distance to look for adjacent stutter products from off-scale peak to tell it from strong pull-up
#'@param off_dist maximum distance to the left and right the off-scale peaks to look for pull-up
#'@return OSIRIS tab-delimited file with pull up removed
remove_pull_up_all <- function(f, ctrl, min_off_rfu, ru_len, off_dist) {
  ctrl_sample <- find_ctrl(f, ctrl)
  exp_sample <- f[!ctrl_sample,]
  sample_name <- unique(exp_sample$File.Name)
  
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


