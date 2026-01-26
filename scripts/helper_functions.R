# helper functions for processing OSIRIS tab-delimited output

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
  
#'@param x a numeric vector or missing value (NA)
#'@param d the maximum difference between two values to be called adjacent values
#'@return integer index of adjacent values if there are adjacent values. Return NA if no adjacent values. 
find_adjacent_values <- function(x, d) {
  index <- which(diff(x) <= d)
  
  if (length(index) == 0) { 
    index2 <- NA
  } else {
    index2 <- unique(sort(c(index,index+1)))
  }
  
  return(index2)
} 

#'@param x a numeric vector or missing value (NA)
#'@param d the maximum difference between two values to be called adjacent values
#'@return adjacent values which difference is less than or equal to d 
get_adjacent_values <- function(x, d) {
  index <- find_adjacent_values(x, d)
  
  if (length(index) == 0) {
    adj_val <- NA
  } else {
    adj_val <- x[index]
  }
  
  return(adj_val)
}

#'@param x a numeric vector
#'@param d the maximum difference between two values to be called adjacent values
#'@return adjacent values which difference is less than or equal to d 
#'@return integer vector of group id for adjacent values
group_adj_values <- function(x, d) {
  stopifnot(!is.na(x))
  
  group <- c() 
  
  id = 1
  
  adj_val <- get_adjacent_values(x, d)
  
  for (i in seq_along(adj_val)) {
    
    if (i == 1 ) { # first value is always put in group 1
      group <- append(group, 1)
    }
    
    if (i > 1 && i < length(adj_val)) { # middle values
      if (adj_val[i] - adj_val[i-1] <= d) {
        group <- append(group, id)
      } else {
        id = id + 1
        group <- append(group, id) 
      }
    }
    
    if (i == length(adj_val)) { # last value
      if (adj_val[i] - adj_val[i-1] <= d) {
        group <- append(group, id)
      } else {
        id = id + 1
        group <- append(group, id) 
      }
    }
  }
  
  out <- list(adjacent_values = adj_val, group = group)
  return(out)
}

#'@param x a numeric vector
#'@param i a numeric vector contained in x
#'@return integer index of where is i in x
get_value_position <- function(x, i) {
  stopifnot(!is.na(x))
  
  index <- which(x %in% i)
  return(index)
}

#'@param x a numeric vector 
#'@param y a numeric vector with the same number of elements as x
#'@param d the maximum difference between two values to be called adjacent values
#'@return integer vector of group id for adjacent x values
#'@return integer index of maximum y values in each group 
find_adj_x_with_max_y <- function(x, y, d) {
  stopifnot(!is.na(x))
  
  adj_x <- group_adj_values(x, d)
  
  adj_y <- y[find_adjacent_values(x, d)]
    
  max_y_by_group <- aggregate(adj_y ~ adj_x[[2]], FUN = max)
  names(max_y_by_group) <- c("group", "max_y")
  
  out <- list(adj_x = adj_x[[1]]
              , adj_y = adj_y
              , max_y = max_y_by_group$max_y)
  
  return(out)
} 

#'@param x a numeric vector
#'@param y a numeric vector with the same number of elements as x
#'@param d the maximum difference between two values to be called adjacent values
#'@return numeric vector of adjacent x after removing x with maximum y values  
#'@return numeric vector of corresponding y for x above
find_adj_x_with_non_max_y <- function(x, y, d) {
  stopifnot(!is.na(x))
  
  max_y <- find_adj_x_with_max_y(x, y, d)
  
  max_index <- get_value_position(max_y[["adj_y"]], max_y[["max_y"]])
  
  adj_x_non_max <- max_y[["adj_x"]][-max_index]
  
  adj_y_non_max <- max_y[["adj_y"]][-max_index]
  
  out <- list(adj_x_non_max = adj_x_non_max
              , adj_y_non_max = adj_y_non_max)
  
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param d the maximum size difference between two peaks to be called adjacent peaks
#'@return a character string with adjacent fragments with non-maximum intensity removed
remove_non_max_adjacent_frag <- function (f, r, d) {
  frag_size <- f$Allele[r] |> chr_to_num() # get fragment size and intensity
  
  frag_intensity <- f$RFU[r] |> chr_to_num()
  
  adj_frag_size <- get_adjacent_values(frag_size, d) # get adjacent fragment sizes
  
  adj_frag_intensity <- frag_intensity[find_adjacent_values(frag_size, d)]
  
  if (unique(is.na(adj_frag_size))) { # no adjacent fragments
    
    new_frag_size <- num_to_chr(frag_size) 
    
    new_frag_intensity <- num_to_chr(frag_intensity)
    
    out <- list(new_frag_size, frag_intensity)
    return(out)
    
  } else { # remove adjacent fragments with non-maximum intensity
    
    non_max_adj_frag <- find_adj_x_with_non_max_y(frag_size, frag_intensity, d) # find adjacent fragments with non-maximum intensity
    
    new_frag_size <- frag_size[! frag_size %in% non_max_adj_frag$adj_x_non_max] |> num_to_chr() # remove adjacent fragments with non-maximum intensity
    
    new_frag_intensity <- frag_intensity[! frag_intensity %in% non_max_adj_frag$adj_y_non_max] |> num_to_chr() # remove corresponding intensity of these fragments
    
    out <- list(new_frag_size, new_frag_intensity)
    return(out)
  }  
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@return a list of fragment size and intensity in each colour channel
split_data_per_sample <- function(f, s) {
  osiris_subset <- subset(osiris_out, File.Name == s)
  
  dat <- split(osiris_subset, seq_len(nrow(osiris_subset)))
  
  size_4channel <- lapply(seq_along(dat), function (i) {chr_to_num(dat[[i]]$Allele)})
  intensity_4channel <- lapply(seq_along(dat), function(i) {chr_to_num(dat[[i]]$RFU)})
  
  per_sample_dat <- list(size_4channel, intensity_4channel)
  names(per_sample_dat) <- c("size", "intensity")
  
  return(per_sample_dat)
}

#'@param per_sample_dat a list of fragment size and intensity per sample by split_data_per_sample function
#'@param min_rfu minimum intensity to find off-scale peaks
#'@return a list of off-scale fragment sizes in each colour channel
get_off_scale_size <- function(per_sample_dat, min_rfu) {
  size_4channel <- per_sample_dat[["size"]]
  
  intensity_4channel <- per_sample_dat[["intensity"]]
  
  off_scale_index <- lapply(seq_along(intensity_4channel), function(i) {which(intensity_4channel[[i]] >= min_rfu)}) 
  
  off_scale_size <- lapply(seq_along(size_4channel), function(i) {size_4channel[[i]][off_scale_index[[i]]]})
  
  return(off_scale_size)
}

#'@param per_sample_dat a list of fragment size and intensity per sample by split_data_per_sample function
#'@param off_scale_size a list made by get_off_scale_size function
#'@param dist_to_stutter maximum distance to look for adjacent stutter products from off-scale peak to tell it from strong pull-up
#'@return a list of off-scale fragment sizes with pull-up removed in each colour channel
remove_strong_pull_up <- function(per_sample_dat, off_scale_size, dist_to_stutter) {
  size_4channel <- per_sample_dat[["size"]]
  intensity_4channel <- per_sample_dat[["intensity"]]
  
  adj_size_4channel <- lapply(size_4channel, get_adjacent_values, dist_to_stutter)
  
  have_adj_frag <- lapply(seq_along(off_scale_size), function(i) {off_scale_size[[i]] %in% adj_size_4channel[[i]]})
  
  off_scale_size2 <- lapply(seq_along(off_scale_size), function(i) {off_scale_size[[i]][have_adj_frag[[i]]]})
  
  return(off_scale_size2)
}

#'@param per_sample_dat a list of fragment size and intensity per sample by split_data_per_sample function
#'@param off_scale_size2 a list made by the remove_strong_pull_up function
#'@param dist_to_pull_up maximum distance to the left and right the off-scale peaks to look for pull-up
#'@return a list of integer index of pull-up fragment sizes in each colour channel
find_pull_up <- function(per_sample_dat, off_scale_size2, dist_to_pull_up) {
  size_4channel <- per_sample_dat[["size"]]

  pull_up_index <- vector(mode = "list", length = 4)
  
  for (i in seq_along(off_scale_size2)) {
    for (j in seq_along(off_scale_size2[[i]])) {
      for (k in seq_along(off_scale_size2)[-i]) {
        if (length(off_scale_size2[[i]]) == 0 | length(size_4channel[[k]]) == 0) {next}
        
        lower_bound = off_scale_size2[[i]][j] - dist_to_pull_up
        upper_bound = off_scale_size2[[i]][j] + dist_to_pull_up
        
        pos <- which(size_4channel[[k]] >= lower_bound & size_4channel[[k]] <= upper_bound)
        pull_up_index[[k]] <- ( pull_up_index[[k]] |> append(pos) |> unique() |> sort() )
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
  
  pull_up_size <- vector(mode = "list", length = 4)
  
  for (i in seq_along(size_4channel)) {
    pull_up_size[[i]] <- size_4channel[[i]][pull_up_index[[i]]]
  }
  return(pull_up_size)
}

#'@param f OSIRIS tab-delimited file
#'@param s sample name
#'@param min_rfu minimum intensity to find off-scale peaks
#'@param dist_to_stutter maximum distance to look for adjacent stutter products from off-scale peak to tell it from strong pull-up
#'@param dist_to_pull_up maximum distance to the left and right the off-scale peaks to look for pull-up
#'@return a character string with pull-up fragments removed
remove_pull_up <- function(f, s, min_rfu, dist_to_stutter, dist_to_pull_up) {
  per_sample_dat <- split_data_per_sample(f, s)
  
  off_scale <- get_off_scale_size(per_sample_dat, min_rfu)
  
  off_scale2 <- remove_strong_pull_up(per_sample_dat, off_scale, dist_to_stutter)
  
  pull_up_index <- find_pull_up(per_sample_dat, off_scale2, dist_to_pull_up)
  
  new_frag_size <- vector(mode = "list", length = 4)
  new_frag_intensity <- vector(mode = "list", length = 4)
  
  for (i in seq_along(per_sample_dat[["size"]])) {
    new_frag_size[[i]] <- per_sample_dat[["size"]][[i]][-pull_up_index[[i]]] |> num_to_chr()
    
    new_frag_intensity[[i]] <- per_sample_dat[["intensity"]][[i]][-pull_up_index[[i]]] |> num_to_chr()
  }
  
  out <- list(new_frag_size, new_frag_intensity)
  names(out) <- c("size", "intensity")
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param p patterns of negative control sample name, e.g. p=c("negative", "NEG").Set p=NA if no negative control.
#'@return logical vector indicating which rows correspond to negative control sample
find_negative_ctrl <- function(f, p) {
  neg_sample <- grepl(pattern = paste(p, collapse='|'), f$File.Name)
  return(neg_sample)
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