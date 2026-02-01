# record allele sizes

source("scripts/allele_calling_helper.R")


# set parameters
num_channels <- 4 # number of colour channels in OSIRIS output. This may not equal to how many colours used
control <- c("negative", "NEG", "Ladder") # negative control and size standard
min_off_rfu <- 27000 # minimum RFU to find an off-scale peak
ru_len <- 4 # maximum microsatellite repeat unit length
pull_up_dist <- 0.5 # base pair distance around off-scale peak to find pull-up
max_dist_A_addition = 1.5 # maximum base pair distance to find non-template addition
max_dist_split_peaks = 1 # maximum base pair distance to find off-scale split peaks


# load data
path_to_files <- "analysis/microsat_data/"
osiris_files <- list.files(path_to_files, pattern = ".tab")
osiris_paths <- paste0(path_to_files, osiris_files)

osiris_out <- lapply(osiris_paths, load_osiris_tab)
names(osiris_out) <- osiris_files


# test the functions on one of the data
osiris_test <- osiris_out[["fragment_analysis_mobix_2024-02-28_1.tab"]]

# remove pull-up 
osiris_test2 <- remove_pull_up_all(osiris_test
                                   , control
                                   , min_off_rfu
                                   , ru_len
                                   , pull_up_dist)

# average off-scale split peaks

average_split(osiris_test2,54, min_off_rfu, max_dist_split_peaks)



# remove fragments found in control samples (negative control, size standard) 
# which shouldn't have any fragments.

# Note: this will remove fragments in any multiplexes. 
# For example, fragments in multiplex 3 of a sample will be removed 
# if these fragments are found in multiplex 1 of negative control.
sample <- remove_cont(osiris_test, control) 









min_rfu = 27000
dist_to_stutter = 4
dist_to_pull_up = 0.5
per_sample <- remove_pull_up(osiris_out, "E01_0728-12--M1_004_5759", 
                             min_rfu, dist_to_stutter, dist_to_pull_up)
print(per_sample)

# To do: loop through all samples in the data frame for removing pull-up
# To do: write another function for adjacent off-scale peaks, replace with mean of these adjacent peaks
# To do: reformat the select_alleles function




# 
# pull_up_intensity <- vector(mode = "list", length = 4)
# 
# for (i in seq_along(intensity_4channel)) {
#   pull_up_intensity[[i]] <- intensity_4channel[[i]][index[[i]]]
# }
# print(pull_up_intensity)

# # remove adjacent peaks which cannot be pull-up
# max_dist_adj_peaks = 1.5
# 
# adj_index <- lapply(pull_up, find_adjacent_values, max_dist_adj_peaks)
# 
# 
# 
# pull_up2 <- vector(mode = "list", length = 4)
# 
# for (i in seq_along(pull_up)) {
#   if (unique(is.na(adj_index[[i]]))) {
#     pull_up2[[i]] <- pull_up[[i]]
#   } else{
#     pull_up2[[i]] <- pull_up[[i]][-adj_index[[i]]]
#   }
# }
# print(pull_up2)
# return a list of 4 numeric vectors of pull-up sizes 

# function 1: remove pull-up signals

# function 2: average split peaks

