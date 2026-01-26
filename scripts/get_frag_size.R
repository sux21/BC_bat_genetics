# record allele sizes

source("scripts/helper_functions.R")

# load FSA data. All data will store in a list once the scripts are working.

osiris_out <- read.table("analysis/microsat_data/fragment_analysis_mobix_2024-02-02_plate1_1.tab"
                         , sep = "\t"
                         , header = TRUE
                         , na.strings = "")

min_rfu = 27000
dist_to_stutter = 4
dist_to_pull_up = 0.5
per_sample <- remove_pull_up(osiris_out, "E01_0728-12--M1_004_5759", min_rfu, dist_to_stutter, dist_to_pull_up)
print(per_sample)

# To do: loop through all samples in the data frame 
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





