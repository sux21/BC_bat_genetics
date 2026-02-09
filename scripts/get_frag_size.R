# record allele sizes

source("scripts/allele_calling_helper.R")


# set parameters
ploidy <- 2 # ploidy of organism, number of allelic peaks to be selected
control <- c("negative", "NEG", "Ladder") # negative control and size standard
min_off_rfu <- 20000 # minimum light intensity to find an off-scale peak
ru_len <- 4 # repeat unit length, base pair distance to look for stutters
off_dist <- 0.5 # base pair distance around off-scale peak to find pull-up
A_dist <- 1.5 # base pair distance to find non-template addition
split_dist <- 1 # base pair distance to find off-scale split peaks
peak_ratio <- 0.1 # minimum ratio of light intensity of heterozygous peaks
stutter_ratio <- 0.5 # maximum ratio between stutter and allelic peak
noise_level <- 100 # no extra peaks above shortest allelic peak minus this number  
cont_dist <- 0.5 # base pair distance around peaks in control to be treated as contamination



# load data
path_to_files <- "analysis/microsat_data/"
osiris_files <- list.files(path_to_files, pattern = ".tab")
osiris_paths <- paste0(path_to_files, osiris_files)

osiris_out <- lapply(osiris_paths, load_osiris_tab)
names(osiris_out) <- osiris_files






# remove pull-up 
osiris_out2 <- lapply(osiris_out
                      , remove_pull_up_all
                      , control
                      , min_off_rfu
                      , ru_len
                      , off_dist)

# Note: when allelic peaks overlap, the shorter peak may be treated as
# pull-up and gets removed.



# replace off-scale split peaks with its mean
osiris_out3 <- lapply(osiris_out2
                      , get_mean_split_all
                      , control
                      , min_off_rfu
                      , split_dist)






# select the highest peak in non-template addition
osiris_out4 <- lapply(osiris_out3
                      , get_high_A_all
                      , control
                      , A_dist)








# select allelic peaks
osiris_out5 <- lapply(osiris_out4
                      , allele_caller_all
                      , control
                      , ploidy
                      , peak_ratio
                      , stutter_ratio
                      , noise_level)






# remove fragments found in control (negative control, size standard) 
osiris_out6 <- lapply(osiris_out5
                      , remove_contamination_all
                      , control
                      , cont_dist)

# Note: this removal is by colours but not loci. Contamination in "green"
# channel is only used to remove sample fragments in "green" channel. However,
# if there are two loci are both "green", contamination shown in control of 
# the first locus will be also used to remove fragments in the second locus. 







# transform the data
osiris_out7 <- lapply(osiris_out6
                      , long_to_wide_data
                      , control)



names(osiris_out7)[30]
View(osiris_out7[[30]])
# People should always inspect the data again and correct any errors created by
# this script. 

# 
# get_off_scale2(osiris_out[[30]], "A11_24SP0175--M3_048_6010", min_off_rfu, ru_len)
# 
# find_pull_up(osiris_out[[30]], "A11_24SP0175--M3_048_6010", min_off_rfu
#              , ru_len, off_dist)
# 
# ( get_sample_dat(osiris_out[[30]], "A11_24SP0175--M3_048_6010")$size[[4]] 
#   |> get_adjacent_values(ru_len))


# 
# Articles to read more about:  
#   - Reliable Genotyping of Samples with Very Low DNA Quantities Using PCR 
# - How to track and assess genotyping errors in population genetics studies
# 
# 
# 
# 
# General allele scoring guideline:
#   - Follow previous scoring practice. Fit into their allele class as best as I can. But if I have enough evidence there is new allele, the new allele can be recorded. Sequence the new allele to confirm. 

