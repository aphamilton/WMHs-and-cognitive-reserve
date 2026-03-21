# Project on white matter hyperintensities and cognitive reserve in older adults
# MRI data has already been through several stages of analysis

# Note to self: Collapse all with Alt+O

# ====== Preliminaries ====== ----

# Set working directory ----

setwd("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Data")
output_folder <- "C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Outputs"
# output_folder <- "../Outputs"

# Install and load libraries ----

# Installation instructions for ADNIMERGE if needed in the future
# install.packages("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Libraries/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")

# Load libraries
library(tidyverse)
library(tidyr)
library(janitor)       # Needed for clean_names function
library(Hmisc)         # Needed to run ADNIMERGE
library(ADNIMERGE)
library(data.table)    # To transpose data frame
library(GGally)        # For correlograms

# Separate way of displaying in brain space
library(RNifti)

# For PLS-PM
library(plspm)

# For harmonization
#library(devtools)
#install_github("jfortin1/neuroCombat_Rpackage")
library(neuroCombat)

# For checking near-zero variance
library(caret)

# For PCA (for second-order constructs)
# library(plsdepot)


# ====== Loading and cleaning data ====== ----

# Create functions to load data ----

# Function to load and clean data frame, which is the basis for the more specific load data functions just below
# Use my_file_name for the name of the .csv file
# is_bial_cog is used to indicate some extra steps needed for the cognitive data from the Bialystok data
# Use id_col to indicate the name of the column containing the participant IDs
# Use id_begin and id_end to indicate the range of characters within each file name that correspond to the participant ID
# Use prefix to indicate a string that should be added to the start of each participant ID
# Use exclusions to indicate a vector listing participants to be excluded from the data
# Use dataset_name to specify a string used to identify what dataset a participant is from once datasets are merged later
import_csv_data <- function(my_file_name, is_bial_cog, id_col, id_begin, id_end,
                            prefix, exclusions, dataset_name)
{
  my_df <- read.csv(file = my_file_name, encoding = "UTF-8")
  if (is_bial_cog == TRUE)
    my_df <- my_df[1:99,]
  my_df <- clean_names(my_df)
  # my_df <- my_df %>% rename("ptid" = id_col)   # Did the same thing as the next line, but stopped working when things got updated
  colnames(my_df)[which(names(my_df) == id_col)] <- "ptid"
  my_df$ptid <- substring(my_df$ptid, id_begin, id_end)
  my_df$ptid <- paste(prefix, my_df$ptid, sep="")
  my_df <- my_df[order(my_df$ptid),]
  if (is_bial_cog == TRUE)
      my_df <- filter(my_df, ptid %in% df_bial_gm_lesions$ptid)
  my_df <- filter(my_df, !ptid %in% exclusions)
  my_df <- my_df %>% mutate(dataset = dataset_name) %>% select(ptid, dataset, everything())
  out <- my_df
}

import_adni_mri <- function(my_file_name, id_col)
{
  out <- import_csv_data(my_file_name, FALSE, id_col, 1, 15, "", adni_to_exclude, "ADNI")
}

import_bial_mri <- function(my_file_name, id_col)
{
  out <- import_csv_data(my_file_name, FALSE, id_col, 1, 13, "", bial_to_exclude, "Bialystok")
}

import_spreng_mri <- function(my_file_name, id_col)
{
  my_df <- import_csv_data(my_file_name, FALSE, id_col, 5, 7, "Spreng_", spreng_to_exclude, "Spreng")
  for (i in 1:nrow(my_df))
    if (substr(my_df[i,"ptid"], 10, 10) == "_")
      my_df[i,"ptid"] <- paste(c(substring(my_df[i,"ptid"], 1, 7), "0", substring(my_df[i,"ptid"], 8, 9)), collapse="")
  my_df <- my_df[order(my_df$ptid),]
  out <- my_df
}

import_bial_cog <- function(my_file_name, good_cols)
{
  my_df <- import_csv_data(my_file_name, TRUE, "id", 1, 3, "Bialystok_", bial_to_exclude, "Bialystok")
  my_df <- my_df[,good_cols]   # Note: Remember to add 1 to the end of the range for good_cols to account for the "dataset" column added earlier
  out <- my_df
}


# Set participants to exclude ----

# Participants to exclude based on visual examination of their lesion images
adni_to_exclude   <- c()
bial_to_exclude   <- c()

# Exclusion list created in Jan 2025 - trying to be strict (err on the side of excluding more)
spreng_to_exclude <- c("Spreng_002", "Spreng_007", "Spreng_010", "Spreng_014", "Spreng_052",
                       "Spreng_053", "Spreng_054", "Spreng_055", "Spreng_056", "Spreng_057",
                       "Spreng_058", "Spreng_061", "Spreng_062", "Spreng_063", "Spreng_064",
                       "Spreng_066", "Spreng_067", "Spreng_068", "Spreng_069", "Spreng_070",
                       "Spreng_071", "Spreng_072", "Spreng_104", "Spreng_105", "Spreng_107",
                       "Spreng_108", "Spreng_109", "Spreng_110", "Spreng_111", "Spreng_112",
                       "Spreng_114", "Spreng_117", "Spreng_118", "Spreng_119", "Spreng_120",
                       "Spreng_122", "Spreng_124", "Spreng_125", "Spreng_184", "Spreng_185",
                       "Spreng_186", "Spreng_187", "Spreng_188", "Spreng_190", "Spreng_196",
                       "Spreng_197", "Spreng_198", "Spreng_200", "Spreng_201", "Spreng_202",
                       "Spreng_203", "Spreng_204", "Spreng_205", "Spreng_272", "Spreng_289",
                       "Spreng_290", "Spreng_293", "Spreng_297")


# Load MRI analysis data ----

# Load data frames containing MRI analysis outputs for ADNI dataset
df_adni_gm_lesions   <- import_adni_mri("ADNI_gray_matter_lesions_20241117.csv", "filename")
df_adni_wm_discon    <- import_adni_mri("ADNI_percent_discon_tracts_20241117.csv", "filename")
df_adni_wm_lesions   <- import_adni_mri("ADNI_lesion_volumes_numbers_20231118.csv", "participant")
df_adni_wm_lesions   <- filter(df_adni_wm_lesions, ptid %in% df_adni_gm_lesions$ptid)

# Load data frames containing MRI analysis outputs for Bialystok dataset
df_bial_gm_lesions   <- import_bial_mri("Bialystok_gray_matter_lesions_20241124.csv", "filename")
df_bial_wm_discon    <- import_bial_mri("Bialystok_percent_discon_tracts_20241124.csv", "filename")
df_bial_wm_lesions   <- import_bial_mri("Bialystok_lesion_volumes_numbers_20231118.csv", "participant")
df_bial_wm_lesions   <- filter(df_bial_wm_lesions, ptid %in% df_bial_gm_lesions$ptid)

# Load data frames containing MRI analysis outputs for Spreng dataset
df_spreng_gm_lesions <- import_spreng_mri("Spreng_gray_matter_lesions_20241124.csv", "filename")
df_spreng_wm_discon  <- import_spreng_mri("Spreng_percent_discon_tracts_20241124.csv", "filename")

# Doing the WM lesion summaries for the Spreng dataset separately because the IDs are in a bad format
df_spreng_wm_lesions <- read.csv(file = "Spreng_lesion_volumes_numbers_20231118.csv", encoding = "UTF-8")
df_spreng_wm_lesions <- clean_names(df_spreng_wm_lesions)
colnames(df_spreng_wm_lesions)[which(names(df_spreng_wm_lesions) == "participant")] <- "ptid"
df_spreng_wm_lesions$ptid <- paste("Spreng_", df_spreng_wm_lesions$ptid, sep="")
df_spreng_wm_lesions <- df_spreng_wm_lesions[order(df_spreng_wm_lesions$ptid),]
df_spreng_wm_lesions <- df_spreng_wm_lesions %>% mutate(dataset = "Spreng") %>% select(ptid, dataset, everything())
for (i in 1:nrow(df_spreng_wm_lesions))
{
  if (nchar(df_spreng_wm_lesions[i,"ptid"]) == 8)
    df_spreng_wm_lesions[i,"ptid"] <- paste(c(substring(df_spreng_wm_lesions[i,"ptid"], 1, 7), "00", substring(df_spreng_wm_lesions[i,"ptid"], 8, 8)), collapse="")
  if (nchar(df_spreng_wm_lesions[i,"ptid"]) == 9)
    df_spreng_wm_lesions[i,"ptid"] <- paste(c(substring(df_spreng_wm_lesions[i,"ptid"], 1, 7), "0", substring(df_spreng_wm_lesions[i,"ptid"], 8, 9)), collapse="")
}
df_spreng_wm_lesions <- df_spreng_wm_lesions[order(df_spreng_wm_lesions$ptid),]
df_spreng_wm_lesions <- filter(df_spreng_wm_lesions, !ptid %in% spreng_to_exclude)
df_spreng_wm_lesions <- filter(df_spreng_wm_lesions, ptid %in% df_spreng_gm_lesions$ptid)

# Import parcel names
parcel_names_wm <- read.csv(file = "Parcel_names_WM.csv", encoding = "UTF-8", header = FALSE)[,1]
parcel_names_gm <- read.csv(file = "Parcel_names_GM.csv", encoding = "UTF-8", header = FALSE)[,1]
# Add parcel names to data frames (and clean the names)
names(df_adni_gm_lesions)[c(3:137)] <- parcel_names_gm
df_adni_gm_lesions <- clean_names(df_adni_gm_lesions)
names(df_bial_gm_lesions)[c(3:137)] <- parcel_names_gm
df_bial_gm_lesions <- clean_names(df_bial_gm_lesions)
names(df_spreng_gm_lesions)[c(3:137)] <- parcel_names_gm
df_spreng_gm_lesions <- clean_names(df_spreng_gm_lesions)
names(df_adni_wm_discon)[c(3:72)] <- parcel_names_wm
df_adni_wm_discon <- clean_names(df_adni_wm_discon)
names(df_bial_wm_discon)[c(3:72)] <- parcel_names_wm
df_bial_wm_discon <- clean_names(df_bial_wm_discon)
names(df_spreng_wm_discon)[c(3:72)] <- parcel_names_wm
df_spreng_wm_discon <- clean_names(df_spreng_wm_discon)


# Load non-MRI data ----

# Load main adnimerge data frame
data(adnimerge)
df_adni_adnimerge <- adnimerge
rm(adnimerge)
df_adni_adnimerge <- clean_names(df_adni_adnimerge)
df_adni_adnimerge <- df_adni_adnimerge[df_adni_adnimerge$colprot == "ADNI3", ]
df_adni_adnimerge$ptid <- paste("ADNI_", df_adni_adnimerge$ptid, sep="")
df_adni_adnimerge <- filter(df_adni_adnimerge, ptid %in% df_adni_gm_lesions$ptid)
df_adni_adnimerge <- df_adni_adnimerge %>%
  group_by(ptid) %>%
  filter(row_number() == 2)
df_adni_adnimerge <- df_adni_adnimerge[order(df_adni_adnimerge$ptid),]
df_adni_adnimerge <- df_adni_adnimerge %>% mutate(dataset = "ADNI") %>% select(ptid, dataset, everything())

# Load uwnpsychsum data frame via adnimerge
# This one has aggregate scores for cognitive tasks in overall domains for the ADNI data
data(uwnpsychsum)
df_adni_uwnpsychsum <- uwnpsychsum
rm(uwnpsychsum)
df_adni_uwnpsychsum <- clean_names(df_adni_uwnpsychsum)
df_adni_uwnpsychsum <- df_adni_uwnpsychsum[df_adni_uwnpsychsum$colprot == "ADNI3", ]

# Joining adnimerge and uwnpsychsum, within df_adni_uwnpsychsum
df_adni_adnimerge$rid <- as.numeric(df_adni_adnimerge$rid)
df_adni_adnimerge$viscode <- as.character(df_adni_adnimerge$viscode)
df_adni_uwnpsychsum$viscode <- as.character(df_adni_uwnpsychsum$viscode)
df_adni_uwnpsychsum <- inner_join(df_adni_uwnpsychsum, df_adni_adnimerge, by = c("rid", "viscode"))
df_adni_uwnpsychsum <- df_adni_uwnpsychsum[order(df_adni_uwnpsychsum$rid),]
df_adni_uwnpsychsum <- df_adni_uwnpsychsum %>% mutate(dataset = "ADNI") %>% select(rid, dataset, everything())

# Loading neuropsychiatric symptoms
data(npi)
df_adni_npi <- npi
rm(npi)
df_adni_npi <- clean_names(df_adni_npi)
df_adni_npi <- df_adni_npi[df_adni_npi$colprot == "ADNI3", ]
df_adni_npi$viscode <- as.character(df_adni_npi$viscode)
df_adni_npi <- inner_join(df_adni_npi, df_adni_adnimerge, by = c("rid", "viscode"))
df_adni_npi <- df_adni_npi[order(df_adni_npi$rid),]
df_adni_npi <- df_adni_npi %>% mutate(dataset = "ADNI") %>% select(rid, dataset, everything())

# Loading and cleaning behavioral data from sheets from the Bialystok .xls file
df_bial_mmse <- import_bial_cog("Bialystok_cog_data_MMSE_20231113.csv",     1:12)
df_bial_lsbq <- import_bial_cog("Bialystok_cog_data_LSBQ_20231113.csv",     1:129)
df_bial_criq <- import_bial_cog("Bialystok_cog_data_CRIQ_20231113.csv",     1:59)
df_bial_ship <- import_bial_cog("Bialystok_cog_data_Shipley2_20231113.csv", 1:11)
df_bial_tmt  <- import_bial_cog("Bialystok_cog_data_TMT_20231113.csv",      1:27)
df_bial_vft  <- import_bial_cog("Bialystok_cog_data_VFT_20231113.csv",      1:57)
df_bial_cwit <- import_bial_cog("Bialystok_cog_data_CWIT_20231113.csv",     1:45)
df_bial_lsbq_factors <- import_bial_cog("Bialystok_cog_data_LSBQ_factors_20240313.csv", 1:5)

# Loading and cleaning cognitive/demographic/etc data from Spreng dataset
df_spreng_cog <- read.csv(file = "Spreng_cog_data_20231016.csv", encoding = "UTF-8")
df_spreng_cog <- clean_names(df_spreng_cog)
colnames(df_spreng_cog)[which(names(df_spreng_cog) == "id")] <- "ptid"
for (i in 1:99)
{
  df_spreng_cog[i,"ptid"] <- paste(c(substring(df_spreng_cog[i,"ptid"], 1, 4), "0", substring(df_spreng_cog[i,"ptid"], 5, 6)), collapse="")
}
df_spreng_cog$ptid <- substring(df_spreng_cog$ptid, 5, 7)
df_spreng_cog$ptid <- paste("Spreng_", df_spreng_cog$ptid, sep="")
df_spreng_cog <- filter(df_spreng_cog, ptid %in% df_spreng_gm_lesions$ptid)
df_spreng_cog <- df_spreng_cog[order(df_spreng_cog$ptid),]
df_spreng_cog <- filter(df_spreng_cog, !ptid %in% spreng_to_exclude)
df_spreng_cog <- df_spreng_cog %>% mutate(dataset = "Spreng") %>% select(ptid, dataset, everything())


# Remove participants who don't meet inclusion criteria  ----
# Shouldn't be necessary since only the good participants went through MRI processing, but can leave in just in case

# "Young" group from the Spreng dataset
df_spreng_cog        <- df_spreng_cog[df_spreng_cog$agegroup != "Y",]
df_spreng_wm_lesions <- filter(df_spreng_wm_lesions, ptid %in% df_spreng_cog$ptid)
df_spreng_gm_lesions <- filter(df_spreng_gm_lesions, ptid %in% df_spreng_cog$ptid)
df_spreng_wm_discon  <- filter(df_spreng_wm_discon,  ptid %in% df_spreng_cog$ptid)

# ADNI participants who don't have dx marked as CN
df_adni_adnimerge   <- df_adni_adnimerge[df_adni_adnimerge$dx == "CN",]
df_adni_adnimerge   <- df_adni_adnimerge[!is.na(df_adni_adnimerge$ptid),]    # Extra check just in case
df_adni_uwnpsychsum <- filter(df_adni_uwnpsychsum, ptid %in% df_adni_adnimerge$ptid)
df_adni_wm_lesions  <- filter(df_adni_wm_lesions,  ptid %in% df_adni_adnimerge$ptid)
df_adni_wm_discon   <- filter(df_adni_wm_discon,   ptid %in% df_adni_adnimerge$ptid)
df_adni_gm_lesions  <- filter(df_adni_gm_lesions,  ptid %in% df_adni_adnimerge$ptid)


# Combining data frames ----

# Creating combined data frame for non-MRI measures available for all three datasets
# Note: For sex, 0 = female and 1 = male; datasets differed on marking sex or gender but only binary values were given

df_adni_cog_to_combine   <- data.frame("ptid" = df_adni_wm_lesions$ptid, "dataset" = "ADNI",
                                       "dx" = df_adni_adnimerge$dx, "sex" = df_adni_adnimerge$ptgender,
                                       "age" = df_adni_adnimerge$age, "yrs_education" = df_adni_adnimerge$pteducat,
                                       "mmse" = df_adni_adnimerge$mmse)
df_adni_cog_to_combine[df_adni_cog_to_combine$sex == "Female",]$sex <- 0
df_adni_cog_to_combine[df_adni_cog_to_combine$sex == "Male",]$sex   <- 1

df_bial_cog_to_combine   <- data.frame("ptid" = df_bial_wm_lesions$ptid, "dataset" = "Bialystok",
                                       "dx" = "CN", "sex" = df_bial_lsbq$sex,
                                       "age" = df_bial_lsbq$total_age, "yrs_education" = df_bial_criq$yrs_of_edu,
                                       "mmse" = df_bial_mmse$mmse_total_score)

df_spreng_cog_to_combine <- data.frame("ptid" = df_spreng_wm_lesions$ptid, "dataset" = "Spreng",
                                       "dx" = "CN", "sex" = df_spreng_cog$gender,
                                       "age" = df_spreng_cog$age, "yrs_education" = df_spreng_cog$education,
                                       "mmse" = df_spreng_cog$mmse)
df_spreng_cog_to_combine[df_spreng_cog_to_combine$sex == "F",]$sex <- 0
df_spreng_cog_to_combine[df_spreng_cog_to_combine$sex == "M",]$sex <- 1

df_combined_cog <- rbind(df_adni_cog_to_combine, df_bial_cog_to_combine, df_spreng_cog_to_combine)
# rm(df_adni_cog_to_combine)
# rm(df_bial_cog_to_combine)
# rm(df_spreng_cog_to_combine)

# One combined data frame for each MRI measure type
df_combined_wm_lesions <- left_join(df_combined_cog,
                                    rbind(df_adni_wm_lesions,
                                          df_bial_wm_lesions,
                                          df_spreng_wm_lesions),
                                    by = c("ptid", "dataset"))
df_combined_gm_lesions <- left_join(df_combined_cog,
                                    rbind(df_adni_gm_lesions,
                                          df_bial_gm_lesions,
                                          df_spreng_gm_lesions),
                                    by = c("ptid", "dataset"))
df_combined_wm_discon  <- left_join(df_combined_cog,
                                    rbind(df_adni_wm_discon,
                                          df_bial_wm_discon,
                                          df_spreng_wm_discon),
                                    by = c("ptid", "dataset"))

# One data frame with all participants showing means of all MRI measure types
df_combined_mri_means <- df_combined_wm_lesions
df_combined_mri_means <- mutate(df_combined_mri_means,
                                gm_lesions_mean = rowMeans(select(df_combined_gm_lesions, 8:142),
                                                           na.rm = TRUE))
df_combined_mri_means <- mutate(df_combined_mri_means,
                                wm_discon_mean = rowMeans(select(df_combined_wm_discon, 8:77),
                                                          na.rm = TRUE))


# ====== Creating dfs for PLS-PM ====== ----

# Create reduced df, from which I will create the PLS-PM dfs ----

combined_df_reduced <- merge(df_combined_cog[,c(1:2,4:7)],
                             df_combined_wm_discon[,c(1,8:77)],
                             by="ptid")
combined_df_reduced <- merge(combined_df_reduced,
                             df_combined_gm_lesions[,c(1,8:142)],
                             by="ptid")
combined_df_reduced <- merge(combined_df_reduced,
                             df_combined_wm_lesions[,c(1,8:9)],
                             by="ptid")
combined_df_reduced$sex <- as.numeric(combined_df_reduced$sex)

combined_df_reduced <- combined_df_reduced[complete.cases(combined_df_reduced[,c("sex","age","yrs_education","mmse")]),]

# # Creating copy for later for calculating means, before any cols removed pre-harmonization
# combined_df_reduced_all_parcels <- combined_df_reduced


# Conduct harmonization for the combined df ----

# List to keep all components used for harmonization together
harmonization_specs <- list()

# Preparing inputs for the harmonization function
harmonization_specs[["combined_df_trans"]] <- transpose(
      combined_df_reduced[,c(7:213)],
      keep.names = NULL)
harmonization_specs[["combined_batch_vector"]] <- combined_df_reduced$dataset
harmonization_specs[["combined_batch_vector"]][harmonization_specs[["combined_batch_vector"]] == "ADNI"]      <- 1
harmonization_specs[["combined_batch_vector"]][harmonization_specs[["combined_batch_vector"]] == "Bialystok"] <- 2
harmonization_specs[["combined_batch_vector"]][harmonization_specs[["combined_batch_vector"]] == "Spreng"]    <- 3
harmonization_specs[["combined_batch_vector"]] <- as.numeric(harmonization_specs[["combined_batch_vector"]])
harmonization_specs[["combined_mod_vector"]]   <- model.matrix(~ combined_df_reduced$sex
                                                               )

# Running the harmonization itself
harmonization_specs[["combined_harmonized_data"]] <- neuroCombat(dat=harmonization_specs[["combined_df_trans"]],
                                                                 batch=harmonization_specs[["combined_batch_vector"]],
                                                                 mod=harmonization_specs[["combined_mod_vector"]]
                                                                 )

# Turning harmonization output into df
harmonization_specs[["combined_harmonized_df"]] <- as.data.frame(harmonization_specs[["combined_harmonized_data"]][[1]])

# Transpose the data frame back
combined_df_post_harm <- transpose(harmonization_specs[["combined_harmonized_df"]],
                                   keep.names = NULL)

# Re-inserting non-MRI variables
combined_df_post_harm$dataset       <- combined_df_reduced$dataset
combined_df_post_harm$age           <- combined_df_reduced$age
combined_df_post_harm$yrs_education <- combined_df_reduced$yrs_education
combined_df_post_harm$mmse          <- combined_df_reduced$mmse

combined_df_post_harm <- combined_df_post_harm %>%
                         select(c("dataset","age","yrs_education","mmse"), everything())

# Adding column titles to harmonized df
names(combined_df_post_harm)[c(5:74)]    <- parcel_names_wm
names(combined_df_post_harm)[c(75:209)]  <- parcel_names_gm[c(1:135)]
names(combined_df_post_harm)[c(210:211)] <- c("lesion_volume_in_m_l","number_of_lesions")
combined_df_post_harm <- clean_names(combined_df_post_harm)


# Functions to remove near-zero var cols and create list with all specs to run PLS-PM ----

# Create a df with near-zero var cols removed, and store the various components in a list
# The other five parameters are for the set of column indices that correspond to each type of measure
prepare_plspm_df_info <- function(my_df, res_cols, cog_cols, wm_cols, gm_cols, wb_cols, my_freq_cut)
{
  # Create the list which will serve as the overall data structure
  my_df_info <- vector("list", length = 6)
  names(my_df_info) <- c("original_df", "original_cols",
                         "nzv_cols", "reduced_df", "reduced_col_nums", "reduced_cols")
  
  # Vector contain the names for the column types (to be used within the function, not as an output)
  col_type_names <- c("res_cols", "cog_cols", "wm_cols", 
                      "gm_cort_cols", "gm_subc_cols", "wb_cols")
  
  # Add the original df
  my_df_info[["original_df"]]  <- my_df
  
  # Add the original column ranges
  my_df_info[["original_cols"]] <- vector("list", length = 6)
  names(my_df_info[["original_cols"]]) <- col_type_names
  my_df_info[["original_cols"]][["res_cols"]]     <- res_cols
  my_df_info[["original_cols"]][["cog_cols"]]     <- cog_cols
  my_df_info[["original_cols"]][["wm_cols"]]      <- wm_cols
  my_df_info[["original_cols"]][["gm_cort_cols"]] <- gm_cols[1:100]
  my_df_info[["original_cols"]][["gm_subc_cols"]] <- gm_cols[101:135]
  my_df_info[["original_cols"]][["wb_cols"]]      <- wb_cols
  
  # Find columns with near-zero variance
  my_df_info[["nzv_cols"]] <- vector("list", length = 6)
  names(my_df_info[["nzv_cols"]]) <- col_type_names
  for (col_type in 1:6)
  {
    my_df_info[["nzv_cols"]][[col_type]] <- nearZeroVar(
            my_df[,my_df_info[["original_cols"]][[col_type]]],
            freqCut = my_freq_cut,
            uniqueCut = 20) + (my_df_info[["original_cols"]][[col_type]][1] - 1)
  }
  
  # Check whether there any cols with near-zero variance
  no_nzv_cols <- TRUE
  for (col_type in 1:6)
  {
    if (length(my_df_info[["nzv_cols"]][[col_type]]) != 0)
      no_nzv_cols <- FALSE
  }
  
  # Remove columns with near-zero variance
  my_df_info[["reduced_df"]] <- my_df_info[["original_df"]]
  # The if statement is due to an error that occurs when trying to remove 0 cols from a df (not sure why)
  if (no_nzv_cols == FALSE)
  {
    my_df_info[["reduced_df"]] <- my_df_info[["reduced_df"]][,-c(my_df_info[["nzv_cols"]][[1]],
                                                                 my_df_info[["nzv_cols"]][[2]],
                                                                 my_df_info[["nzv_cols"]][[3]],
                                                                 my_df_info[["nzv_cols"]][[4]],
                                                                 my_df_info[["nzv_cols"]][[5]],
                                                                 my_df_info[["nzv_cols"]][[6]])]
  }
  
  # Number of columns by measure type after near-zero variance columns removed
  # Includes an extra column before 5 from above, for columns that aren't in those groups
  my_df_info[["reduced_col_nums"]] <- vector("list", length = 7)
  names(my_df_info[["reduced_col_nums"]]) <- c("early_cols", col_type_names)
  my_df_info[["reduced_col_nums"]][["early_cols"]] <- my_df_info[["original_cols"]][["res_cols"]][1] - 1
  for (col_type in 1:6)
  {
    my_df_info[["reduced_col_nums"]][[col_type_names[col_type]]] <- length(my_df_info[["original_cols"]][[col_type_names[col_type]]]) - length(my_df_info[["nzv_cols"]][[col_type_names[col_type]]])
  }
  
  # Cumulative number of columns by measure type
  cumulative_col_nums <- vector("list", length = 7)
  names(cumulative_col_nums) <- c("early_cols", col_type_names)
  cumulative_col_nums[1] <- my_df_info[["reduced_col_nums"]][1]
  for (i in 2:7)
  {
    cumulative_col_nums[[i]] <- cumulative_col_nums[[i-1]] + as.numeric(my_df_info[["reduced_col_nums"]][i])
  }
  
  # Create reduced column ranges
  my_df_info[["reduced_cols"]] <- vector("list", length = 6)
  names(my_df_info[["reduced_cols"]]) <- c(col_type_names)
  for (i in 1:6)
  {
    my_df_info[["reduced_cols"]][[i]] <- c((cumulative_col_nums[[i]] + 1):
                                           (cumulative_col_nums[[i+1]]))
  }
  
  return(my_df_info)
}

# Now create a list to contain the lists of PLS-PM preparation info for all dfs
plspm_prep_all_dfs <- list()


# Create the dfs for models from the combined df and its subsets ----
# This section creates dfs for combined dataset and for individual datasets
# But only includes variables shared between the three datasets
# The versions with additional variables are ADNI_only, Bial_only, and Spreng_only

# Combined df (using the data after harmonization)
df_plspm_combined <- combined_df_post_harm

# Creating copy of the reduced df (before harmonization)
combined_df_reduced_no_harm <- combined_df_reduced

# Creating reduced data frames for the individual datasets
df_plspm_adni   <- combined_df_reduced_no_harm[combined_df_reduced_no_harm$dataset == "ADNI",] %>% select(-c(ptid, sex))
df_plspm_bial   <- combined_df_reduced_no_harm[combined_df_reduced_no_harm$dataset == "Bialystok",] %>% select(-c(ptid, sex))
df_plspm_spreng <- combined_df_reduced_no_harm[combined_df_reduced_no_harm$dataset == "Spreng",] %>% select(-c(ptid, sex))

# Process those dfs and add the results to the top-level lists
plspm_prep_all_dfs[["Combined"]] <- prepare_plspm_df_info(
       my_df = df_plspm_combined,
       res_cols = 3, cog_cols = 4,
       wm_cols = 5:74, gm_cols = 75:209, wb_cols = 210:211,
       my_freq_cut = 95/5)
plspm_prep_all_dfs[["ADNI"]] <- prepare_plspm_df_info(
       my_df = df_plspm_adni,
       res_cols = 3, cog_cols = 4,
       wm_cols = 5:74, gm_cols = 75:209, wb_cols = 210:211,
       my_freq_cut = 95/5)
plspm_prep_all_dfs[["Bialystok"]] <- prepare_plspm_df_info( #_2(
       my_df = df_plspm_bial,
       res_cols = 3, cog_cols = 4,
       wm_cols = 5:74, gm_cols = 75:209, wb_cols = 210:211,
       my_freq_cut = 85/15)   # More stringent exclusion of near-zero var cols b/c PLS-PM could not converge
plspm_prep_all_dfs[["Spreng"]] <- prepare_plspm_df_info( #_2(
       my_df = df_plspm_spreng,
       res_cols = 3, cog_cols = 4,
       wm_cols = 5:74, gm_cols = 75:209, wb_cols = 210:211,
       my_freq_cut = 85/15)   # More stringent exclusion of near-zero var cols b/c PLS-PM could not converge


# Create df for separate ADNI model ----

# Merge relevant behavioral and MRI data for ADNI dataset
df_plspm_adni_only <- merge(df_adni_uwnpsychsum[,c("ptid","age")],
                            df_adni_npi[,c("ptid","npitotal")],
                            by="ptid")
df_plspm_adni_only <- merge(df_plspm_adni_only,
                            df_adni_uwnpsychsum[,c("ptid","pteducat",
                                                   "adni_mem","adni_ef2","adni_lan","adni_vs")],
                            by="ptid")
df_plspm_adni_only <- merge(df_plspm_adni_only,
                            df_adni_wm_discon[,c(1,3:72)],
                            by="ptid")
df_plspm_adni_only <- merge(df_plspm_adni_only,
                            df_adni_gm_lesions[,c(1,3:137)],
                            by="ptid")
df_plspm_adni_only <- merge(df_plspm_adni_only,
                            df_adni_wm_lesions[,c(1,3:4)],
                            by="ptid")

# Other df cleaning
colnames(df_plspm_adni_only)[which(names(df_plspm_adni_only) == "pteducat")] <- "yrs_education"
df_plspm_adni_only <- df_plspm_adni_only[complete.cases(df_plspm_adni_only[,2:8]),]
df_plspm_adni_only$adni_mem <- as.numeric(df_plspm_adni_only$adni_mem)
df_plspm_adni_only$adni_ef2 <- as.numeric(df_plspm_adni_only$adni_ef2)
df_plspm_adni_only$adni_lan <- as.numeric(df_plspm_adni_only$adni_lan)
df_plspm_adni_only$adni_vs  <- as.numeric(df_plspm_adni_only$adni_vs)

# Process the df and add the results to the top-level list
plspm_prep_all_dfs[["ADNI_only"]] <- prepare_plspm_df_info(
       my_df = df_plspm_adni_only,
       res_cols = 4, cog_cols = 5:8,
       wm_cols = 9:78, gm_cols = 79:213, wb_cols = 214:215,
       my_freq_cut = 95/5)


# Create df for separate Bialystok model ----

# Merge relevant behavioral data for bial dataset
# Note that the LSBQ data didn't end up being used in the analyses
df_plspm_bial_only <- merge(df_bial_cog_to_combine[,c("ptid","age","sex","yrs_education")],
                            df_bial_lsbq[,c("ptid","group")],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_lsbq_factors[,c("ptid",
                                                    "use","proficiency","english")],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_criq[,c("ptid",
                                            "cri_education",
                                            "cri_working_activity",
                                            "cri_leisure_time")],
                            by="ptid")

# Alternate version to be implemented, just uses scaled scores
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_ship[,c("ptid",
                                            "shipley_composite_b_standard")],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_tmt[,c("ptid",
                                           "tmt_scaled_score_difference")],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_vft[,c("ptid",
                                           "vft_c3_category_switching_total_correct_scaled_score")],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_cwit[,c("ptid",
                                            "cwit_c3_inhibition_scaled_score")],
                            by="ptid")

# Renaming columns to be simpler (and specify the battery they are from where applicable)
# Bilingual status (1=monolingual, 2=bilingual)
colnames(df_plspm_bial_only)[which(names(df_plspm_bial_only) == "group")] <- "bilingual_status"
df_plspm_bial_only$bilingual_status <- as.factor(df_plspm_bial_only$bilingual_status)
# Cognitive reserve
cri_col_names <- c("cri_edu", "cri_work", "cri_leisure")
colnames(df_plspm_bial_only)[9:11] <- cri_col_names
rm(cri_col_names)
# Cognitive tasks (scaled scores)
cog_task_names <- c("shipley_scaled","tmt_scaled","vft_scaled","cwit_scaled")
colnames(df_plspm_bial_only)[12:15] <- cog_task_names
rm(cog_task_names)

# Set columns to numeric
for (my_col in 6:15)
{
  df_plspm_bial_only[,my_col] <- as.numeric(df_plspm_bial_only[,my_col])
}

# Remove rows with missing values
df_plspm_bial_only <- df_plspm_bial_only[complete.cases(df_plspm_bial_only[,c(2,9:15)]),]

# Merge relevant MRI data for bial dataset
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_wm_discon[,c(1,3:72)],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_gm_lesions[,c(1,3:137)],
                            by="ptid")
df_plspm_bial_only <- merge(df_plspm_bial_only,
                            df_bial_wm_lesions[,c(1,3:4)],
                            by="ptid")

# Process the df and add the results to the top-level list
plspm_prep_all_dfs[["Bialystok_only"]] <- prepare_plspm_df_info( #_2(
       my_df = df_plspm_bial_only,
       res_cols = 9:11, cog_cols = 12:15,
       wm_cols = 16:85, gm_cols = 86:220, wb_cols = 221:222,
       my_freq_cut = 85/15)   # More stringent exclusion of near-zero var cols b/c PLS-PM could not converge


# Create df for separate Spreng model ----

# Merge relevant behavioral and MRI data for spreng dataset
df_plspm_spreng_only <- merge(df_spreng_cog_to_combine[,c("ptid","age","yrs_education")],
                              df_spreng_cog[,c("ptid",
                                               "episodic_index",
                                               "semantic_index",
                                               "executive_index")],
                              by="ptid")
df_plspm_spreng_only <- merge(df_plspm_spreng_only,
                              df_spreng_wm_discon[,c(1,3:72)],
                              by="ptid")
df_plspm_spreng_only <- merge(df_plspm_spreng_only,
                              df_spreng_gm_lesions[,c(1,3:137)],
                              by="ptid")
df_plspm_spreng_only <- merge(df_plspm_spreng_only,
                              df_spreng_wm_lesions[,c(1,3:4)],
                              by="ptid")

# Remove rows with missing data
df_plspm_spreng_only <- df_plspm_spreng_only[complete.cases(df_plspm_spreng_only[,2:6]),]

# Process the df and add the results to the top-level list
plspm_prep_all_dfs[["Spreng_only"]] <- prepare_plspm_df_info( #_2(
       my_df = df_plspm_spreng_only,
       res_cols = 3, cog_cols = 4:6,
       wm_cols = 7:76, gm_cols = 77:211, wb_cols = 212:213,
       my_freq_cut = 85/15)   # More stringent exclusion of near-zero var cols b/c PLS-PM could not converge


# ====== Making PLS-PM models ====== ----

# Functions for running PLS-PM ----

# Function to insert the correct cols by measure type for this dataset
select_measure_cols <- function(x, my_measure_cols)
{
  ifelse (x == "res_cols",
          return(my_measure_cols[["res_cols"]]),
          x)
  ifelse (x == "cog_cols",
          return(my_measure_cols[["cog_cols"]]),
          x)
  ifelse (x == "wm_cols",
          return(my_measure_cols[["wm_cols"]]),
          x)
  ifelse (x == "gm_cols",
          return(c(my_measure_cols[["gm_cort_cols"]],
                   my_measure_cols[["gm_subc_cols"]])),
          x)
  ifelse (x == "gm_cort_cols",
          return(my_measure_cols[["gm_cort_cols"]]),
          x)
  ifelse (x == "gm_subc_cols",
          return(my_measure_cols[["gm_subc_cols"]]),
          x)
  ifelse (x == "wb_cols",
          return(my_measure_cols[["wb_cols"]]),
          x)
  ifelse (x == "parcel_cols",
          return(c(my_measure_cols[["wm_cols"]],
                   my_measure_cols[["gm_cort_cols"]],
                   my_measure_cols[["gm_subc_cols"]])),
          x)
  ifelse (x == "brain_cols",
          return(c(my_measure_cols[["wm_cols"]],
                   my_measure_cols[["gm_cort_cols"]],
                   my_measure_cols[["gm_subc_cols"]],
                   my_measure_cols[["wb_cols"]])),
          x)
}

# Create a list w/ the configuration for a PLS-PM model
# And then a model for each dataset
create_plspms <- function(my_variables,    my_matrix,
                          my_blocks,       my_modes,
                          my_datasets)
{
  # Create the list which will serve as the overall data structure
  my_plspms <- vector("list", length = 6)
  names(my_plspms) <- c("Variables", "Matrix", "Blocks",
                        "Modes", "Datasets", "Models")
  
  # List the variables included in this PLS-PM configuration
  my_plspms[["Variables"]] <- my_variables
  
  # Create the path matrix for this PLS-PM configuration
  my_plspms[["Matrix"]] <- matrix(my_matrix,
                                  nrow=length(my_variables),
                                  ncol=length(my_variables),
                                  byrow=TRUE)
  colnames(my_plspms[["Matrix"]]) <- my_variables
  rownames(my_plspms[["Matrix"]]) <- my_variables
  
  # List the data columns that correspond to each variable in the model
  my_plspms[["Blocks"]] <- my_blocks
  
  # List the mode to be used for each variable in the model
  my_plspms[["Modes"]] <- my_modes
  
  # List the datasets for which models will be created
  my_plspms[["Datasets"]] <- my_datasets
  
  # Create a PLS-PM model with each dataset using the above configuration
  my_plspm_models <- vector("list", length = length(my_datasets))
  names(my_plspm_models) <- my_datasets
  for (my_dataset in 1:(length(my_datasets)))
  {
    my_temp_blocks <- lapply(my_plspms[["Blocks"]],
                             select_measure_cols,
                             plspm_prep_all_dfs[[my_datasets[my_dataset]]][["reduced_cols"]])
    my_plspm_models[[my_dataset]] <- plspm(plspm_prep_all_dfs[[my_datasets[my_dataset]]][["reduced_df"]],
                                           my_plspms[["Matrix"]],
                                           my_temp_blocks,
                                           modes=my_plspms[["Modes"]],
                                           boot.val = TRUE, br = 500)
  }
  my_plspms[["Models"]] <- my_plspm_models
  
  # Return the completed data structure containing the configuration and the models
  return(my_plspms)
}

# Check the properties of a PLS-PM model and plot the inner model and loadings
analyze_plspm_model <- function(plspm_model)
{
  
  # General summary of characteristics
  summary(plspm_model)
  
  # Check measures of unidimensionality
  plspm_model$unidim
  # plot(plspm_model, what = "loadings")
  
  # Show the outer model
  plspm_model$outer_model
  
  # Check cross loadings for "traitors"
  plspm_model$crossloadings
  
  # # Check the overall coefficients between the factors
  # plot(plspm_model)
  # Matrix of path coefficients (same information but not visually)
  plspm_model$path_coefs
  
  # Check the inner model
  plspm_model$inner_model
  
  # Need to review what this step does
  plspm_model$effects
  
  # Need to review what this step does
  plspm_model$inner_summary
  
  # Goodness of fit
  plspm_model$gof
  
  # Bootstrap results
  plspm_model$boot
  
  # Return the ggplot or else it won't print (even without <- )
  # It also worked to have the ggplot as the last thing in the function
  # A forum suggested using the print command around the ggplot but that didn't help
  # return(loadings_barchart)

}

# Create a barplot of the loadings of barcharts
plot_plspm_loadings <- function(plspm_model)
{
  loadings_barchart <- ggplot(data = plspm_model$outer_model, aes(x = name, y = loading, fill = block)) +
                       geom_bar(stat = "identity", position = "dodge") +
                       theme_bw() + theme(plot.title=element_text(hjust=0.5)) +
                       theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank()) +
                       ggtitle("Barchart of Loadings") + xlab("Observed Variable") + ylab("Loading")
  return(loadings_barchart)
}


# Create the models ----

# Create top-level list for PLS-PM models
plspms_all_datasets <- list()
# Create list of all datasets
all_dataset_names <- c("Combined",
                       "ADNI", "Bialystok", "Spreng",
                       "ADNI_only", "Bialystok_only", "Spreng_only")

# Age, cog, GM disruption, cog res
plspms_all_datasets[["cog_gmp_age_res"]] <- create_plspms(
      my_variables=c("Age",        "Cog_Res",
                     "Lesions_GM", "Cognition"),
      my_matrix   =c(0,0,0,0,    # age
                     0,0,0,0,    # cog_reserve
                     1,0,0,0,    # lesions_gm_parcels
                     1,1,1,0),   # cog_ability
      my_blocks   =list(2,             # age
                        "res_cols",    # cog_reserve
                        "gm_cols",     # lesions_gm_parcels
                        "cog_cols"),   # cog_ability
      my_modes    =rep("A",4),
      my_datasets =all_dataset_names)

# Age, cog, WM disruption, cog res
plspms_all_datasets[["cog_wmp_age_res"]] <- create_plspms(
      my_variables=c("Age",        "Cog_Res",
                     "Lesions_WM", "Cognition"),
      my_matrix   =c(0,0,0,0,    # age
                     0,0,0,0,    # cog_reserve
                     1,0,0,0,    # lesions_wm_parcels
                     1,1,1,0),   # cog_ability
      my_blocks   =list(2,             # age
                        "res_cols",    # cog_reserve
                        "wm_cols",     # lesions_wm_parcels
                        "cog_cols"),   # cog_ability
      my_modes    =rep("A",4),
      my_datasets =all_dataset_names)

# Age, cog, WB measures, cog res
plspms_all_datasets[["cog_wb_age_res"]] <- create_plspms(
      my_variables=c("Age",        "Cog_Res",
                     "WB_Lesions", "Cognition"),
      my_matrix   =c(0,0,0,0,    # age
                     0,0,0,0,    # cog_reserve
                     1,0,0,0,    # whole_brain_measures
                     1,1,1,0),   # cog_ability
      my_blocks   =list(2,             # age
                        "res_cols",    # cog_reserve
                        "wb_cols",     # whole_brain_measures
                        "cog_cols"),   # cog_ability
      my_modes    =rep("A",4),
      my_datasets =all_dataset_names)


# Create df for export of PLS-PM info ----

# Create top-level list for PLS-PM models
plspm_outputs_for_export <- list()

# Now loop through model types adding info
plspm_model_types <- c("cog_gmp_age_res", "cog_wmp_age_res", "cog_wb_age_res")
dataset_names <- c("Combined",
                   "ADNI",      "Bialystok",      "Spreng",
                   "ADNI_only", "Bialystok_only", "Spreng_only")
for (model_type in 1:3)
{
  cur_model_type <- plspm_model_types[model_type]
  plspm_outputs_for_export[[cur_model_type]] <- as.data.frame(matrix(NA, nrow=7, ncol=40))
  colnames(plspm_outputs_for_export[[cur_model_type]]) <- c(
        "dataset",
        "path_AgeToBrain_ci_lower", "path_AgeToBrain_coef", "path_AgeToBrain_ci_upper",
        "path_AgeToCog_ci_lower",   "path_AgeToCog_coef",   "path_AgeToCog_ci_upper",
        "path_ResToCog_ci_lower",   "path_ResToCog_coef",   "path_ResToCog_ci_upper",
        "path_BrainToCog_ci_lower", "path_BrainToCog_coef", "path_BrainToCog_ci_upper",
        "reg_brain_r_sq_ci_lower",  "reg_brain_r_sq_coef",  "reg_brain_r_sq_ci_upper",
        "reg_cog_r_sq_ci_lower",    "reg_cog_r_sq_coef",    "reg_cog_r_sq_ci_upper",
        "reg_AgeToBrain_coef", "reg_AgeToBrain_se", "reg_AgeToBrain_t_val", "reg_AgeToBrain_p_val",
        "reg_AgeToCog_coef",   "reg_AgeToCog_se",   "reg_AgeToCog_t_val",   "reg_AgeToCog_p_val",
        "reg_ResToCog_coef",   "reg_ResToCog_se",   "reg_ResToCog_t_val",   "reg_ResToCog_p_val",
        "reg_BrainToCog_coef", "reg_BrainToCog_se", "reg_BrainToCog_t_val", "reg_BrainToCog_p_val",
        "unidim_brain_cron_alpha", "unidim_brain_dg_rho", "unidim_brain_eig1st", "unidim_brain_eig2nd",
        "goodness_of_fit")
  rownames(plspm_outputs_for_export[[cur_model_type]]) <- dataset_names
  for (dataset in 1:7)
  {
    # Path coefficients
    for (path in 1:4)
    {
      my_ci_lower <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["boot"]][["paths"]][["perc.025"]][path]
      my_coef     <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["boot"]][["paths"]][["Original"]][path]
      my_ci_upper <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["boot"]][["paths"]][["perc.975"]][path]
      plspm_outputs_for_export[[cur_model_type]][dataset, 2 + ((path-1) * 3)] <- my_ci_lower
      plspm_outputs_for_export[[cur_model_type]][dataset, 3 + ((path-1) * 3)] <- my_coef
      plspm_outputs_for_export[[cur_model_type]][dataset, 4 + ((path-1) * 3)] <- my_ci_upper
    }
    
    # Regression results - R-squared
    for (regr in 1:2)
    {
      my_ci_lower <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["boot"]][["rsq"]][["perc.025"]][regr]
      my_coef     <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["boot"]][["rsq"]][["Original"]][regr]
      my_ci_upper <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["boot"]][["rsq"]][["perc.975"]][regr]
      plspm_outputs_for_export[[cur_model_type]][dataset, 14 + ((regr-1) * 3)] <- my_ci_lower
      plspm_outputs_for_export[[cur_model_type]][dataset, 15 + ((regr-1) * 3)] <- my_coef
      plspm_outputs_for_export[[cur_model_type]][dataset, 16 + ((regr-1) * 3)] <- my_ci_upper
    }
    
    # Regression results - info for individual model parameters
    # Regression for brain LV
    for (regr_param in 1:1)
    {
      my_coef  <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[1]][1 + (regr_param)]
      my_se    <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[1]][3 + (regr_param)]
      my_t_val <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[1]][5 + (regr_param)]
      my_p_val <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[1]][7 + (regr_param)]
      plspm_outputs_for_export[[cur_model_type]][dataset, 20] <- my_coef
      plspm_outputs_for_export[[cur_model_type]][dataset, 21] <- my_se
      plspm_outputs_for_export[[cur_model_type]][dataset, 22] <- my_t_val
      plspm_outputs_for_export[[cur_model_type]][dataset, 23] <- my_p_val
    }
    # Regression for cognition LV
    for (regr_param in 1:3)
    {
      my_coef  <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[2]][1  + (regr_param)]
      my_se    <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[2]][5  + (regr_param)]
      my_t_val <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[2]][9  + (regr_param)]
      my_p_val <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["inner_model"]][[2]][13 + (regr_param)]
      plspm_outputs_for_export[[cur_model_type]][dataset, 24 + ((regr_param-1) * 4)] <- my_coef
      plspm_outputs_for_export[[cur_model_type]][dataset, 25 + ((regr_param-1) * 4)] <- my_se
      plspm_outputs_for_export[[cur_model_type]][dataset, 26 + ((regr_param-1) * 4)] <- my_t_val
      plspm_outputs_for_export[[cur_model_type]][dataset, 27 + ((regr_param-1) * 4)] <- my_p_val
    }
    
    # Unidimensionality indicators
    my_cron_alpha <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["unidim"]][["C.alpha"]][3]
    my_dg_rho     <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["unidim"]][["DG.rho"]][3]
    my_eig1st     <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["unidim"]][["eig.1st"]][3]
    my_eig2nd     <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]][["unidim"]][["eig.2nd"]][3]
    plspm_outputs_for_export[[cur_model_type]][dataset, 36] <- my_cron_alpha
    plspm_outputs_for_export[[cur_model_type]][dataset, 37] <- my_dg_rho
    plspm_outputs_for_export[[cur_model_type]][dataset, 38] <- my_eig1st
    plspm_outputs_for_export[[cur_model_type]][dataset, 39] <- my_eig2nd
    
    # Goodness of fit
    my_gof <- plspms_all_datasets[[cur_model_type]][["Models"]][[dataset]]["gof"]
    plspm_outputs_for_export[[cur_model_type]][dataset, 40] <- my_gof
  }
}

# Export the PLS-PM results
setwd(output_folder)
write.csv(plspm_outputs_for_export[["cog_gmp_age_res"]], "PLSPM_Results_GM_Models.csv")
write.csv(plspm_outputs_for_export[["cog_wmp_age_res"]], "PLSPM_Results_WM_Models.csv")
write.csv(plspm_outputs_for_export[["cog_wb_age_res"]],  "PLSPM_Results_WB_Models.csv")
setwd("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Data")


# Listing GM parcels in order of model loadings, by dataset ----

# Combined model
gm_loadings_combined <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["outer_model"]][["name"]][3:137],
                                   loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["outer_model"]][["loading"]][3:137])
gm_loadings_combined <- gm_loadings_combined[order(gm_loadings_combined$loading),]
gm_loadings_combined <- apply(gm_loadings_combined, 2, rev)   # Sort from highest to lowest instead of vice versa

# ADNI model
gm_loadings_adni <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI"]][["outer_model"]][["name"]][3:28],
                               loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI"]][["outer_model"]][["loading"]][3:28])
gm_loadings_adni <- gm_loadings_adni[order(gm_loadings_adni$loading),]
gm_loadings_adni <- apply(gm_loadings_adni, 2, rev)   # Sort from highest to lowest instead of vice versa

# Bialystok model
gm_loadings_bial <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok"]][["outer_model"]][["name"]][3:23],
                               loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok"]][["outer_model"]][["loading"]][3:23])
gm_loadings_bial <- gm_loadings_bial[order(gm_loadings_bial$loading),]
gm_loadings_bial <- apply(gm_loadings_bial, 2, rev)   # Sort from highest to lowest instead of vice versa

# Spreng model
gm_loadings_spreng <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng"]][["outer_model"]][["name"]][3:24],
                                 loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng"]][["outer_model"]][["loading"]][3:24])
gm_loadings_spreng <- gm_loadings_spreng[order(gm_loadings_spreng$loading),]
gm_loadings_spreng <- apply(gm_loadings_spreng, 2, rev)   # Sort from highest to lowest instead of vice versa

# ADNI only model
gm_loadings_adni_only <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI_only"]][["outer_model"]][["name"]][3:28],
                               loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI_only"]][["outer_model"]][["loading"]][3:28])
gm_loadings_adni_only <- gm_loadings_adni_only[order(gm_loadings_adni_only$loading),]
gm_loadings_adni_only <- apply(gm_loadings_adni_only, 2, rev)   # Sort from highest to lowest instead of vice versa

# Bialystok only model
gm_loadings_bial_only <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok_only"]][["outer_model"]][["name"]][5:25],
                               loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok_only"]][["outer_model"]][["loading"]][5:25])
gm_loadings_bial_only <- gm_loadings_bial_only[order(gm_loadings_bial_only$loading),]
gm_loadings_bial_only <- apply(gm_loadings_bial_only, 2, rev)   # Sort from highest to lowest instead of vice versa

# Spreng only model
gm_loadings_spreng_only <- data.frame(parcel=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng_only"]][["outer_model"]][["name"]][3:24],
                               loading=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng_only"]][["outer_model"]][["loading"]][3:24])
gm_loadings_spreng_only <- gm_loadings_spreng_only[order(gm_loadings_spreng_only$loading),]
gm_loadings_spreng_only <- apply(gm_loadings_spreng_only, 2, rev)   # Sort from highest to lowest instead of vice versa

# Export the sorted lists
setwd(output_folder)
write.csv(plspm_outputs_for_export[["gm_loadings_combined"]],    "GM_Loadings_Combined.csv")
write.csv(plspm_outputs_for_export[["gm_loadings_adni"]],        "GM_Loadings_ADNI.csv")
write.csv(plspm_outputs_for_export[["gm_loadings_bial"]],        "GM_Loadings_Bialystok.csv")
write.csv(plspm_outputs_for_export[["gm_loadings_spreng"]],      "GM_Loadings_Spreng.csv")
write.csv(plspm_outputs_for_export[["gm_loadings_adni_only"]],   "GM_Loadings_ADNI_only.csv")
write.csv(plspm_outputs_for_export[["gm_loadings_bial_only"]],   "GM_Loadings_Bialystok_only.csv")
write.csv(plspm_outputs_for_export[["gm_loadings_spreng_only"]], "GM_Loadings_Spreng_only.csv")
setwd("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Data")


# ====== Descriptive stats and visualizations ====== ----

# Check descriptive stats by dataset ----

# Data frame for descriptive stats
descr_stats_by_dataset <- as.data.frame(matrix(NA, nrow=4, ncol=10))
rownames(descr_stats_by_dataset) <- c("Combined", "ADNI", "Bialystok", "Spreng")
colnames(descr_stats_by_dataset) <- c("N_after_vis_excl", "N_final",
                                      "age_mean",     "age_sd",     "sex_f_n",   "sex_m_n",
                                      "edu_yrs_mean", "edu_yrs_sd", "mmse_mean", "mmse_sd")

# Number of participants
descr_stats_by_dataset["Combined", "N_after_vis_excl"] <- nrow(df_combined_mri_means)
descr_stats_by_dataset["Combined", "N_final"]          <- nrow(df_plspm_combined)
descr_stats_by_dataset["ADNI",     "N_after_vis_excl"] <- nrow(df_combined_mri_means[df_combined_mri_means$dataset == "ADNI",])
descr_stats_by_dataset["ADNI",     "N_final"]          <- nrow(df_plspm_combined[df_plspm_combined$dataset == "ADNI",])
descr_stats_by_dataset["Bialystok","N_after_vis_excl"] <- nrow(df_combined_mri_means[df_combined_mri_means$dataset == "Bialystok",])
descr_stats_by_dataset["Bialystok","N_final"]          <- nrow(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",])
descr_stats_by_dataset["Spreng",   "N_after_vis_excl"] <- nrow(df_combined_mri_means[df_combined_mri_means$dataset == "Spreng",])
descr_stats_by_dataset["Spreng",   "N_final"]          <- nrow(df_plspm_combined[df_plspm_combined$dataset == "Spreng",])

# Age
descr_stats_by_dataset["Combined", "age_mean"] <- mean(df_plspm_combined$age, na.rm=TRUE)
descr_stats_by_dataset["Combined", "age_sd"]   <- sd(df_plspm_combined$age, na.rm=TRUE)
descr_stats_by_dataset["ADNI",     "age_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "ADNI",]$age, na.rm=TRUE)
descr_stats_by_dataset["ADNI",     "age_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "ADNI",]$age, na.rm=TRUE)
descr_stats_by_dataset["Bialystok","age_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",]$age, na.rm=TRUE)
descr_stats_by_dataset["Bialystok","age_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",]$age, na.rm=TRUE)
descr_stats_by_dataset["Spreng",   "age_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Spreng",]$age, na.rm=TRUE)
descr_stats_by_dataset["Spreng",   "age_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "Spreng",]$age, na.rm=TRUE)

# Sex
descr_stats_by_dataset["Combined", "sex_f_n"]  <- sum(as.numeric(combined_df_reduced$sex == 0))    # Female
descr_stats_by_dataset["Combined", "sex_m_n"]  <- sum(as.numeric(combined_df_reduced$sex == 1))    # Male
descr_stats_by_dataset["ADNI",     "sex_f_n"]  <- sum(as.numeric(combined_df_reduced[combined_df_reduced$dataset == "ADNI",]$sex == 0))       # Female
descr_stats_by_dataset["ADNI",     "sex_m_n"]  <- sum(as.numeric(combined_df_reduced[combined_df_reduced$dataset == "ADNI",]$sex == 1))       # Male
descr_stats_by_dataset["Bialystok","sex_f_n"]  <- sum(as.numeric(combined_df_reduced[combined_df_reduced$dataset == "Bialystok",]$sex == 0))  # Female
descr_stats_by_dataset["Bialystok","sex_m_n"]  <- sum(as.numeric(combined_df_reduced[combined_df_reduced$dataset == "Bialystok",]$sex == 1))  # Male
descr_stats_by_dataset["Spreng",   "sex_f_n"]  <- sum(as.numeric(combined_df_reduced[combined_df_reduced$dataset == "Spreng",]$sex == 0))     # Female
descr_stats_by_dataset["Spreng",   "sex_m_n"]  <- sum(as.numeric(combined_df_reduced[combined_df_reduced$dataset == "Spreng",]$sex == 1))     # Male

# Education (yrs)
descr_stats_by_dataset["Combined", "edu_yrs_mean"] <- mean(df_plspm_combined$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["Combined", "edu_yrs_sd"]   <- sd(df_plspm_combined$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["ADNI",     "edu_yrs_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "ADNI",]$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["ADNI",     "edu_yrs_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "ADNI",]$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["Bialystok","edu_yrs_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",]$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["Bialystok","edu_yrs_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",]$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["Spreng",   "edu_yrs_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Spreng",]$yrs_education, na.rm=TRUE)
descr_stats_by_dataset["Spreng",   "edu_yrs_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "Spreng",]$yrs_education, na.rm=TRUE)

# MMSE
descr_stats_by_dataset["Combined", "mmse_mean"] <- mean(df_plspm_combined$mmse, na.rm=TRUE)
descr_stats_by_dataset["Combined", "mmse_sd"]   <- sd(df_plspm_combined$mmse, na.rm=TRUE)
descr_stats_by_dataset["ADNI",     "mmse_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "ADNI",]$mmse, na.rm=TRUE)
descr_stats_by_dataset["ADNI",     "mmse_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "ADNI",]$mmse, na.rm=TRUE)
descr_stats_by_dataset["Bialystok","mmse_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",]$mmse, na.rm=TRUE)
descr_stats_by_dataset["Bialystok","mmse_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "Bialystok",]$mmse, na.rm=TRUE)
descr_stats_by_dataset["Spreng",   "mmse_mean"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Spreng",]$mmse, na.rm=TRUE)
descr_stats_by_dataset["Spreng",   "mmse_sd"]   <- sd(df_plspm_combined[df_plspm_combined$dataset == "Spreng",]$mmse, na.rm=TRUE)

# Export the descriptive stats
setwd(output_folder)
write.csv(descr_stats_by_dataset, "DescriptiveStatsByDataset.csv")
setwd("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Data")


# ====== Display findings in brain space ====== ----

# Load files for Schaefer atlas ----

# GM atlas: Load the NIFTI file
schaefer_atlas_nifti <- RNifti::readNifti("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/NIFTI/100Parcels7Networks.nii")

# Load in the .csv file for the Schaefer atlas
schaefer_atlas_csv <- read_csv("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/NIFTI/100Parcels7Networks.csv")
schaefer_atlas_csv <- subset(schaefer_atlas_csv, select = -1)


# Function for plotting GM findings in brain space ----
# Output title should be a string ending in .nii
create_gm_nifti <- function(gm_parcel_names, gm_parcel_values, output_title)
{
  # Setting up standard numbering for the parcels
  gm_region_nums <- data.frame(parcel_num = 1:135, parcel_name = parcel_names_gm)
  gm_region_nums$parcel_name <- make_clean_names(gm_region_nums$parcel_name)
  
  # Importing the data to be visualized
  gm_names_values <- data.frame(parcel_name  = gm_parcel_names,
                                parcel_value = gm_parcel_values)
  
  # Adding the standard region numbers to the imported data
  joined_gm_df <- left_join(gm_region_nums, gm_names_values, by="parcel_name")

  # Create duplicate nifti array to update
  schaefer_atlas_nifti_updated <- schaefer_atlas_nifti

  # Add values from data to updated nifti file
  for (i in 1:(length(joined_gm_df$parcel_num)))
  {
    my_region_num <- joined_gm_df[i, "parcel_num"]
    region_value <- joined_gm_df[joined_gm_df$parcel_num == my_region_num, "parcel_value"]
    schaefer_atlas_nifti_updated[schaefer_atlas_nifti == my_region_num] <- region_value
  }

  # Create and save the new nifti file
  new_file_name <- paste0(output_folder, "/", output_title, ".nii")
  RNifti::writeNifti(schaefer_atlas_nifti_updated, new_file_name)

  # Confirmation message
  print(paste("Successfully created", new_file_name))
}


# Function to check against 2.5 and 97.5 percentiles ----

# Determine if a given finding is statistically reliable based on whether its 95% bootstrap CI includes 0
check_item_reliability <- function(low, high)
{
  if ((low > 0 && high > 0) || (low < 0 && high < 0))
  {
    return(1)
  }
  else
  {
    return(0)
  }
}

check_block_reliability <- function(block, cols)
{
  low_perc  <- block[["perc.025"]][cols]
  high_perc <- block[["perc.975"]][cols]
  is_reliable <- rep(NA, length(cols))
  for (i in 1:(length(cols)))
  {
    is_reliable[i] <- check_item_reliability(low_perc[i], high_perc[i]) 
  }
  return(is_reliable)
}


# Create .nii files of GM disruption per parcel, by dataset ----

# Create df of the mean disruption by parcel, with a row for each dataset
gm_means_by_dataset <- as.data.frame(matrix(NA, ncol=4, nrow=135))
colnames(gm_means_by_dataset) <- c("Combined", "ADNI", "Bialystok", "Spreng")
rownames(gm_means_by_dataset) <- parcel_names_gm
for (parcel in 1:135)
{
  gm_means_by_dataset[parcel, "Combined"]  <- mean(df_plspm_combined[, parcel+76])
  gm_means_by_dataset[parcel, "ADNI"]      <- mean(df_plspm_combined[df_plspm_combined$dataset == "ADNI", parcel+76])
  gm_means_by_dataset[parcel, "Bialystok"] <- mean(df_plspm_combined[df_plspm_combined$dataset == "Bialystok", parcel+76])
  gm_means_by_dataset[parcel, "Spreng"]    <- mean(df_plspm_combined[df_plspm_combined$dataset == "Spreng", parcel+76])
}

# GM disruption by parcel -- by dataset
create_gm_nifti(gm_parcel_names =parcel_names_gm,
                gm_parcel_values=as.numeric(as.vector(gm_means_by_dataset$Combined)),
                output_title="Combined_GM_Disruption")
create_gm_nifti(gm_parcel_names =parcel_names_gm,
                gm_parcel_values=as.numeric(as.vector(gm_means_by_dataset$ADNI)),
                output_title="ADNI_GM_Disruption")
create_gm_nifti(gm_parcel_names =parcel_names_gm,
                gm_parcel_values=as.numeric(as.vector(gm_means_by_dataset$Bialystok)),
                output_title="Bialystok_GM_Disruption")
create_gm_nifti(gm_parcel_names =parcel_names_gm,
                gm_parcel_values=as.numeric(as.vector(gm_means_by_dataset$Spreng)),
                output_title="Spreng_GM_Disruption")


# Create .nii files of model loadings by GM parcel ----

# Loadings for GM parcels -- Combined dataset
combined_gm_cols <- 3:137   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["outer_model"]][["name"]][combined_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["outer_model"]][["loading"]][combined_gm_cols],
                output_title="Combined_GM_Loadings")
combined_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["boot"]][["loadings"]],
                                                cols = combined_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["outer_model"]][["name"]][combined_gm_cols],
                gm_parcel_values=combined_gm_reliable,
                output_title="Combined_GM_Loadings_Reliable")

# Loadings for GM parcels -- ADNI
adni_gm_cols <- 3:28   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI"]][["outer_model"]][["name"]][adni_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI"]][["outer_model"]][["loading"]][adni_gm_cols],
                output_title="ADNI_GM_Loadings")
adni_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI"]][["boot"]][["loadings"]],
                                            cols = adni_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI"]][["outer_model"]][["name"]][adni_gm_cols],
                gm_parcel_values=adni_gm_reliable,
                output_title="ADNI_GM_Loadings_Reliable")

# Loadings for GM parcels -- Bialystok
bial_gm_cols <- 3:23   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok"]][["outer_model"]][["name"]][bial_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok"]][["outer_model"]][["loading"]][bial_gm_cols],
                output_title="Bialystok_GM_Loadings")
bial_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok"]][["boot"]][["loadings"]],
                                            cols = bial_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok"]][["outer_model"]][["name"]][bial_gm_cols],
                gm_parcel_values=bial_gm_reliable,
                output_title="Bialystok_GM_Loadings_Reliable")

# Loadings for GM parcels -- Spreng
spreng_gm_cols <- 3:24   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng"]][["outer_model"]][["name"]][spreng_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng"]][["outer_model"]][["loading"]][spreng_gm_cols],
                output_title="Spreng_GM_Loadings")
spreng_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng"]][["boot"]][["loadings"]],
                                              cols = spreng_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng"]][["outer_model"]][["name"]][spreng_gm_cols],
                gm_parcel_values=spreng_gm_reliable,
                output_title="Spreng_GM_Loadings_Reliable")

# Loadings for GM parcels -- ADNI only
adni_only_gm_cols <- 3:28   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI_only"]][["outer_model"]][["name"]][adni_only_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI_only"]][["outer_model"]][["loading"]][adni_only_gm_cols],
                output_title="ADNI_Only_GM_Loadings")
adni_only_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI_only"]][["boot"]][["loadings"]],
                                                 cols = adni_only_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["ADNI_only"]][["outer_model"]][["name"]][adni_only_gm_cols],
                gm_parcel_values=adni_only_gm_reliable,
                output_title="ADNI_Only_GM_Loadings_Reliable")

# Loadings for GM parcels -- Bialystok only
bial_only_gm_cols    <- 5:25   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok_only"]][["outer_model"]][["name"]][bial_only_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok_only"]][["outer_model"]][["loading"]][bial_only_gm_cols],
                output_title="Bialystok_Only_GM_Loadings")
bial_only_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Combined"]][["boot"]][["loadings"]],
                                                 cols = bial_only_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Bialystok_only"]][["outer_model"]][["name"]][bial_only_gm_cols],
                gm_parcel_values=bial_only_gm_reliable,
                output_title="Bialystok_Only_GM_Loadings_Reliable")

# Loadings for GM parcels -- Spreng only
spreng_only_gm_cols <- 3:24   # Must check this
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng_only"]][["outer_model"]][["name"]][spreng_only_gm_cols],
                gm_parcel_values=plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng_only"]][["outer_model"]][["loading"]][spreng_only_gm_cols],
                output_title="Spreng_Only_GM_Loadings")
spreng_only_gm_reliable <- check_block_reliability(block = plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng_only"]][["boot"]][["loadings"]],
                                                   cols = spreng_only_gm_cols)
create_gm_nifti(gm_parcel_names =plspms_all_datasets[["cog_gmp_age_res"]][["Models"]][["Spreng_only"]][["outer_model"]][["name"]][spreng_only_gm_cols],
                gm_parcel_values=spreng_only_gm_reliable,
                output_title="Spreng_Only_GM_Loadings_Reliable")



