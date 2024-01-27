# Project on white matter hyperintensities and cognitive reserve in older adults
# MRI data has already been through several stages of analysis

# Note to self: Collapse all with Alt+O

# Preliminaries ----

setwd("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_CR_Project/Data")

# Installation instructions for ADNIMERGE if needed in the future
# install.packages("C:/Users/Arthur/OneDrive - Carleton University/Documents/R/WMH_Project_LQT_Output_Analysis/Libraries/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")

# Load libraries
library(tidyverse)
library(tidyr)
library(janitor)       # Needed for clean_names function
library(Hmisc)         # Needed to run ADNIMERGE
library(ADNIMERGE)
library(data.table)    # To transpose data frame
library(grid)          # To display multiple plots in a grid
library(gridExtra)     # Extension to the grid library
library(GGally)        # For correlogram

# For displaying brain images (though may replace with others)
# library(ggseg)
# library(ggsegYeo2011)
# Code for ggseg that doesn't appear to be necessary
# options(repos = c(
#     ggseg = 'https://ggseg.r-universe.dev',
#     CRAN = 'https://cloud.r-project.org'))
# install.packages("ggsegYeo2011")

# For PLS
# install.packages("TExPosition")
# install.packages("remotes")   # Used for the installation
# library(remotes)
# remotes::install_github("HerveAbdi/data4PCCAR")
# remotes::install_github("HerveAbdi/PTCA4CATA")
# library(TExPosition)
# library(data4PCCAR)
# library(PTCA4CATA)
# Note: John also had library(TInPosition) but it doesn't appear to be needed

# For PLS-PM
library(plspm)

# For harmonization
# library(devtools)
# install_github("jfortin1/neuroCombatData")
# install_github("jfortin1/neuroCombat_Rpackage")
library(neuroCombatData)
library(neuroCombat)

# For checking near-zero variance
library(caret)

# Loading and cleaning data: Functions ----

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
  out <- import_csv_data(my_file_name, FALSE, id_col, 5, 7, "Spreng_", spreng_to_exclude, "Spreng")
}

import_bial_cog <- function(my_file_name, good_cols)
{
  my_df <- import_csv_data(my_file_name, TRUE, "id", 1, 3, "Bialystok_", bial_to_exclude, "Bialystok")
  my_df <- my_df[,good_cols]   # Note: Remember to add 1 to the end of the range for good_cols to account for the "dataset" column added earlier
  out <- my_df
}


# Loading and cleaning data: MRI analysis data ----

# Participants to exclude based on visual examination of their lesion images
adni_to_exclude   <- c("ADNI_003_S_6258", "ADNI_003_S_6260", "ADNI_005_S_4185", "ADNI_009_S_4324",
                       "ADNI_019_S_4835", "ADNI_021_S_6312", "ADNI_021_S_6896", "ADNI_022_S_6796",
                       "ADNI_022_S_6863", "ADNI_027_S_2219", "ADNI_027_S_4919", "ADNI_027_S_6733",
                       "ADNI_027_S_6788", "ADNI_029_S_2395", "ADNI_029_S_4290", "ADNI_029_S_5219",
                       "ADNI_029_S_6726", "ADNI_033_S_4176", "ADNI_033_S_7066", "ADNI_057_S_6869",
                       "ADNI_098_S_0896", "ADNI_098_S_6593", "ADNI_099_S_6038", "ADNI_099_S_6396",
                       "ADNI_099_S_6632", "ADNI_114_S_6597", "ADNI_123_S_0106", "ADNI_123_S_6891",
                       "ADNI_126_S_0680", "ADNI_126_S_4507", "ADNI_126_S_4514", "ADNI_126_S_5243",
                       "ADNI_126_S_6683", "ADNI_126_S_6724", "ADNI_126_S_7015", "ADNI_126_S_7060",
                       "ADNI_127_S_0259", "ADNI_127_S_1427", "ADNI_127_S_4197", "ADNI_127_S_5200",
                       "ADNI_127_S_6330", "ADNI_127_S_6512", "ADNI_128_S_2002", "ADNI_128_S_2220",
                       "ADNI_128_S_4607", "ADNI_129_S_4422", "ADNI_129_S_6288", "ADNI_129_S_6763",
                       "ADNI_129_S_6784", "ADNI_129_S_6852", "ADNI_130_S_4294", "ADNI_130_S_4417",
                       "ADNI_130_S_6604", "ADNI_135_S_6509", "ADNI_135_S_6544", "ADNI_135_S_6545",
                       "ADNI_135_S_6586", "ADNI_135_S_6622", "ADNI_168_S_6591", "ADNI_168_S_6828",
                       "ADNI_301_S_6615", "ADNI_305_S_6850", "ADNI_941_S_6345", "ADNI_941_S_6962")
bial_to_exclude   <- c()
spreng_to_exclude <- c("Spreng_109")

# Load data frames containing MRI analysis outputs for ADNI dataset
df_adni_wm_lesions   <- import_adni_mri("ADNI_lesion_volumes_numbers_20231118.csv", "participant")
df_adni_gm_lesions   <- import_adni_mri("ADNI_gray_matter_lesions_20231116.csv", "filename")
df_adni_wm_discon    <- import_adni_mri("ADNI_percent_discon_tracts_20231116.csv", "filename")
# df_adni_gm_discon    <- import_adni_mri("ADNI_gray_matter_discon_20231116.csv", "filename")
# df_adni_delta_sspl   <- import_adni_mri("ADNI_delta_SSPL_20231116.csv", "filename")

# Load data frames containing MRI analysis outputs for Bialystok dataset
df_bial_wm_lesions   <- import_bial_mri("Bialystok_lesion_volumes_numbers_20231118.csv", "participant")
df_bial_gm_lesions   <- import_bial_mri("Bialystok_gray_matter_lesions_20231113.csv", "filename")
df_bial_wm_discon    <- import_bial_mri("Bialystok_percent_discon_tracts_20231113.csv", "filename")
# df_bial_gm_discon    <- import_bial_mri("Bialystok_gray_matter_discon_20231113.csv", "filename")
# df_bial_delta_sspl   <- import_bial_mri("Bialystok_delta_SSPL_20231113.csv", "filename")

# Load data frames containing MRI analysis outputs for Spreng dataset
df_spreng_gm_lesions <- import_spreng_mri("Spreng_gray_matter_lesions_20231113.csv", "filename")
df_spreng_wm_discon  <- import_spreng_mri("Spreng_percent_discon_tracts_20231113.csv", "filename")
# df_spreng_gm_discon  <- import_spreng_mri("Spreng_gray_matter_discon_20231113.csv", "filename")
# df_spreng_delta_sspl <- import_spreng_mri("Spreng_delta_SSPL_20231113.csv", "filename")


# Doing the WM lesion summaries for the Spreng dataset separately because the IDs are in a bad format
df_spreng_wm_lesions <- read.csv(file = "Spreng_lesion_volumes_numbers_20231118.csv", encoding = "UTF-8")
df_spreng_wm_lesions <- clean_names(df_spreng_wm_lesions)
colnames(df_spreng_wm_lesions)[which(names(df_spreng_wm_lesions) == "participant")] <- "ptid"
# df_spreng_wm_lesions <- df_spreng_wm_lesions %>% rename("ptid" = "participant")
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

# Importing parcel names
parcel_names_wm <- read.csv(file = "Parcel_names_WM.csv", encoding = "UTF-8", header = FALSE)[,1]
parcel_names_gm <- read.csv(file = "Parcel_names_GM.csv", encoding = "UTF-8", header = FALSE)[,1]
# Adding parcel names to data frames (and cleaning the names)
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


# Loading and cleaning data: Non-MRI data ----

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
# This one has aggregate scores for cognitive tasks in overall domains
data(uwnpsychsum)
df_adni_uwnpsychsum <- uwnpsychsum
rm(uwnpsychsum)
df_adni_uwnpsychsum <- clean_names(df_adni_uwnpsychsum)
df_adni_uwnpsychsum <- df_adni_uwnpsychsum[df_adni_uwnpsychsum$colprot == "ADNI3", ]
df_adni_adnimerge$rid <- as.numeric(df_adni_adnimerge$rid)
df_adni_adnimerge$viscode <- as.character(df_adni_adnimerge$viscode)
df_adni_uwnpsychsum$viscode <- as.character(df_adni_uwnpsychsum$viscode)
df_adni_uwnpsychsum <- inner_join(df_adni_uwnpsychsum, df_adni_adnimerge, by = c("rid", "viscode"))
df_adni_uwnpsychsum <- df_adni_uwnpsychsum[order(df_adni_uwnpsychsum$rid),]
df_adni_uwnpsychsum <- df_adni_uwnpsychsum %>% mutate(dataset = "ADNI") %>% select(rid, dataset, everything())

# Loading and cleaning cognitive/demographic/etc data from sheets from the Bialystok .xls file
df_bial_mmse <- import_bial_cog("Bialystok_cog_data_MMSE_20231113.csv", 1:12)
df_bial_lsbq <- import_bial_cog("Bialystok_cog_data_LSBQ_20231113.csv", 1:129)
df_bial_criq <- import_bial_cog("Bialystok_cog_data_CRIQ_20231113.csv", 1:59)

# Loading and cleaning cognitive/demographic/etc data from Spreng dataset
df_spreng_cog <- read.csv(file = "Spreng_cog_data_20231016.csv", encoding = "UTF-8")
df_spreng_cog <- clean_names(df_spreng_cog)
# df_spreng_cog <- df_spreng_cog %>% rename("ptid" = "id")
colnames(df_spreng_cog)[which(names(df_spreng_cog) == "id")] <- "ptid"
for (i in 1:99) {
  df_spreng_cog[i,"ptid"] <- paste(c(substring(df_spreng_cog[i,"ptid"], 1, 4), "0", substring(df_spreng_cog[i,"ptid"], 5, 6)), collapse="")
}
df_spreng_cog$ptid <- substring(df_spreng_cog$ptid, 5, 7)
df_spreng_cog$ptid <- paste("Spreng_", df_spreng_cog$ptid, sep="")
df_spreng_cog <- filter(df_spreng_cog, ptid %in% df_spreng_gm_lesions$ptid)
df_spreng_cog <- df_spreng_cog[order(df_spreng_cog$ptid),]
df_spreng_cog <- filter(df_spreng_cog, !ptid %in% spreng_to_exclude)
df_spreng_cog <- df_spreng_cog %>% mutate(dataset = "Spreng") %>% select(ptid, dataset, everything())


# Loading and cleaning data: Remove participants from the "Young" group from the Spreng dataset  ----
df_spreng_cog        <- df_spreng_cog[df_spreng_cog$agegroup != "Y",]
df_spreng_wm_lesions <- filter(df_spreng_wm_lesions, ptid %in% df_spreng_cog$ptid)
df_spreng_gm_lesions <- filter(df_spreng_gm_lesions, ptid %in% df_spreng_cog$ptid)
df_spreng_wm_discon  <- filter(df_spreng_wm_discon,  ptid %in% df_spreng_cog$ptid)
# df_spreng_gm_discon  <- filter(df_spreng_gm_discon,  ptid %in% df_spreng_cog$ptid)
# df_spreng_delta_sspl <- filter(df_spreng_delta_sspl, ptid %in% df_spreng_cog$ptid)

# Maybe this part should be integrated into the Spreng loading function above?


# Loading and cleaning data: Keep only CN participants from the ADNI dataset ----

df_adni_adnimerge   <- df_adni_adnimerge[df_adni_adnimerge$dx == "CN",]
df_adni_adnimerge   <- df_adni_adnimerge[!is.na(df_adni_adnimerge$ptid),]
df_adni_uwnpsychsum <- filter(df_adni_uwnpsychsum, ptid %in% df_adni_adnimerge$ptid)
df_adni_wm_lesions  <- filter(df_adni_wm_lesions, ptid %in% df_adni_adnimerge$ptid)
df_adni_wm_discon   <- filter(df_adni_wm_discon, ptid %in% df_adni_adnimerge$ptid)
df_adni_gm_lesions  <- filter(df_adni_gm_lesions, ptid %in% df_adni_adnimerge$ptid)
# df_adni_gm_discon   <- filter(df_adni_gm_discon, ptid %in% df_adni_adnimerge$ptid)
# df_adni_delta_sspl  <- filter(df_adni_delta_sspl, ptid %in% df_adni_adnimerge$ptid)

# Maybe this part should be integrated into the ADNI loading function above?


# Loading and cleaning data: Combining data frames ----

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
rm(df_adni_cog_to_combine)
rm(df_bial_cog_to_combine)
rm(df_spreng_cog_to_combine)

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
# df_combined_gm_discon  <- left_join(df_combined_cog,
#                                     rbind(df_adni_gm_discon,
#                                           df_bial_gm_discon,
#                                           df_spreng_gm_discon),
#                                     by = c("ptid", "dataset"))
# df_combined_delta_sspl <- left_join(df_combined_cog,
#                                     rbind(df_adni_delta_sspl,
#                                           df_bial_delta_sspl,
#                                           df_spreng_delta_sspl),
#                                     by = c("ptid", "dataset"))

# One data frame with all participants showing means of all MRI measure types
df_combined_mri_means <- df_combined_wm_lesions
df_combined_mri_means <- mutate(df_combined_mri_means,
                                gm_lesions_mean = rowMeans(select(df_combined_gm_lesions, 8:142),
                                                           na.rm = TRUE))
df_combined_mri_means <- mutate(df_combined_mri_means,
                                wm_discon_mean = rowMeans(select(df_combined_wm_discon, 8:77),
                                                          na.rm = TRUE))
# df_combined_mri_means <- mutate(df_combined_mri_means,
#                                 gm_discon_mean = rowMeans(select(df_combined_gm_discon, 8:9052),
#                                                           na.rm = TRUE))
# df_combined_mri_means <- mutate(df_combined_mri_means,
#                                 delta_sspl_mean = rowMeans(select(df_combined_delta_sspl, 8:9052),
#                                                            na.rm = TRUE))

# df_combined_all <- left_join(df_combined_cog, df_combined_mri_means, by = c("ptid", "dataset"))


# Visualizing descriptives: Creating functions for creating plots in a common style ----

# Function to create histogram in this style
# Need to use "binwidth" parameter to improve appearance
create_histogram <- function(my_data, my_variable, my_title, my_facet)
{
  ggplot(data=my_data, aes(x={{my_variable}}, rm.na = TRUE)) +
         geom_histogram(fill='blue', color='lightgrey') + facet_wrap(vars({{my_facet}})) +
         ylab("Number of participants") + ggtitle(my_title) +
         theme_bw() + theme(plot.title=element_text(hjust=0.5))
}

# Function to create scatterplot in this style (needs to be fixed -- I think double curly braces)
create_scatterplot <- function(my_data, my_x, my_y, my_facet, my_title)
{
ggplot(data=my_data, aes(x={{my_x}}, y={{my_y}}, col={{my_facet}}, rm.na=TRUE)) +
       geom_point(alpha=0.5) + geom_smooth(method=lm) + facet_wrap(vars({{my_facet}})) +
       ggtitle(my_title) + theme_bw() + theme(plot.title=element_text(hjust=0.5))
}


# # Visualizing descriptives: Plotting descriptive statistics and various relationships without considering parcel-level differences ----
#
# # Histograms for non-MRI measures
# create_histogram(df_combined_cog, age, "Age by dataset", dataset)
# create_histogram(df_combined_cog, yrs_education, "Education (yrs) by dataset", dataset)
# create_histogram(df_combined_cog, mmse, "MMSE scores by dataset", dataset)
# 
# # Histograms for MRI measures
# create_histogram(df_combined_mri_means, number_of_lesions, "Number of WM lesions", dataset)
# create_histogram(df_combined_mri_means, lesion_volume_in_m_l, "WM lesion volume (mL)", dataset)
# create_histogram(df_combined_mri_means, gm_lesions_mean, "GM mean lesion volume (inferred)", dataset)
# create_histogram(df_combined_mri_means, wm_discon_mean, "WM mean disconectivity", dataset)
# create_histogram(df_combined_mri_means, gm_discon_mean, "GM mean disconnectivity", dataset)
# create_histogram(df_combined_mri_means, delta_sspl_mean, "Mean change in SSPL", dataset)
# 
# # Correlogram for all MRI measure means
# ggpairs(df_combined_mri_means[,8:13], title="Correlogram of MRI Measure Means",
#         upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
#         lower = list(continuous = wrap("smooth", color='blue', alpha=0.25), combo = "facethist", discrete = "facetbar", na = "na"),
#         diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
#         columnLabels = c("WM Les Vol", "WM Les Num", "GM Les Vol", "WM Discon", "GM Discon", "SSPL"), axisLabels = "none") +
#         theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# # Correlations between first four MRI measures (but not the two using pairs of parcels)
# create_scatterplot(df_combined_mri_means, log(number_of_lesions), log(lesion_volume_in_m_l),
#                    dataset, "Number of lesions vs. lesion volume")
# create_scatterplot(df_combined_mri_means, log(number_of_lesions), log(gm_lesions_mean),
#                    dataset, "Number of lesions vs. mean GM lesions")
# create_scatterplot(df_combined_mri_means, log(number_of_lesions), log(wm_discon_mean),
#                    dataset, "Number of lesions vs. mean WM disconnection")
# create_scatterplot(df_combined_mri_means, log(lesion_volume_in_m_l), log(gm_lesions_mean),
#                    dataset, "Lesion volume vs. mean GM lesions")
# create_scatterplot(df_combined_mri_means, log(lesion_volume_in_m_l), log(wm_discon_mean),
#                    dataset, "Lesion volume vs. mean WM disconnection")
# create_scatterplot(df_combined_mri_means, log(gm_lesions_mean), log(wm_discon_mean),
#                    dataset, "Mean GM lesions vs. mean WM disconnection")
# 
# # Histograms grouped by sex rather than dataset
# create_histogram(df_combined_mri_means, lesion_volume_in_m_l, "WM lesion volume (mL) by sex", sex)
# create_histogram(df_combined_mri_means, gm_lesions_mean, "GM mean lesion volume by sex", sex)
# create_histogram(df_combined_mri_means, wm_discon_mean, "WM mean disconection by sex", sex)
# create_histogram(df_combined_mri_means, mmse, "MMSE scores by sex", sex)
# 
# # Correlations between age and MRI measures
# create_scatterplot(df_combined_mri_means, age, log(number_of_lesions),
#                    dataset, "Age vs. number of lesions")
# create_scatterplot(df_combined_mri_means, age, log(lesion_volume_in_m_l),
#                    dataset, "Age vs. WM lesion volume (mL)")
# create_scatterplot(df_combined_mri_means, age, log(gm_lesions_mean),
#                    dataset, "Age vs. mean GM lesion volume")
# create_scatterplot(df_combined_mri_means, age, log(wm_discon_mean),
#                    dataset, "Age vs. mean WM disconnection")
# create_scatterplot(df_combined_mri_means, age, log(gm_discon_mean),
#                    dataset, "Age vs. mean GM disconnection")
# create_scatterplot(df_combined_mri_means, age, log(delta_sspl_mean),
#                    dataset, "Age vs. change in SSPL")
# 
# # Correlations between education and MRI measures
# create_scatterplot(df_combined_mri_means, yrs_education, log(number_of_lesions),
#                    dataset, "Education vs. number of lesions")
# create_scatterplot(df_combined_mri_means, yrs_education, log(lesion_volume_in_m_l),
#                    dataset, "Education vs. WM lesion volume (mL)")
# create_scatterplot(df_combined_mri_means, yrs_education, log(gm_lesions_mean),
#                    dataset, "Education vs. mean GM lesion volume")
# create_scatterplot(df_combined_mri_means, yrs_education, log(wm_discon_mean),
#                    dataset, "Education vs. mean WM disconnection")
# create_scatterplot(df_combined_mri_means, yrs_education, log(gm_discon_mean),
#                    dataset, "Education vs. mean GM disconnection")
# create_scatterplot(df_combined_mri_means, yrs_education, log(delta_sspl_mean),
#                    dataset, "Education vs. change in SSPL")
# 
# # Correlations between MMSE and MRI measures
# create_scatterplot(df_combined_mri_means, mmse, log(number_of_lesions),
#                    dataset, "MMSE vs. number of lesions")
# create_scatterplot(df_combined_mri_means, mmse, log(lesion_volume_in_m_l),
#                    dataset, "MMSE vs. WM lesion volume (mL)")
# create_scatterplot(df_combined_mri_means, mmse, log(gm_lesions_mean),
#                    dataset, "MMSE vs. mean GM lesion volume")
# create_scatterplot(df_combined_mri_means, mmse, log(wm_discon_mean),
#                    dataset, "MMSE vs. mean WM disconnection")
# create_scatterplot(df_combined_mri_means, mmse, log(gm_discon_mean),
#                    dataset, "MMSE vs. mean GM disconnection")
# create_scatterplot(df_combined_mri_means, mmse, log(delta_sspl_mean),
#                    dataset, "MMSE vs. change in SSPL")
# 
# # Correlations between demographic measures and MMSE
# create_scatterplot(df_combined_mri_means, age, mmse,
#                    dataset, "Age vs. MMSE")
# create_scatterplot(df_combined_mri_means, yrs_education, mmse,
#                    dataset, "Age vs. education")


# PLS-PM: Creating reduced df, from which I will create the PLS-PM dfs ----

reduced_df_combined <- merge(df_combined_cog[,c(1:2,4:7)],
                             df_combined_wm_discon[,c(1,8:77)],
                             by="ptid")
reduced_df_combined <- merge(reduced_df_combined,
                             df_combined_gm_lesions[,c(1,8:142)],
                             by="ptid")
reduced_df_combined <- merge(reduced_df_combined,
                             df_combined_wm_lesions[,c(1,8:9)],
                             by="ptid")
reduced_df_combined$sex <- as.numeric(reduced_df_combined$sex)


# PLS-PM: Conducting harmonization for the combined df ----

# Preparing inputs for the harmonization function
trans_df_combined <- transpose(reduced_df_combined[,c(7:213)],
                               keep.names = NULL)
combined_batch_vector <- reduced_df_combined$dataset
combined_batch_vector[combined_batch_vector == "ADNI"]      <- 1
combined_batch_vector[combined_batch_vector == "Bialystok"] <- 2
combined_batch_vector[combined_batch_vector == "Spreng"]    <- 3
combined_batch_vector <- as.numeric(combined_batch_vector)
combined_mod_vector   <- reduced_df_combined$sex

# Running the harmonization itself
combined_harmonized_data <- neuroCombat(dat=trans_df_combined,
                                        batch=combined_batch_vector,
                                        mod=combined_mod_vector)

# Turning harmonization output into df
# combined_standardized_df <- as.data.frame(combined_harmonized_data[[5]])  # This one is standardized but not harmonized
combined_harmonized_df <- as.data.frame(combined_harmonized_data[[1]])

# Transpose the data frame back
combined_df_post_harm <- transpose(combined_harmonized_df, keep.names = NULL)

# Re-inserting non-MRI variables
combined_df_post_harm$dataset       <- reduced_df_combined$dataset
combined_df_post_harm$age           <- reduced_df_combined$age
combined_df_post_harm$yrs_education <- reduced_df_combined$yrs_education
combined_df_post_harm$mmse          <- reduced_df_combined$mmse
combined_df_post_harm <- combined_df_post_harm %>%
                         select(c("dataset","age","yrs_education","mmse"), everything())

# Adding column titles to harmonized df
names(combined_df_post_harm)[c(5:74)]    <- parcel_names_wm
names(combined_df_post_harm)[c(75:209)]  <- parcel_names_gm
names(combined_df_post_harm)[c(210:211)] <- c("lesion_volume_in_m_l","number_of_lesions")
combined_df_post_harm <- clean_names(combined_df_post_harm)

# # Adding means across parcels for the WM and GM parcels
# combined_df_post_harm <- mutate(combined_df_post_harm,
#                                 wm_discon_mean = rowMeans(select(combined_df_post_harm, 1:70), na.rm = TRUE))
# combined_df_post_harm <- mutate(combined_df_post_harm,
#                                 gm_lesions_mean = rowMeans(select(combined_df_post_harm, 71:205), na.rm = TRUE))

# # Testing that nothing went badly wrong
# plot(reduced_df$wm_discon_mean, test_df_post_harm$wm_discon_mean)
# cor.test(reduced_df$wm_discon_mean, test_df_post_harm$wm_discon_mean)

# Some ggplots (WM vs MMSE)
# ggplot(data=reduced_df,
#        aes(x=log(wm_discon_mean), y=mmse, col=dataset, rm.na=TRUE)) +
#        geom_point(alpha=0.5) + geom_smooth(method=lm) + facet_wrap(vars(dataset)) + theme_bw()
# ggplot(data=test_df_post_harm,
#        aes(x=log(wm_discon_mean), y=mmse, col=dataset, rm.na=TRUE)) +
#        geom_point(alpha=0.5) + geom_smooth(method=lm) + facet_wrap(vars(dataset)) + theme_bw()
# ggplot(data=test_df_post_harm2,
#        aes(x=log(wm_discon_mean), y=mmse, col=dataset, rm.na=TRUE)) +
#        geom_point(alpha=0.5) + geom_smooth(method=lm) + facet_wrap(vars(dataset)) + theme_bw()


# PLS-PM: Preparing the dfs (combined and for individual datasets) ----

# Combined df (using the data after harmonization)
df_plspm_combined <- combined_df_post_harm
df_plspm_combined <- filter(df_plspm_combined, !is.na(df_plspm_combined$age))
df_plspm_combined <- filter(df_plspm_combined, !is.na(df_plspm_combined$yrs_education))
df_plspm_combined <- filter(df_plspm_combined, !is.na(df_plspm_combined$mmse))

# Creating copy of the reduced df (before harmonization),
# with sex removed since this variable was for the harmonization process
reduced_df_combined_no_harm <- reduced_df_combined %>% select(-ptid,-sex)
reduced_df_combined_no_harm <- filter(reduced_df_combined_no_harm,
                                      !is.na(reduced_df_combined_no_harm$age))
reduced_df_combined_no_harm <- filter(reduced_df_combined_no_harm,
                                      !is.na(reduced_df_combined_no_harm$yrs_education))
reduced_df_combined_no_harm <- filter(reduced_df_combined_no_harm,
                                      !is.na(reduced_df_combined_no_harm$mmse))

# Creating reduced data frames for the individual datasets
df_plspm_adni   <- reduced_df_combined_no_harm[reduced_df_combined_no_harm$dataset == "ADNI",]
df_plspm_bial   <- reduced_df_combined_no_harm[reduced_df_combined_no_harm$dataset == "Bialystok",]
df_plspm_spreng <- reduced_df_combined_no_harm[reduced_df_combined_no_harm$dataset == "Spreng",]

# Removing columns with near-zero variance -- combined dataset
combined_nearzero_wm_parcels  <- nearZeroVar(df_plspm_combined[,5:74])
combined_nearzero_gm_parcels  <- nearZeroVar(df_plspm_combined[,75:209])
combined_nearzero_whole_brain <- nearZeroVar(df_plspm_combined[,210:211])
df_plspm_combined <- df_plspm_combined[,-c(combined_nearzero_wm_parcels,
                                           combined_nearzero_gm_parcels,
                                           combined_nearzero_whole_brain)]

# Removing columns with near-zero variance -- ADNI dataset
adni_nearzero_wm_parcels  <- nearZeroVar(df_plspm_adni[,5:74])
adni_nearzero_gm_parcels  <- nearZeroVar(df_plspm_adni[,75:209])
adni_nearzero_whole_brain <- nearZeroVar(df_plspm_adni[,210:211])
df_plspm_adni <- df_plspm_adni[,-c(adni_nearzero_wm_parcels,
                                   adni_nearzero_gm_parcels,
                                   adni_nearzero_whole_brain)]

# Removing columns with near-zero variance -- Bialystok dataset
bial_nearzero_wm_parcels  <- nearZeroVar(df_plspm_bial[,5:74])
bial_nearzero_gm_parcels  <- nearZeroVar(df_plspm_bial[,75:209])
bial_nearzero_whole_brain <- nearZeroVar(df_plspm_bial[,210:211])
df_plspm_bial <- df_plspm_bial[,-c(bial_nearzero_wm_parcels,
                                   bial_nearzero_gm_parcels,
                                   bial_nearzero_whole_brain)]

# Removing columns with near-zero variance -- Spreng dataset
spreng_nearzero_wm_parcels  <- nearZeroVar(df_plspm_spreng[,5:74])
spreng_nearzero_gm_parcels  <- nearZeroVar(df_plspm_spreng[,75:209])
spreng_nearzero_whole_brain <- nearZeroVar(df_plspm_spreng[,210:211])
df_plspm_spreng <- df_plspm_spreng[,-c(spreng_nearzero_wm_parcels,
                                       spreng_nearzero_gm_parcels,
                                       spreng_nearzero_whole_brain)]


# PLS-PM: Functions ----

analyze_plspm_model <- function(plspm_model)
{

# General summary of characteristics
summary(plspm_model)

# Check measures of unidimensionality
plspm_model$unidim
# plot(plspm_model, what = "loadings")

# Show the outer model
plspm_model$outer_model

# Barchart of loadings
loadings_barchart <- ggplot(data = plspm_model$outer_model, aes(x = name, y = loading, fill = block)) +
                     geom_bar(stat = "identity", position = "dodge") +
                     # Threshold line (based on the rule of thumb that acceptable loadings are above 0.7)
                     geom_hline(yintercept = 0.7, color = "black") +
                     # Other improvements to appearance
                     ggtitle("Barchart of Loadings") + theme_bw() + theme(plot.title=element_text(hjust=0.5)) +
                     # Rotate x-axis names
                     theme(axis.text.x = element_text(angle = 90))

# Check cross loadings for "traitors"
plspm_model$crossloadings

# Check the overall coefficients between the factors
plot(plspm_model)
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

# There may be a bit more that can be added from the Sanchez book

# Return the ggplot or else it won't print (even without <- )
# It also worked to have the ggplot as the last thing in the function
# A forum suggested using the print command around the ggplot but that didn't help
return(loadings_barchart)

}


# PLS-PM: Prepare matrices, blocks, and modes ----

# Model with 6 latent variables

# Path matrix
plspm_6_matrix <- matrix(c(0,0,0,0,0,0,    # age
                           0,0,0,0,0,0,    # lesions_wm_parcels
                           0,0,0,0,0,0,    # lesions_gm_parcels
                           0,0,0,0,0,0,    # lesions_whole_brain
                           1,1,1,1,0,0,    # lesions
                           1,0,0,0,1,0),   # cog_ability
                           nrow = 6, ncol = 6, byrow = TRUE)
colnames(plspm_6_matrix) <- c("age", "lesions_wm_parcels", "lesions_gm_parcels",
                              "lesions_whole_brain",  "lesions",  "cog_ability")
rownames(plspm_6_matrix) <- colnames(plspm_6_matrix)
# Modes
plspm_6_modes  <- rep("A",6)


# PLS-PM: Run models ----

# Combined dataset, 6 latent variables
plspm_6_blocks_combined <- list(2,
                                5:74-length(combined_nearzero_wm_parcels),
                                75-length(combined_nearzero_wm_parcels):209-sum(length(combined_nearzero_wm_parcels),
                                                                                length(combined_nearzero_gm_parcels)),
                                210-sum(length(combined_nearzero_wm_parcels),
                                        length(combined_nearzero_gm_parcels)):211-sum(length(combined_nearzero_wm_parcels),
                                                                                      length(combined_nearzero_gm_parcels),
                                                                                      length(combined_nearzero_whole_brain)),
                                5:211-sum(length(combined_nearzero_wm_parcels),
                                          length(combined_nearzero_gm_parcels),
                                          length(combined_nearzero_whole_brain)),
                                4)

stop()

plspm_model_combined <- plspm(df_plspm_combined, plspm_6_matrix,
                              plspm_6_blocks_combined, modes=plspm_6_modes,
                              boot.val = TRUE, br = 200)
analyze_plspm_model(plspm_model_combined)

stop()

# ADNI dataset, 6 latent variables

# Blocks
plspm_6_blocks <- list(209,1:70,71:205,206:207,1:207,211)
# Modes
plspm_6_modes  <- rep("A",6)

plspm_model_adni <- plspm(df_plspm_adni, plspm_6_matrix,
                          plspm_6_blocks, modes=plspm_6_modes,
                          boot.val = TRUE, br = 200)
analyze_plspm_model(plspm_model_adni)

stop()

# Bialystok dataset, 6 latent variables
plspm_model_bial <- plspm(df_plspm_bial, plspm_6_matrix,
                          plspm_6_blocks, modes=plspm_6_modes,
                          boot.val = TRUE, br = 200)
analyze_plspm_model(plspm_model_bial)

# Spreng dataset, 6 latent variables
plspm_model_spreng <- plspm(df_plspm_spreng, plspm_6_matrix,
                            plspm_6_blocks, modes=plspm_6_modes,
                            boot.val = TRUE, br = 200)
analyze_plspm_model(plspm_model_spreng)

stop()


# # PLS-PM: Preparing alternate df for Bialystok dataset ----
# 
# # Need to change df name if I use this
# 
# # Creating data frame
# df_plspm_bial <- merge(df_bial_wm_discon[,c(1,3:72)],
#                        df_bial_gm_lesions[,c(1,3:137)],
#                        by="ptid")
# df_plspm_bial <- merge(df_plspm_bial,
#                        df_bial_wm_lesions[,c(1,3:4)],
#                        by="ptid")
# df_plspm_bial <- merge(df_plspm_bial,
#                        df_bial_lsbq[,c("ptid","total_age")],
#                        by="ptid")
# df_plspm_bial <- merge(df_plspm_bial,
#                        df_bial_criq[,c("ptid","yrs_of_edu")],
#                        by="ptid")
# df_plspm_bial <- merge(df_plspm_bial,
#                        df_bial_mmse[,c("ptid","mmse_total_score")],
#                        by="ptid")
# 
# # Renaming columns to match the combined data frame used above
# names(df_plspm_bial)[names(df_plspm_bial) == "total_age"]        <- "age"
# names(df_plspm_bial)[names(df_plspm_bial) == "yrs_of_edu"]       <- "yrs_education"
# names(df_plspm_bial)[names(df_plspm_bial) == "mmse_total_score"] <- "mmse"
# 
# # # Removing ptid column
# # df_plspm_bial <- select(df_plspm_bial, -"ptid")
# # 
# # all_0_cols <- c()
# # for (i in 1:ncol(df_plspm_bial))
# # {
# #   if (!all(df_plspm_bial[,i] == 0) == FALSE)
# #     all_0_cols <- append(all_0_cols, i)
# # }
# 
# # Eliminating rows with missing values
# df_plspm_bial <- filter(df_plspm_bial, !is.na(df_plspm_bial$age))
# df_plspm_bial <- filter(df_plspm_bial, !is.na(df_plspm_bial$yrs_education))
# df_plspm_bial <- filter(df_plspm_bial, !is.na(df_plspm_bial$mmse))


# # PLS-PM: Model for Bialystok dataset ----
# 
# # # Rows for the path matrix
# # wm_discon_b            <- c(0,0,0,0,0,0,0,0)
# # gm_lesions_b           <- c(0,0,0,0,0,0,0,0)
# # lesions_tracts_b       <- c(1,1,0,0,0,0,0,0)
# # lesions_whole_brain_b  <- c(0,0,0,0,0,0,0,0)
# # age_b                  <- c(0,0,0,0,0,0,0,0)
# # cog_reserve_b          <- c(0,0,0,0,0,0,0,0)
# # lesions_b              <- c(0,0,1,1,1,1,0,0)
# # cog_ability_b          <- c(0,0,0,0,1,1,1,0)
# 
# # Rows
# cog_reserve_b <- c(0,0,0)
# lesions_b     <- c(1,0,0)
# cog_ability_b <- c(1,1,0)
# 
# # Create path matrix
# # path_matrix_b <- rbind(wm_discon_b, gm_lesions_b, lesions_tracts_b, lesions_whole_brain_b,
# #                        age_b, cog_reserve_b, lesions_b, cog_ability_b)
# path_matrix_b <- rbind(cog_reserve_b, lesions_b, cog_ability_b)
# colnames(path_matrix_b) <- rownames(path_matrix_b)
# 
# # Plot the path matrix
# innerplot(path_matrix_b, box.size = 0.1)
# 
# # Blocks and modes
# # plspm_blocks_b <- list(1:70,71:205,1:205,206:207,208,209,1:207,210)
# plspm_blocks_b <- list(209,207:208,210)
# plspm_modes_b  <- rep("A",3)
# 
# # Run PLS-PM
# plspm_model_bial <- plspm(df_plspm_bial, path_matrix_b, plspm_blocks_b, modes=plspm_modes_b,
#                           boot.val = TRUE, br = 200)
# 
# analyze_plspm_model(plspm_model_bial)


# # PLS-PM: Extra unused code ----
#
# # Plot the path matrix
# innerplot(plspm_matrix, box.size = 0.1)
#
# # Rows for the path matrix -- old version with cognitive reserve
# wm_discon            <- c(0,0,0,0,0,0,0,0)
# gm_lesions           <- c(0,0,0,0,0,0,0,0)
# lesions_tracts       <- c(1,1,0,0,0,0,0,0)
# lesions_whole_brain  <- c(0,0,0,0,0,0,0,0)
# age                  <- c(0,0,0,0,0,0,0,0)
# cog_reserve          <- c(0,0,0,0,0,0,0,0)
# lesions              <- c(0,0,1,1,1,1,0,0)
# cog_ability          <- c(0,0,0,0,1,1,1,0)
#
# # Simplest PLS-PM model code
# plspm_model_combined <- plspm(df_plspm_combined, path_matrix, plspm_blocks)


# # Correlations between cognitive measures and WM discon averaged across tracts ----
# 
# # This section needs to be updated to use the better organized data frames
# 
# # Raw correlations -- better if I could find a way to save output for all of these to a single data frame
#
# # With white matter disconnectivity
# cor(df_adni_wm$tracts_discon_mean, df_adni_wm$adni_mem, use='complete.obs')
# cor(df_adni_wm$tracts_discon_mean, df_adni_wm$adni_ef,  use='complete.obs')
# cor(df_adni_wm$tracts_discon_mean, df_adni_wm$adni_lan, use='complete.obs')
# cor(df_adni_wm$tracts_discon_mean, df_adni_wm$adni_vs,  use='complete.obs')
# cor(df_adni_wm$tracts_discon_mean, df_adni_wm$adni_ef2, use='complete.obs')
# 
# # ggplots of cog scores and white matter disconnectivity
# # Important note: consider adding the following code to the plots to make groups separate panels:
# # + facet_wrap(vars(dx))
# testing_this <- ggplot(data=df_adni_wm[df_adni_wm$dx!="na",], aes(x=tracts_discon_mean, y=adni_mem, col=dx, rm.na=TRUE)) +
#        geom_point() + geom_smooth(method=lm) + theme_bw()
# ggplot(data=df_adni_wm[df_adni_wm$dx!="na",], aes(x=tracts_discon_mean, y=adni_ef, col=dx, rm.na=TRUE)) +
#        geom_point() + geom_smooth(method=lm) + theme_bw()
# ggplot(data=df_adni_wm[df_adni_wm$dx!="na",], aes(x=tracts_discon_mean, y=adni_lan, col=dx, rm.na=TRUE)) +
#        geom_point() + geom_smooth(method=lm) + theme_bw()
# ggplot(data=df_adni_wm[df_adni_wm$dx!="na",], aes(x=tracts_discon_mean, y=adni_vs, col=dx, rm.na=TRUE)) +
#        geom_point() + geom_smooth(method=lm) + theme_bw()
# ggplot(data=df_adni_wm[df_adni_wm$dx!="na",], aes(x=tracts_discon_mean, y=adni_ef2, col=dx, rm.na=TRUE)) +
#        geom_point() + geom_smooth(method=lm) + theme_bw()
# 
# # Sample code for a ggplot for correlation between a cognitive domain and disconnectivity for a single tract
# ggplot(data=df_tracts_cog[df_tracts_cog$dx != 'NA',], aes(x=adni_ef, y=cs_l, col=dx, rm.na=TRUE)) +
#        geom_point() + geom_smooth(method=lm) + theme_bw() + facet_wrap(vars(dx))
# 
# 
# # Combined data frame with demographic, cognitive, and WM disconnectivity data
# df_adni_wm <- df_adni_adnimerge
# df_adni_wm <- left_join(df_adni_wm, df_adni_uwnpsychsum)
# df_adni_wm <- df_adni_wm %>% mutate_at('ptid', as.character)
# df_adni_wm <- left_join(df_adni_wm, df_tracts_discon)
# df_adni_wm <- df_adni_wm %>% mutate_at(c(116:120), as.numeric)
# df_adni_wm <- df_adni_wm %>% mutate_at('dx', as.factor)
# # Creating separate data frame with subset of columns (diagnosis, 5 cognitive totals, tract discon by parcel)
# df_adni_wm_reduced <- df_adni_wm[,c(61,116:190)]


# # PLS ----
# # Haven't done much adaptation of this part, still fairly similar to Dr. Anderson's code
#
# # Needs to be adapted to use the new better organized data frames
# 
# # Preparing data frames for PLS
# df_adni_wm_pls       <- df_adni_wm[,c(1,116:190)]
# df_adni_wm_pls       <- df_adni_wm_pls[complete.cases(df_adni_wm_pls),]
# df_adni_wm_pls_cog   <- df_adni_wm_pls[,c(1,2:6)]
# df_adni_wm_pls_tract <- df_adni_wm_pls[,c(1,7:76)]
# 
# # Running PLS
# pls <- tepPLS(df_adni_wm_pls_cog[,-1], df_adni_wm_pls_tract[,-1], graphs=TRUE)
# 
# # Permutation tests by columns
# # Note: For the perm.bycol function, the default type is byMat which permute by labels of observations, 
# # whereas the byColumns option permutes all cols of each data matrix independently
# my_nIter <- 1000   # Number of permutations for the permutation tests
# perm.bycol <- perm4PLSC(df_adni_wm_pls_cog[,-1], df_adni_wm_pls_tract[,-1],
#                         permType = 'byColumns', nIter = my_nIter)
# scree      <- PlotScree(ev = pls$TExPosition.Data$eigs,
#                         p.ev = perm.bycol$pEigenvalues,
#                         title = "Explained Variance per Dimension + Permutation Tests",
#                         plotKaiser = TRUE)
# 
# # Computing bootstrap ratios
# resBoot4PLSC <- Boot4PLSC(df_adni_wm_pls_cog[,-1],
#                           df_adni_wm_pls_tract[,-1],
#                           nIter = my_nIter,
#                       Fi = pls$TExPosition.Data$fi,
#                       Fj = pls$TExPosition.Data$fj,
#                       nf2keep = 3,
#                       critical.value = 2,
#                       # To be implemented later
#                       # has no effect currently
#                       alphaLevel = .05)
# 
# BR.I <- resBoot4PLSC$bootRatios.i %>%
#   as.data.frame() %>%
#   add_rownames(var = "rowname")
# BR.J <- resBoot4PLSC$bootRatios.j%>%
#   as.data.frame() %>%
#   add_rownames(var = "rowname")
# 
# cog_plot <- BR.I %>%
#   ggplot(aes(x = rowname, y = `Dimension 1`)) + geom_bar(stat="identity") +
#   geom_hline(yintercept = -2) + geom_hline(yintercept = 2) + coord_flip() + theme_bw()
# 
# tracts_plot <- BR.J %>%
#   ggplot(aes(x = rowname, y = `Dimension 1`)) + geom_bar(stat="identity") +
#   geom_hline(yintercept = -2) + geom_hline(yintercept = 2) + coord_flip() + theme_bw()
# 
# cog_plot
# tracts_plot


# # Find parcels in which lesions are associated with cognitive functions ----
# 
# # Function to create correlations by parcel
# create_df_cors <- function(my_df, non_brain_cols, cog_measure) {
#   out <- my_df %>%
#     gather(key = "region", value = "brain_value", -non_brain_cols) %>%
#     group_by(region) %>%
#     do(cor_test = cor.test(.[[cog_measure]], .[["brain_value"]])) %>%
#     rowwise(region) %>%
#     transmute(
#       region,
#       correlation_coefficient = cor_test$estimate,
#       p_value = cor_test$p.value,
#       confidence_lower = cor_test$conf.int[1],
#       confidence_upper = cor_test$conf.int[2]
#     )
# }
# 
# # Creating empty list to hold data frames with correlation results
# cor_dfs_list <- list()
# datasets <- c("ADNI", "Bialystok", "Spreng")
# 
# # Calculating correlation results for GM lesions
# cor_dfs_list[[1]] <- create_df_cors(df_combined_gm_lesions, c(1:7), 7)
# for (i in 1:3)
# {
#   cor_dfs_list[[i+1]] <-
#       create_df_cors(df_combined_gm_lesions[df_combined_gm_lesions$dataset == datasets[i],], c(1:7), 7)
# }
# 
# # Calculating correlation results for WM disconnectivity
# cor_dfs_list[[5]] <- create_df_cors(df_combined_wm_discon, c(1:7), 7)
# for (i in 1:3)
# {
#   cor_dfs_list[[i+5]] <-
#       create_df_cors(df_combined_wm_discon[df_combined_wm_discon$dataset == datasets[i],], c(1:7), 7)
# }
# 
# # Calculating correlation results for GM disconnectivity
# cor_dfs_list[[9]] <- create_df_cors(df_combined_gm_discon, c(1:7), 7)
# for (i in 1:3)
# {
#   cor_dfs_list[[i+9]] <-
#       create_df_cors(df_combined_gm_discon[df_combined_gm_discon$dataset == datasets[i],], c(1:7), 7)
# }
# 
# # Calculating correlation results for change in SSPL
# cor_dfs_list[[13]] <- create_df_cors(df_combined_delta_sspl, c(1:7), 7)
# for (i in 1:3)
# {
#   cor_dfs_list[[i+13]] <-
#       create_df_cors(df_combined_delta_sspl[df_combined_delta_sspl$dataset == datasets[i],], c(1:7), 7)
# }
# 
# # Sorting all data frames in correlations list by correlation coefficient
# # Order of LQT measures if GM lesions, WM discon, GM discon (direct), delta SSPL
# # Within each of those, order of datasets is All, ADNI, Bialystok, Spreng
# cor_dfs_sorted_list        <- cor_dfs_list
# for (i in 1:16)
# {
#   cor_dfs_sorted_list[[i]] <- cor_dfs_sorted_list[[i]][order(cor_dfs_sorted_list[[i]]$correlation_coefficient),]
# }
# 
# # Creating a new list based on the sorted list that only includes the top three parcels/connections per df
# cor_dfs_list_top5        <- cor_dfs_sorted_list
# for (i in 1:16)
# {
#   cor_dfs_list_top5[[i]] <- cor_dfs_list_top5[[i]][1:5,]
#   # Make the regions an ordered factor so that it doesn't alphabetize them for the plots
#   cor_dfs_list_top5[[i]]$region <- factor(cor_dfs_list_top5[[i]]$region, cor_dfs_list_top5[[i]]$region)
# }
# 
# # Creating figure equivalent to Panel 2 of Figure 7 in Griffis et al. (2021)
# cor_plot_titles <- c("GM Lesions: All Datasets", "GM Lesions: ADNI", "GM Lesions: Bialystok", "GM Lesions: Spreng",
#                      "WM Discon: All Datasets",  "WM Discon: ADNI",  "WM Discon: Bialystok",  "WM Discon: Spreng",
#                      "GM Discon: All Datasets",  "GM Discon: ADNI",  "GM Discon: Bialystok",  "GM Discon: Spreng",
#                      "Change in SSPL: All Datasets",  "Change in SSPL: ADNI",  "Change in SSPL: Bialystok",  "Change in SSPL: Spreng")
# top5_plots <- vector("list", 16)
# for (i in 1:16)
# {
#   top5_plots[[i]] <- ggplot(data=cor_dfs_list_top5[[i]],
#                             aes(x=correlation_coefficient, y=region)) +
#                      geom_bar(stat="identity", fill="blue", alpha=.6, width=.4) +
#                      scale_x_reverse() +
#                      xlab("Pearson's r") +
#                      ylab("Parcel") +
#                      ggtitle(cor_plot_titles[i]) +
#                      theme_bw()
# }
# 
# # Create a figure showing the top 5 GM parcels associated with MMSE score in each dataset
# figure_top_gm_lesion_parcels <- grid.arrange(top5_plots[[1]], top5_plots[[2]],
#                                              top5_plots[[3]], top5_plots[[4]],
#                                              nrow = 2)
# figure_top_gm_lesion_parcels
# 
# # Create a figure showing the top 5 GM parcels associated with MMSE score in each dataset
# figure_top_wm_discon_parcels <- grid.arrange(top5_plots[[5]], top5_plots[[6]],
#                                              top5_plots[[7]], top5_plots[[8]],
#                                              nrow = 2)
# figure_top_wm_discon_parcels
# 
# 
# # Create a data frame with correlations between the MMSE and GM parcels in different datasets
# df_gm_parcels_mmse_cors <- data.frame("All"       <- cor_dfs_list[[1]]$correlation_coefficient,
#                                       "ADNI"      <- cor_dfs_list[[2]]$correlation_coefficient,
#                                       "Bialystok" <- cor_dfs_list[[3]]$correlation_coefficient,
#                                       "Spreng"    <- cor_dfs_list[[4]]$correlation_coefficient)
# # Create a data frame with correlations between the MMSE and WM parcels in different datasets
# df_wm_discon_mmse_cors  <- data.frame("All"       <- cor_dfs_list[[5]]$correlation_coefficient,
#                                       "ADNI"      <- cor_dfs_list[[6]]$correlation_coefficient,
#                                       "Bialystok" <- cor_dfs_list[[7]]$correlation_coefficient,
#                                       "Spreng"    <- cor_dfs_list[[8]]$correlation_coefficient)
# 
# # Correlogram for correlations by parcel with MMSE, gray matter lesions
# ggpairs(df_gm_parcels_mmse_cors, title="GM/MMSE Corrs. by Dataset",
#         upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
#         lower = list(continuous = wrap("smooth", color='blue', alpha=0.25), combo = "facethist", discrete = "facetbar", na = "na"),
#         diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
#         columnLabels = c("All", "ADNI", "Bialystok", "Spreng"), axisLabels = "none") +
#         theme_bw() + theme(plot.title=element_text(hjust=0.5))
# 
# # Correlogram for correlations by parcel with MMSE, white matter disconnectivity
# ggpairs(df_wm_discon_mmse_cors,  title="WM/MMSE Corrs. by Dataset",
#         upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na"),
#         lower = list(continuous = wrap("smooth", color='blue', alpha=0.25), combo = "facethist", discrete = "facetbar", na = "na"),
#         diag = list(continuous = "blankDiag", discrete = "barDiag", na = "naDiag"),
#         columnLabels = c("All", "ADNI", "Bialystok", "Spreng"), axisLabels = "none") +
#         theme_bw() + theme(plot.title=element_text(hjust=0.5))


# # Working on displaying on a brain image ----
# 
# # Starting to set up ggseg images
# data(yeo7)
# ggplot() + geom_brain(atlas=ggsegYeo2011::yeo7)
# plot(yeo7) +
#   theme(legend.position = "bottom",
#         legend.text = element_text(size = 9)) +
#   guides(fill = guide_legend(ncol = 3))
# 
# # Extra stuff I haven't used / got working for ggseg
#
# library(ggsegExtra)
# yeo_repo <- ggseg_atlas_repos("yeo")
# ggseg_atlas_repos("yeo")
# plot(yeo7$data)
#
# # coldcuts package would let us plot results onto a brain image if we have a NIFTI file of the atlas


# Notes for future reference (probably out of date) ----

# Sources:
# Check Gibbons et al. (2012) for overall ADNI executive function measure
# Check Crane et al. (2012) for memory

# Note: Use FDR for MC


# Possible libraries for plotting brain maps:
# ggseg
# ggbrain
# Method in paper (voxel-wise)



