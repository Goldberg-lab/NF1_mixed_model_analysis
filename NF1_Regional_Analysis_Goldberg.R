# Set the working directory for your data files

setwd("E:/RGC_Survival_data")


# Install needed packages (if not done already).

install.packages(c("dplyr", "tidyr", "readr", "tidyverse", "lme4", "ggplot2"))

# Prep the packages necessisary for analysis. Install first if needed.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

# Check the current working directory
getwd()

# Import while preserving the original column names
dat_6DC <- read.csv(
  "Data_RGC_6D_C.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  strip.white = TRUE
)

# Repair columns with missing column names

# Inspect invalid column names
which(is.na(names(dat_6DC)) | names(dat_6DC) == "")
names(dat_6DC)[is.na(names(dat_6DC)) | names(dat_6DC) == ""]

# Repair blank and duplicate names while preserving valid names
names(dat_6DC) <- vctrs::vec_as_names(
  names(dat_6DC),
  repair = "unique"
)

names(dat_6DC)

#Verfiy the import

names(dat_6DC)
dim(dat_6DC)
head(dat_6DC)

# Check that regional columns are preserved

regional_columns <- c(
  "0,0", "0,1", "0,2", "0,3",
  "1,0", "1,1", "1,2", "1,3",
  "2,0", "2,1", "2,2", "2,3"
)

setdiff(regional_columns, names(dat_6DC))

# Define columns

emm_options(lmer.df = "satterthwaite")

ctrl_label <- "NF1 fl/fl, Cre:- Ctrl"
gli_label  <- "NF1 KO/fl, Cre:+, Gli"

regional_columns <- c(
  "0,0", "0,1", "0,2", "0,3",
  "1,0", "1,1", "1,2", "1,3",
  "2,0", "2,1", "2,2", "2,3"
)

# Confirm the regional columns exist, this should return character(0)

missing_regional_columns <- setdiff(
  regional_columns,
  names(dat_6DC)
)

missing_regional_columns


# Prepare the Goldberg lab data

goldberg_region_wide <- dat_6DC %>%
  mutate(
    ID_original = trimws(
      as.character(`ID, Ear Mark`)
    ),
    
    Raised_original = trimws(
      as.character(Raised)
    ),
    
    Genotype_original = trimws(
      as.character(Genotype)
    ),
    
    Age_weeks = parse_number(
      as.character(`Age at Death (weeks)`)
    ),
    
    Sex = case_when(
      str_to_lower(trimws(as.character(Sex))) %in%
        c("female", "f") ~ "female",
      
      str_to_lower(trimws(as.character(Sex))) %in%
        c("male", "m") ~ "male",
      
      TRUE ~ NA_character_
    ),
    
    Eye = case_when(
      str_to_upper(trimws(as.character(Eye))) %in%
        c("L", "LEFT") ~ "Left",
      
      str_to_upper(trimws(as.character(Eye))) %in%
        c("R", "RIGHT") ~ "Right",
      
      TRUE ~ NA_character_
    ),
    
    Genotype_group = case_when(
      Genotype_original == ctrl_label ~ "Control",
      Genotype_original == gli_label  ~ "Gli",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Convert the regional columns to numeric
  mutate(
    across(
      all_of(regional_columns),
      ~ parse_number(as.character(.x))
    )
  ) %>%
  
  filter(
    str_to_lower(Raised_original) == "goldberg",
    Genotype_original %in% c(ctrl_label, gli_label),
    !is.na(ID_original),
    ID_original != "",
    !is.na(Age_weeks),
    !is.na(Sex),
    !is.na(Eye),
    !is.na(Genotype_group)
  ) %>%
  
  mutate(
    Genotype_group = factor(
      Genotype_group,
      levels = c("Control", "Gli")
    ),
    
    Sex = factor(
      Sex,
      levels = c("female", "male")
    ),
    
    Eye = factor(
      Eye,
      levels = c("Left", "Right")
    ),
    
    Animal_ID = factor(ID_original),
    
    # Unique identifier for each eye nested within an animal
    Animal_Eye_ID = interaction(
      Animal_ID,
      Eye,
      sep = "_",
      drop = TRUE
    )
  ) %>%
  droplevels()

# Check the imported regional measurements

regional_missingness <- goldberg_region_wide %>%
  summarise(
    across(
      all_of(regional_columns),
      ~ sum(is.na(.x))
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Position",
    values_to = "N_missing"
  )

print(
  as.data.frame(regional_missingness),
  row.names = FALSE
)

# Convert the 12 positions to three regional retina regions: 
# Center, Middle, Peripheral

goldberg_region_long <- goldberg_region_wide %>%
  pivot_longer(
    cols = all_of(regional_columns),
    names_to = "Retinal_position",
    values_to = "Position_density"
  ) %>%
  mutate(
    Retinal_region = case_when(
      str_starts(Retinal_position, "0,") ~ "Center",
      str_starts(Retinal_position, "1,") ~ "Middle",
      str_starts(Retinal_position, "2,") ~ "Peripheral",
      TRUE ~ NA_character_
    ),
    
    Retinal_region = factor(
      Retinal_region,
      levels = c(
        "Center",
        "Middle",
        "Peripheral"
      )
    )
  ) %>%
  filter(!is.na(Retinal_region))

# Check the position and region assignment

goldberg_region_long %>%
  distinct(Retinal_position, Retinal_region) %>%
  arrange(Retinal_region, Retinal_position)

# Average the three retional density measurements per eye

goldberg_region_dat <- goldberg_region_long %>%
  group_by(
    Animal_ID,
    Animal_Eye_ID,
    Genotype_group,
    Sex,
    Age_weeks,
    Eye,
    Retinal_region
  ) %>%
  summarise(
    N_positions_available = sum(
      !is.na(Position_density)
    ),
    
    Regional_density = if (
      sum(!is.na(Position_density)) > 0
    ) {
      mean(
        Position_density,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  filter(!is.na(Regional_density)) %>%
  droplevels()

# Check the number of positions used

position_availability <- goldberg_region_dat %>%
  count(
    Retinal_region,
    N_positions_available,
    name = "N_regional_observations"
  ) %>%
  arrange(
    Retinal_region,
    N_positions_available
  )

print(
  as.data.frame(position_availability),
  row.names = FALSE
)

# Complete region analysis

goldberg_region_complete <- goldberg_region_dat %>%
  filter(N_positions_available == 4) %>%
  droplevels()

# Center the age

mean_goldberg_age <- mean(
  goldberg_region_complete$Age_weeks,
  na.rm = TRUE
)

goldberg_region_complete <- goldberg_region_complete %>%
  mutate(
    Age_centered = Age_weeks - mean_goldberg_age
  )

# Check sample sizes

animal_counts <- goldberg_region_complete %>%
  distinct(
    Animal_ID,
    Genotype_group,
    Sex,
    Age_weeks
  ) %>%
  count(
    Genotype_group,
    Sex,
    name = "N_animals"
  ) %>%
  arrange(
    Genotype_group,
    Sex
  )

print(
  as.data.frame(animal_counts),
  row.names = FALSE
)

# Count regional observations

regional_observation_counts <- goldberg_region_complete %>%
  count(
    Genotype_group,
    Sex,
    Eye,
    Retinal_region,
    name = "N_observations"
  ) %>%
  arrange(
    Genotype_group,
    Sex,
    Eye,
    Retinal_region
  )

print(
  as.data.frame(regional_observation_counts),
  row.names = FALSE
)

# Check how many animals contribute complete left and right eyes

eye_availability <- goldberg_region_complete %>%
  distinct(
    Animal_ID,
    Eye
  ) %>%
  count(
    Animal_ID,
    name = "N_eyes"
  ) %>%
  count(
    N_eyes,
    name = "N_animals"
  )

eye_availability

# Fit models to test genotype by retinal region, try Additive or interaction

# Additive model

region_additive_model_ML <- lmer(
  Regional_density ~
    Genotype_group +
    Retinal_region +
    Age_centered +
    Sex +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = FALSE
)

# Genotype by region interaction model

region_interaction_model_ML <- lmer(
  Regional_density ~
    Genotype_group * Retinal_region +
    Age_centered +
    Sex +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = FALSE
)

# Compare the models

genotype_by_region_test <- anova(
  region_additive_model_ML,
  region_interaction_model_ML
)

genotype_by_region_test

# Additive model had best fit
# Use the addititve model with REML

region_additive_model <- lmer(
  Regional_density ~
    Genotype_group +
    Retinal_region +
    Age_centered +
    Sex +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = TRUE
)

# Check the model

summary(region_additive_model)

isSingular(
  region_additive_model,
  tol = 1e-4
)

as.data.frame(
  VarCorr(region_additive_model)
)

# Type III tests for the additive model

region_additive_type3 <- anova(
  region_additive_model,
  type = 3,
  ddf = "Satterthwaite"
)

region_additive_type3

# Overall Glioma versus Control difference

library(emmeans)

overall_genotype_emmeans <- emmeans(
  region_additive_model,
  ~ Genotype_group,
  weights = "equal",
  lmer.df = "satterthwaite"
)

overall_genotype_means <- as.data.frame(
  summary(
    overall_genotype_emmeans,
    infer = c(TRUE, TRUE)
  )
)

overall_genotype_result <- as.data.frame(
  summary(
    pairs(
      overall_genotype_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

overall_genotype_means
overall_genotype_result

# Retinal region comparisions

region_emmeans <- emmeans(
  region_additive_model,
  ~ Retinal_region,
  lmer.df = "satterthwaite"
)

region_mean_results <- as.data.frame(
  summary(
    region_emmeans,
    infer = c(TRUE, TRUE)
  )
)

region_pairwise_results <- as.data.frame(
  summary(
    pairs(
      region_emmeans,
      adjust = "tukey"
    ),
    infer = c(TRUE, TRUE)
  )
)

region_mean_results
region_pairwise_results


# Keep regional genotype estimates as descriptive results

region_interaction_model_REML <- lmer(
  Regional_density ~
    Genotype_group * Retinal_region +
    Age_centered +
    Sex +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = TRUE
)

genotype_by_region_emmeans <- emmeans(
  region_interaction_model_REML,
  ~ Genotype_group | Retinal_region,
  lmer.df = "satterthwaite"
)

genotype_by_region_results <- as.data.frame(
  summary(
    pairs(
      genotype_by_region_emmeans,
      reverse = TRUE,
      adjust = "none"
    ),
    infer = c(TRUE, TRUE)
  )
) %>%
  mutate(
    p_value_Holm = p.adjust(
      p.value,
      method = "holm"
    )
  )

genotype_by_region_results

# Directly compare the male and female genotype effects

genotype_sex_emmeans_all <- emmeans(
  genotype_sex_interaction_REML,
  ~ Genotype_group * Sex,
  weights = "equal",
  lmer.df = "satterthwaite"
)

genotype_sex_emmeans_all

genotype_effect_sex_difference <- contrast(
  genotype_sex_emmeans_all,
  method = list(
    "Male genotype effect - Female genotype effect" =
      c(1, -1, -1, 1)
  )
)

genotype_effect_sex_difference_result <- as.data.frame(
  summary(
    genotype_effect_sex_difference,
    infer = c(TRUE, TRUE)
  )
)

genotype_effect_sex_difference_result

# Test if the female mice have regional difference

female_region_dat <- goldberg_region_complete %>%
  filter(Sex == "female") %>%
  droplevels()

# Common genotype effect across regions
female_region_additive_ML <- lmer(
  Regional_density ~
    Genotype_group +
    Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = female_region_dat,
  REML = FALSE
)

# Genotype effect allowed to vary by region
female_region_interaction_ML <- lmer(
  Regional_density ~
    Genotype_group * Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = female_region_dat,
  REML = FALSE
)

female_genotype_by_region_test <- anova(
  female_region_additive_ML,
  female_region_interaction_ML
)

female_genotype_by_region_test

# Create the female-only dataset

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

female_region_dat <- goldberg_region_complete %>%
  filter(
    tolower(as.character(Sex)) == "female"
  ) %>%
  filter(
    complete.cases(
      Regional_density,
      Genotype_group,
      Retinal_region,
      Age_centered,
      Eye,
      Animal_ID,
      Animal_Eye_ID
    )
  ) %>%
  droplevels()

# Obtain female region specific means without imposing equal differences

female_region_interaction_REML <- lmer(
  Regional_density ~
    Genotype_group * Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = female_region_dat,
  REML = TRUE
)

  # Calculate adjusted means

female_region_emmeans <- emmeans(
  female_region_interaction_REML,
  ~ Genotype_group | Retinal_region,
  lmer.df = "satterthwaite"
)

female_region_mean_results <- as.data.frame(
  summary(
    female_region_emmeans,
    infer = c(TRUE, TRUE)
  )
)

# Calculate the female Gli-Control differences separately in each region

female_region_contrast_results <- as.data.frame(
  summary(
    pairs(
      female_region_emmeans,
      reverse = TRUE,
      adjust = "none"
    ),
    infer = c(TRUE, TRUE)
  )
) %>%
  mutate(
    p_value_Holm = p.adjust(
      p.value,
      method = "holm"
    )
  )

female_region_mean_results
female_region_contrast_results

# Save the female only data

write.csv(
  female_region_mean_results,
  "6DC_Goldberg_female_adjusted_means_by_region.csv",
  row.names = FALSE
)

write.csv(
  female_region_contrast_results,
  "6DC_Goldberg_female_Gli_vs_Control_by_region.csv",
  row.names = FALSE
)

# Directly test regional loss in males

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

male_region_dat <- goldberg_region_complete %>%
  filter(Sex == "male") %>%
  droplevels()

# Fit a common-effect model and region specific model using maximum likelihood

male_region_additive_ML <- lmer(
  Regional_density ~
    Genotype_group +
    Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = male_region_dat,
  REML = FALSE
)

male_region_interaction_ML <- lmer(
  Regional_density ~
    Genotype_group * Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = male_region_dat,
  REML = FALSE
)

male_genotype_by_region_test <- anova(
  male_region_additive_ML,
  male_region_interaction_ML
)

male_genotype_by_region_test

# Number of mice by sex and genotype

sex_genotype_counts <- goldberg_region_complete %>%
  distinct(Animal_ID, Sex, Genotype_group) %>%
  count(
    Sex,
    Genotype_group,
    name = "N_mice"
  ) %>%
  arrange(Sex, Genotype_group)

sex_genotype_counts

# Save the results

write.csv(
  male_genotype_by_region_test,
  "6DC_Goldberg_regional_male_genotype_by_region_test.csv",
  row.names = FALSE
)

write.csv(
  regional_means_by_sex_results,
  "6DC_Goldberg_regional_adjusted_means_by_sex_and_genotype.csv",
  row.names = FALSE
)

region_additive_type3_results <-
  tibble::rownames_to_column(
    as.data.frame(region_additive_type3),
    var = "Model_term"
  )

genotype_by_region_test_results <-
  tibble::rownames_to_column(
    as.data.frame(genotype_by_region_test),
    var = "Model"
  )

write.csv(
  genotype_by_region_test_results,
  "6DC_Goldberg_genotype_by_region_interaction_test.csv",
  row.names = FALSE
)

write.csv(
  region_additive_type3_results,
  "6DC_Goldberg_additive_regional_model_Type3.csv",
  row.names = FALSE
)

write.csv(
  overall_genotype_means,
  "6DC_Goldberg_overall_genotype_adjusted_means.csv",
  row.names = FALSE
)

write.csv(
  overall_genotype_result,
  "6DC_Goldberg_overall_Gli_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  region_mean_results,
  "6DC_Goldberg_retinal_region_adjusted_means.csv",
  row.names = FALSE
)

write.csv(
  region_pairwise_results,
  "6DC_Goldberg_retinal_region_pairwise_comparisons.csv",
  row.names = FALSE
)

write.csv(
  genotype_by_region_results,
  "6DC_Goldberg_descriptive_Gli_vs_Control_by_region.csv",
  row.names = FALSE
)

saveRDS(
  region_additive_model,
  "6DC_Goldberg_final_additive_regional_model.rds"
)


# Forest plot of genotype effects

regional_difference_plot_dat <- genotype_by_region_results %>%
  mutate(
    Retinal_region = factor(
      Retinal_region,
      levels = c(
        "Peripheral",
        "Middle",
        "Center"
      )
    )
  )

# Create the plot

regional_difference_plot <- ggplot(
  regional_difference_plot_dat,
  aes(
    x = estimate,
    y = Retinal_region
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey45",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(
      xmin = lower.CL,
      xmax = upper.CL
    ),
    orientation = "y",
    width = 0.15,
    linewidth = 0.7,
    color = "#D55E00"
  ) +
  geom_point(
    size = 3,
    color = "#D55E00"
  ) +
  scale_x_continuous(
    breaks = seq(-350, 100, by = 50),
    minor_breaks = seq(-350, 100, by = 25)
  ) +
  labs(
    x = expression(
      "Glioma - Control density (cells/mm"^2*")"
    ),
    y = NULL,
    title = "Genotype-associated density difference",
    subtitle = "Estimates with 95% confidence intervals"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold")
  )

regional_difference_plot

# Regional heat map bar chart
# Prepare the eye level dots

library(dplyr)
library(ggplot2)
library(scales)

eye_region_age_dat <- goldberg_region_complete %>%
  mutate(
    Plot_group = case_when(
      Sex == "female" & Genotype_group == "Control" ~ "Female Control",
      Sex == "female" & Genotype_group == "Gli"     ~ "Female Glioma",
      Sex == "male"   & Genotype_group == "Control" ~ "Male Control",
      Sex == "male"   & Genotype_group == "Gli"     ~ "Male Glioma",
      TRUE ~ NA_character_
    ),
    
    Plot_group = factor(
      Plot_group,
      levels = c(
        "Female Control",
        "Female Glioma",
        "Male Control",
        "Male Glioma"
      )
    ),
    
    Retinal_region = factor(
      Retinal_region,
      levels = c("Center", "Middle", "Peripheral")
    )
  ) %>%
  filter(
    !is.na(Plot_group),
    !is.na(Regional_density),
    !is.na(Age_weeks),
    !is.na(Eye)
  ) %>%
  droplevels()

# Calculate the summaries at the animal level

animal_region_summary_dat <- eye_region_age_dat %>%
  group_by(
    Animal_ID,
    Genotype_group,
    Sex,
    Age_weeks,
    Plot_group,
    Retinal_region
  ) %>%
  summarise(
    Animal_region_density = mean(
      Regional_density,
      na.rm = TRUE
    ),
    N_eyes = n_distinct(Eye),
    .groups = "drop"
  )

# Calculate the 95% confidence interval

bar_summary <- animal_region_summary_dat %>%
  group_by(
    Retinal_region,
    Plot_group
  ) %>%
  summarise(
    N_animals = n(),
    
    Mean_density = mean(
      Animal_region_density,
      na.rm = TRUE
    ),
    
    SD_density = sd(
      Animal_region_density,
      na.rm = TRUE
    ),
    
    SE_density = SD_density / sqrt(N_animals),
    
    CI_multiplier = if_else(
      N_animals > 1,
      qt(0.975, df = N_animals - 1),
      NA_real_
    ),
    
    Lower_95_CI =
      Mean_density - CI_multiplier * SE_density,
    
    Upper_95_CI =
      Mean_density + CI_multiplier * SE_density,
    
    .groups = "drop"
  )

print(as.data.frame(bar_summary), row.names = FALSE)

# Create the graph

dodge_width <- 0.82

max_y <- ceiling(
  max(
    c(
      eye_region_age_dat$Regional_density,
      bar_summary$Upper_95_CI
    ),
    na.rm = TRUE
  ) / 500
) * 500

# Create the plot

regional_age_eye_plot <- ggplot() +
  
  # Bars: animal-level group means
  geom_col(
    data = bar_summary,
    aes(
      x = Retinal_region,
      y = Mean_density,
      fill = Plot_group,
      group = Plot_group
    ),
    position = position_dodge(width = dodge_width),
    width = 0.72,
    color = "black",
    linewidth = 0.35,
    alpha = 0.60
  ) +
  
  # Error bars: animal-level 95% confidence intervals
  geom_errorbar(
    data = bar_summary,
    aes(
      x = Retinal_region,
      ymin = Lower_95_CI,
      ymax = Upper_95_CI,
      group = Plot_group
    ),
    position = position_dodge(width = dodge_width),
    width = 0.15,
    linewidth = 0.65,
    color = "black"
  ) +
  
  # Individual left- and right-eye observations
  geom_point(
    data = eye_region_age_dat,
    aes(
      x = Retinal_region,
      y = Regional_density,
      color = Age_weeks,
      group = Plot_group
    ),
    position = position_jitterdodge(
      jitter.width = 0.10,
      jitter.height = 0,
      dodge.width = dodge_width,
      seed = 123
    ),
    size = 2.2,
    alpha = 0.85
  ) +
  
  # Bar colors represent sex and genotype
  scale_fill_manual(
    values = c(
      "Female Control" = "#56B4E9",
      "Female Glioma"  = "#D55E00",
      "Male Control"   = "#009E73",
      "Male Glioma"    = "#CC79A7"
    ),
    drop = FALSE
  ) +
  
  # Dot colors represent age
  scale_color_viridis_c(
    option = "plasma",
    direction = -1,
    name = "Age (weeks)"
  ) +
  
  scale_y_continuous(
    breaks = seq(0, max_y, by = 500),
    minor_breaks = seq(0, max_y, by = 100),
    labels = label_comma(),
    expand = expansion(mult = c(0, 0.03)),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  
  coord_cartesian(
    ylim = c(0, max_y)
  ) +
  
  labs(
    x = "Retinal region",
    y = expression(
      "Regional density (cells/mm"^2*")"
    ),
    fill = "Group",
    title = "Regional retinal density by genotype, sex, and age",
    subtitle = paste(
      "Each dot represents one eye;",
      "left and right eyes from the same mouse have the same age color"
    )
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

regional_age_eye_plot

# Save the graph


# Test wether the genotype effect differs by sex

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

# Model assuming the same genotype effect in males and females
genotype_sex_additive_ML <- lmer(
  Regional_density ~
    Genotype_group +
    Sex +
    Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = FALSE
)

# Model allowing the genotype effect to differ by sex
genotype_sex_interaction_ML <- lmer(
  Regional_density ~
    Genotype_group * Sex +
    Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = FALSE
)

genotype_by_sex_test <- anova(
  genotype_sex_additive_ML,
  genotype_sex_interaction_ML
)

genotype_by_sex_test


# Estimate Gli versus Control separately in males and females

genotype_sex_interaction_REML <- lmer(
  Regional_density ~
    Genotype_group * Sex +
    Retinal_region +
    Age_centered +
    Eye +
    (1 | Animal_ID) +
    (1 | Animal_Eye_ID),
  data = goldberg_region_complete,
  REML = TRUE
)

summary(genotype_sex_interaction_REML)

isSingular(
  genotype_sex_interaction_REML,
  tol = 1e-4
)

# Calculate adjusted means

genotype_by_sex_emmeans <- emmeans(
  genotype_sex_interaction_REML,
  ~ Genotype_group | Sex,
  weights = "equal",
  lmer.df = "satterthwaite"
)

genotype_by_sex_means <- as.data.frame(
  summary(
    genotype_by_sex_emmeans,
    infer = c(TRUE, TRUE)
  )
)

genotype_by_sex_means

# Calculate Glioma minus control seperately by sex:

genotype_by_sex_contrasts <- pairs(
  genotype_by_sex_emmeans,
  reverse = TRUE,
  adjust = "none"
)

genotype_by_sex_results <- as.data.frame(
  summary(
    genotype_by_sex_contrasts,
    infer = c(TRUE, TRUE)
  )
) %>%
  mutate(
    p_value_Holm = p.adjust(
      p.value,
      method = "holm"
    )
  )

genotype_by_sex_results


# Save sex adjusted data

write.csv(
  genotype_by_sex_means,
  "6DC_Goldberg_adjusted_genotype_means_by_sex.csv",
  row.names = FALSE
)

write.csv(
  genotype_by_sex_results,
  "6DC_Goldberg_Gli_vs_Control_by_sex.csv",
  row.names = FALSE
)

write.csv(
  sex_difference_in_genotype_effect_results,
  "6DC_Goldberg_genotype_effect_difference_by_sex.csv",
  row.names = FALSE
)

write.csv(
  sex_genotype_counts,
  "6DC_Goldberg_animal_counts_by_sex_genotype.csv",
  row.names = FALSE
)

write.csv(
  sex_genotype_age_summary,
  "6DC_Goldberg_age_summary_by_sex_genotype.csv",
  row.names = FALSE
)

write.csv(
  tibble::rownames_to_column(
    as.data.frame(genotype_by_sex_test),
    var = "Model"
  ),
  "6DC_Goldberg_genotype_by_sex_interaction_test.csv",
  row.names = FALSE
)

saveRDS(
  genotype_sex_interaction_REML,
  "6DC_Goldberg_genotype_by_sex_mixed_model.rds"
)

write.csv(
  genotype_effect_sex_difference_result,
  "6DC_Goldberg_male_minus_female_genotype_effect.csv",
  row.names = FALSE
)
