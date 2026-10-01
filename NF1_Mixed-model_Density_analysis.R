# Set the working directory for your data files

setwd("E:/RGC_Survival_data")


# Install needed packages (if not done already).

install.packages(c("dplyr", "tidyr", "readr", "tidyverse", "lme4", "ggplot2"))

# Prep the packages necessisary for analysis. Install first if needed.

library(dplyr)
library(tidyr)
library(lme4)
library(readr)
library(tidyverse)
library(lme4)
library(ggplot2)

# Import the data without altering column names
dat <- read.csv(
  "Data_RGC_6_A_B_E_F.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify the density column
density_column <- intersect(
  c("Density", "Density (cells/mm^2)"),
  names(dat)
)[1]

if (is.na(density_column)) {
  stop("A Density column was not found. Run names(dat) to inspect the column names.")
}

# Create a consistently named analysis variable
dat <- dat %>%
  mutate(
    ID = factor(ID),
    Eye = toupper(trimws(Eye)),
    Density_value = as.numeric(.data[[density_column]])
  ) %>%
  filter(
    !is.na(ID),
    Eye %in% c("L", "R"),
    !is.na(Density_value)
  )

# Check that each animal has one left-eye and one right-eye measurement

dat %>%
  count(ID, Eye) %>%
  arrange(ID, Eye)

# Identify IDs with exactly one measurement from each eye

paired_ids <- dat %>%
  count(ID, Eye) %>%
  filter(n == 1) %>%
  count(ID) %>%
  filter(n == 2) %>%
  pull(ID)

paired_dat <- dat %>%
  filter(ID %in% paired_ids)

length(paired_ids)  # Number of complete animal pairs

# Estimate between- and within-animal variablility

variance_model <- lmer(
  Density_value ~ Eye + (1 | ID),
  data = paired_dat,
  REML = TRUE
)

summary(variance_model)

# Extract the variance components

variance_results <- as.data.frame(VarCorr(variance_model))

between_animal_variance <- variance_results$vcov[
  variance_results$grp == "ID"
]

within_animal_variance <- variance_results$vcov[
  variance_results$grp == "Residual"
]

variance_ratio <- between_animal_variance / within_animal_variance

ICC <- between_animal_variance /
  (between_animal_variance + within_animal_variance)

data.frame(
  Between_animal_variance = between_animal_variance,
  Within_animal_variance = within_animal_variance,
  Between_to_within_ratio = variance_ratio,
  ICC = ICC
)

# Directly inspect left-right agreement

wide_dat <- paired_dat %>%
  select(ID, Eye, Density_value) %>%
  pivot_wider(
    names_from = Eye,
    values_from = Density_value
  ) %>%
  mutate(
    eye_difference = R - L,
    animal_mean = (R + L) / 2
  )

wide_dat

# Calculate the correlation and summarize the paired differences

cor.test(
  wide_dat$L,
  wide_dat$R,
  method = "pearson"
)

wide_dat %>%
  summarise(
    number_of_animals = n(),
    mean_left = mean(L),
    mean_right = mean(R),
    mean_difference_R_minus_L = mean(eye_difference),
    SD_of_differences = sd(eye_difference)
  )

# Visualize the paired dataset

ggplot(paired_dat, aes(x = Eye, y = Density_value, group = ID)) +
  geom_line(alpha = 0.4, color = "grey50") +
  geom_point(size = 2) +
  labs(
    x = "Eye",
    y = "Cell density",
    title = "Left and right eye density within each animal"
  ) +
  theme_classic()

# Preparing the plotting data for left-right comparision

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

eye_plot_dat <- dat %>%
  mutate(
    Density_value = parse_number(
      as.character(`Density`)
    ),
    
    Eye_standard = case_when(
      toupper(trimws(as.character(Eye))) %in% c("L", "LEFT")  ~ "Left",
      toupper(trimws(as.character(Eye))) %in% c("R", "RIGHT") ~ "Right",
      TRUE ~ NA_character_
    ),
    
    Raised = trimws(as.character(Raised)),
    Genotype = trimws(as.character(Genotype)),
    ID = trimws(as.character(ID))
  ) %>%
  filter(
    !is.na(ID),
    ID != "",
    !is.na(Eye_standard),
    !is.na(Density_value),
    !is.na(Raised),
    Raised != "",
    !is.na(Genotype),
    Genotype != ""
  ) %>%
  mutate(
    Eye_standard = factor(
      Eye_standard,
      levels = c("Left", "Right")
    ),
    
    Raised_Genotype = interaction(
      Raised,
      Genotype,
      sep = " | ",
      drop = TRUE
    ),
    
    # Ensures animals with the same ID in different groups are distinct
    Animal_ID = interaction(
      Raised,
      Genotype,
      ID,
      sep = "_",
      drop = TRUE
    )
  ) %>%
  group_by(
    Animal_ID,
    ID,
    Raised,
    Genotype,
    Raised_Genotype,
    Eye_standard
  ) %>%
  summarise(
    Density_value = mean(Density_value, na.rm = TRUE),
    .groups = "drop"
  )

# Check the groups

table(eye_plot_dat$Raised_Genotype)

# Paired left–right eye plot
# Each line connects the left and right eyes of one mouse.

library(ggplot2)
library(grid)

max_y <- ceiling(
  max(eye_plot_dat$Density_value, na.rm = TRUE) / 500
) * 500

major_y <- seq(0, max_y, by = 500)
minor_y <- setdiff(
  seq(0, max_y, by = 100),
  major_y
)

paired_eye_plot <- ggplot(
  eye_plot_dat,
  aes(
    x = Eye_standard,
    y = Density_value,
    group = Animal_ID,
    color = Raised_Genotype
  )
) +
  geom_line(
    linewidth = 0.6,
    alpha = 0.55
  ) +
  geom_point(
    size = 2.2,
    alpha = 0.85
  ) +
  scale_color_viridis_d(
    option = "turbo",
    end = 0.9
  ) +
  scale_y_continuous(
    limits = c(0, max_y),
    breaks = major_y,
    minor_breaks = minor_y,
    expand = expansion(mult = c(0, 0.02)),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  labs(
    x = "Eye",
    y = expression("Cell density (cells/mm"^2*")"),
    color = "Raised | Genotype",
    title = "Paired left- and right-eye cell density",
    subtitle = "Each line represents one mouse"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    
    # Minor ticks on the left side of the y-axis
    axis.minor.ticks.y.left = element_line(
      color = "black",
      linewidth = 0.4
    ),
    
    axis.minor.ticks.length.y = grid::unit(2, "mm")
  )

paired_eye_plot

# Save the graph

ggsave(
  filename = "Paired_left_right_eye_density.pdf",
  plot = paired_eye_plot,
  device = cairo_pdf,
  width = 10,
  height = 7,
  units = "in"
)

# Plot the right-minus-left difference for each mouse
# This version directly displays the paired difference.

eye_difference_dat <- eye_plot_dat %>%
  select(
    Animal_ID,
    ID,
    Raised,
    Genotype,
    Raised_Genotype,
    Eye_standard,
    Density_value
  ) %>%
  pivot_wider(
    names_from = Eye_standard,
    values_from = Density_value
  ) %>%
  filter(
    !is.na(Left),
    !is.na(Right)
  ) %>%
  mutate(
    Difference_R_minus_L = Right - Left
  )

# Create the plot

difference_plot <- ggplot(
  eye_difference_dat,
  aes(
    x = Raised_Genotype,
    y = Difference_R_minus_L,
    color = Raised_Genotype
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.3,
    alpha = 0.85
  ) +
  scale_color_viridis_d(
    option = "turbo",
    end = 0.9
  ) +
  labs(
    x = "Raised and genotype",
    y = expression("Right - Left density (cells/mm"^2*")"),
    color = "Raised | Genotype",
    title = "Within-animal right-minus-left density differences"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

difference_plot

# Paired statistical analysis

# Paired t-test
t.test(wide_dat$R, wide_dat$L, paired = TRUE)

# Nonparametric alternative
wilcox.test(wide_dat$R, wide_dat$L, paired = TRUE, exact = FALSE)

# Calculare the direction of the left-right variance

wide_dat <- wide_dat %>%
  filter(!is.na(R), !is.na(L)) %>%
  mutate(Difference_R_minus_L = R - L)

wide_dat %>%
  summarise(
    N_animals = n(),
    Mean_right = mean(R),
    Mean_left = mean(L),
    Mean_R_minus_L = mean(Difference_R_minus_L),
    Median_R_minus_L = median(Difference_R_minus_L),
    N_right_higher = sum(Difference_R_minus_L > 0),
    N_left_higher = sum(Difference_R_minus_L < 0),
    N_ties = sum(Difference_R_minus_L == 0)
  )

# Positive median difference: right-eye density tends to be higher.
# Negative median difference: left-eye density tends to be higher.

# Estimate the paired location shift and confidence interval

wilcox_result <- wilcox.test(
  wide_dat$R,
  wide_dat$L,
  paired = TRUE,
  alternative = "two.sided",
  exact = FALSE,
  correct = TRUE,
  conf.int = TRUE,
  conf.level = 0.95
)

wilcox_result

# Extract the results

wilcox_summary <- data.frame(
  N_animals = nrow(wide_dat),
  V = unname(wilcox_result$statistic),
  P_value = wilcox_result$p.value,
  Estimated_location_shift = unname(wilcox_result$estimate),
  Lower_95_CI = wilcox_result$conf.int[1],
  Upper_95_CI = wilcox_result$conf.int[2],
  Median_observed_R_minus_L =
    median(wide_dat$Difference_R_minus_L)
)

wilcox_summary

# Export the results of the wilcox test

write.csv(
  wilcox_summary,
  "Overall_paired_Wilcoxon_right_vs_left.csv",
  row.names = FALSE
)

# Split data set into age and sex groups

# Prepare the data

library(dplyr)
library(tidyr)
library(lme4)
library(readr)

# Inspect the data table labels

library(dplyr)

sort(unique(trimws(as.character(dat$Raised))))
sort(unique(trimws(as.character(dat$Genotype))))

# Prepare the complete analysis dataset

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)

ctrl_label <- "NF1 fl/fl, Cre:- Ctrl"
gli_label  <- "NF1 KO/fl, Cre:+, Gli"
charles_genotype_label <- "C57Blk/6"

analysis_dat <- dat %>%
  mutate(
    Density_value = parse_number(
      as.character(`Density`)
    ),
    
    Age_weeks = parse_number(
      as.character(`Age at Death (weeks)`)
    ),
    
    ID = trimws(as.character(ID)),
    
    Gender = tolower(
      trimws(as.character(Gender))
    ),
    
    Raised_original = trimws(
      as.character(Raised)
    ),
    
    Genotype_original = trimws(
      as.character(Genotype)
    ),
    
    Eye = case_when(
      toupper(trimws(as.character(Eye))) %in%
        c("L", "LEFT") ~ "Left",
      
      toupper(trimws(as.character(Eye))) %in%
        c("R", "RIGHT") ~ "Right",
      
      TRUE ~ NA_character_
    ),
    
    # Standardize the raised-location names
    Raised_group = case_when(
      str_detect(
        str_to_lower(Raised_original),
        "goldberg"
      ) ~ "Goldberg",
      
      str_detect(
        str_to_lower(Raised_original),
        "gutmann"
      ) ~ "Gutmann",
      
      str_detect(
        str_to_lower(Raised_original),
        "charles"
      ) ~ "Charles River",
      
      TRUE ~ NA_character_
    ),
    
    # Standardize genotype/treatment labels
    Genotype_group = case_when(
      Genotype_original == ctrl_label ~ "Control",
      Genotype_original == gli_label ~ "Gli",
      
      # Charles River C57Blk/6 animals serve as controls
      Raised_group == "Charles River" &
        Genotype_original == charles_genotype_label ~ "Control",
      
      # This also handles Charles River animals if the genotype
      # label contains C57 but is not written exactly as above
      Raised_group == "Charles River" &
        str_detect(
          str_to_lower(Genotype_original),
          "c57"
        ) ~ "Control",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(ID),
    ID != "",
    !is.na(Density_value),
    !is.na(Age_weeks),
    !is.na(Eye),
    !is.na(Gender),
    !is.na(Raised_group),
    !is.na(Genotype_group)
  )

# Create source-specific age categories

analysis_dat <- analysis_dat %>%
  mutate(
    Age_group = case_when(
      # Goldberg females
      Raised_group == "Goldberg" &
        Gender == "female" &
        Age_weeks >= 11.5 &
        Age_weeks < 19 ~ "12-18",
      
      Raised_group == "Goldberg" &
        Gender == "female" &
        Age_weeks >= 23.5 &
        Age_weeks < 33 ~ "24-32",
      
      Raised_group == "Goldberg" &
        Gender == "female" &
        Age_weeks >= 51.5 &
        Age_weeks < 60 ~ "52-59",
      
      # Goldberg males
      Raised_group == "Goldberg" &
        Gender == "male" &
        Age_weeks >= 11.5 &
        Age_weeks < 19 ~ "12-18",
      
      Raised_group == "Goldberg" &
        Gender == "male" &
        Age_weeks >= 24.5 &
        Age_weeks < 39 ~ "25-38",
      
      # Gutmann groups
      Raised_group == "Gutmann" &
        between(Age_weeks, 16, 18) ~ "17",
      
      Raised_group == "Gutmann" &
        between(Age_weeks, 23, 25) ~ "24",
      
      # Charles River 24-week controls
      Raised_group == "Charles River" &
        between(Age_weeks, 23, 25) ~ "24",
      
      TRUE ~ NA_character_
    )
  )

# Create the comparision categories

analysis_dat <- analysis_dat %>%
  mutate(
    Comparison_group = case_when(
      # Goldberg comparisons
      Raised_group == "Goldberg" &
        Gender == "female" &
        Age_group == "12-18" ~
        "Goldberg_Female_12-18",
      
      Raised_group == "Goldberg" &
        Gender == "male" &
        Age_group == "12-18" ~
        "Goldberg_Male_12-18",
      
      Raised_group == "Goldberg" &
        Gender == "female" &
        Age_group == "24-32" ~
        "Goldberg_Female_24-32",
      
      Raised_group == "Goldberg" &
        Gender == "male" &
        Age_group == "25-38" ~
        "Goldberg_Male_25-38",
      
      Raised_group == "Goldberg" &
        Gender == "female" &
        Age_group == "52-59" ~
        "Goldberg_Female_52-59",
      
      # Gutmann 17-week comparisons, separated by gender
      Raised_group == "Gutmann" &
        Gender == "female" &
        Age_group == "17" ~
        "Gutmann_Female_17",
      
      Raised_group == "Gutmann" &
        Gender == "male" &
        Age_group == "17" ~
        "Gutmann_Male_17",
      
      # Gutmann Gli versus Charles River Control at 24 weeks
      Raised_group %in% c("Gutmann", "Charles River") &
        Gender == "female" &
        Age_group == "24" ~
        "GutmannGli_vs_CharlesControl_Female_24",
      
      Raised_group %in% c("Gutmann", "Charles River") &
        Gender == "male" &
        Age_group == "24" ~
        "GutmannGli_vs_CharlesControl_Male_24",
      
      TRUE ~ NA_character_
    ),
    
    # Identify which observations are allowed in each comparison
    Include_in_comparison = case_when(
      # Goldberg Control and Gli
      Raised_group == "Goldberg" &
        Genotype_group %in% c("Control", "Gli") ~ TRUE,
      
      # Gutmann Control and Gli at 17 weeks
      Raised_group == "Gutmann" &
        Age_group == "17" &
        Genotype_group %in% c("Control", "Gli") ~ TRUE,
      
      # At 24 weeks, use Charles River controls
      Raised_group == "Charles River" &
        Age_group == "24" &
        Genotype_group == "Control" ~ TRUE,
      
      # At 24 weeks, use Gutmann Gli mice
      Raised_group == "Gutmann" &
        Age_group == "24" &
        Genotype_group == "Gli" ~ TRUE,
      
      # Gutmann 24-week controls are intentionally omitted
      TRUE ~ FALSE
    )
  ) %>%
  filter(
    Include_in_comparison,
    !is.na(Comparison_group)
  ) %>%
  mutate(
    Genotype_group = factor(
      Genotype_group,
      levels = c("Control", "Gli")
    ),
    
    Eye = factor(
      Eye,
      levels = c("Left", "Right")
    ),
    
    Gender = factor(
      Gender,
      levels = c("female", "male")
    ),
    
    Raised_group = factor(
      Raised_group,
      levels = c(
        "Goldberg",
        "Gutmann",
        "Charles River"
      )
    ),
    
    Comparison_group = factor(
      Comparison_group,
      levels = c(
        "Goldberg_Female_12-18",
        "Goldberg_Male_12-18",
        "Goldberg_Female_24-32",
        "Goldberg_Male_25-38",
        "Goldberg_Female_52-59",
        "Gutmann_Female_17",
        "Gutmann_Male_17",
        "GutmannGli_vs_CharlesControl_Female_24",
        "GutmannGli_vs_CharlesControl_Male_24"
      )
    ),
    
    # Protects against repeated ID labels across raised locations
    Animal_ID = interaction(
      Raised_group,
      ID,
      sep = "_",
      drop = TRUE
    )
  ) %>%
  droplevels()

# Verfiy age assignments

age_assignment_check <- analysis_dat %>%
  distinct(
    Animal_ID,
    Raised_group,
    Genotype_group,
    Gender,
    Age_weeks,
    Age_group,
    Comparison_group
  ) %>%
  arrange(
    Raised_group,
    Comparison_group,
    Genotype_group,
    Age_weeks
  )

print(
  as.data.frame(age_assignment_check),
  row.names = FALSE
)

# Summarize the actual age range included in each category

age_range_check <- analysis_dat %>%
  distinct(
    Animal_ID,
    Comparison_group,
    Age_weeks
  ) %>%
  group_by(Comparison_group) %>%
  summarise(
    Minimum_age = min(Age_weeks),
    Maximum_age = max(Age_weeks),
    N_animals = n_distinct(Animal_ID),
    .groups = "drop"
  )

print(
  as.data.frame(age_range_check),
  row.names = FALSE
)


# Check samples sizes, count independant animals

animal_counts <- analysis_dat %>%
  distinct(
    Animal_ID,
    Comparison_group,
    Raised_group,
    Genotype_group,
    Gender
  ) %>%
  count(
    Comparison_group,
    Raised_group,
    Genotype_group,
    name = "N_animals"
  ) %>%
  arrange(
    Comparison_group,
    Genotype_group
  )

print(
  as.data.frame(animal_counts),
  row.names = FALSE
)

# Count eye observations

eye_counts <- analysis_dat %>%
  count(
    Comparison_group,
    Raised_group,
    Genotype_group,
    Eye,
    name = "N_eye_observations"
  ) %>%
  arrange(
    Comparison_group,
    Genotype_group,
    Eye
  )

print(
  as.data.frame(eye_counts),
  row.names = FALSE
)

# Check which comparision categories contain both groups

comparison_eligibility <- analysis_dat %>%
  distinct(
    Animal_ID,
    Comparison_group,
    Genotype_group
  ) %>%
  count(
    Comparison_group,
    Genotype_group,
    name = "N_animals"
  ) %>%
  group_by(Comparison_group) %>%
  summarise(
    N_control = sum(
      N_animals[Genotype_group == "Control"],
      na.rm = TRUE
    ),
    N_gli = sum(
      N_animals[Genotype_group == "Gli"],
      na.rm = TRUE
    ),
    Can_compare = N_control > 0 & N_gli > 0,
    .groups = "drop"
  )

print(
  as.data.frame(comparison_eligibility),
  row.names = FALSE
)

# Keep only estimable comparison categories

valid_comparison_groups <- comparison_eligibility %>%
  filter(Can_compare) %>%
  pull(Comparison_group)

model_dat <- analysis_dat %>%
  filter(
    Comparison_group %in% valid_comparison_groups
  ) %>%
  droplevels()

# Inspect any categories that were omitted:

comparison_eligibility %>%
filter(!Can_compare)

# Check dublicate eye measurements

duplicate_eye_check <- model_dat %>%
  count(
    Animal_ID,
    Comparison_group,
    Eye,
    name = "N_records"
  ) %>%
  filter(N_records > 1)

duplicate_eye_check
# IF this returns a 0 proceed directly, if not, average the repeated measurement

# Fit the mixed model

combined_mixed_model <- lmer(
  Density_value ~
    Genotype_group * Comparison_group +
    Eye * Comparison_group +
    (1 | Animal_ID),
  data = model_dat,
  REML = TRUE
)

summary(combined_mixed_model)

# Check the random effect

isSingular(
  combined_mixed_model,
  tol = 1e-4
)

as.data.frame(
  VarCorr(combined_mixed_model)
)

# Overall model tests

type3_tests <- anova(
  combined_mixed_model,
  type = 3,
  ddf = "Satterthwaite"
)

type3_tests

# Estimate Control-versus-Gli differences

emm_options(
  lmer.df = "satterthwaite"
)

# Group the Goldberg data

goldberg_model_dat <- model_dat %>%
  filter(
    Raised_group == "Goldberg",
    Genotype_group %in% c("Control", "Gli")
  ) %>%
  droplevels()



# Group the Gutmann 17-week data

gutmann17_model_dat <- model_dat %>%
  filter(
    Raised_group == "Gutmann",
    Age_group == "17",
    Genotype_group %in% c("Control", "Gli")
  ) %>%
  droplevels()

# Gutmann Gli versus Charles River Control at 24 weeks

week24_model_dat <- model_dat %>%
  filter(
    Age_group == "24",
    (
      Raised_group == "Gutmann" &
        Genotype_group == "Gli"
    ) |
      (
        Raised_group == "Charles River" &
          Genotype_group == "Control"
      )
  ) %>%
  mutate(
    Cohort = case_when(
      Raised_group == "Charles River" &
        Genotype_group == "Control" ~
        "CharlesRiver_Control",
      
      Raised_group == "Gutmann" &
        Genotype_group == "Gli" ~
        "Gutmann_Gli",
      
      TRUE ~ NA_character_
    ),
    Cohort = factor(
      Cohort,
      levels = c(
        "CharlesRiver_Control",
        "Gutmann_Gli"
      )
    )
  ) %>%
  filter(!is.na(Cohort)) %>%
  droplevels()

# Verfiy sample sizes seperately

# Goldberg counts

goldberg_counts <- goldberg_model_dat %>%
  distinct(
    Animal_ID,
    Comparison_group,
    Genotype_group
  ) %>%
  count(
    Comparison_group,
    Genotype_group,
    name = "N_animals"
  ) %>%
  arrange(
    Comparison_group,
    Genotype_group
  )

print(
  as.data.frame(goldberg_counts),
  row.names = FALSE
)

# Gutmann 17-week counts

gutmann17_counts <- gutmann17_model_dat %>%
  distinct(
    Animal_ID,
    Gender,
    Genotype_group
  ) %>%
  count(
    Gender,
    Genotype_group,
    name = "N_animals"
  ) %>%
  arrange(
    Gender,
    Genotype_group
  )

print(
  as.data.frame(gutmann17_counts),
  row.names = FALSE
)

# 24-week cohort counts

week24_counts <- week24_model_dat %>%
  distinct(
    Animal_ID,
    Gender,
    Cohort
  ) %>%
  count(
    Gender,
    Cohort,
    name = "N_animals"
  ) %>%
  arrange(
    Gender,
    Cohort
  )

print(
  as.data.frame(week24_counts),
  row.names = FALSE
)

# Check eye counts

goldberg_eye_counts <- goldberg_model_dat %>%
  count(
    Comparison_group,
    Genotype_group,
    Eye,
    name = "N_eye_observations"
  )

gutmann17_eye_counts <- gutmann17_model_dat %>%
  count(
    Gender,
    Genotype_group,
    Eye,
    name = "N_eye_observations"
  )

week24_eye_counts <- week24_model_dat %>%
  count(
    Gender,
    Cohort,
    Eye,
    name = "N_eye_observations"
  )


# Goldberg only mixed model

goldberg_model <- lmer(
  Density_value ~
    Genotype_group * Comparison_group +
    Eye * Comparison_group +
    (1 | Animal_ID),
  data = goldberg_model_dat,
  REML = TRUE
)

summary(goldberg_model)

isSingular(
  goldberg_model,
  tol = 1e-4
)

as.data.frame(
  VarCorr(goldberg_model)
)

# Goldberg Type III tests

goldberg_type3 <- anova(
  goldberg_model,
  type = 3,
  ddf = "Satterthwaite"
)

goldberg_type3

# Goldberg Control-versus-Gli comparisons

goldberg_genotype_emmeans <- emmeans(
  goldberg_model,
  ~ Genotype_group | Comparison_group,
  lmer.df = "satterthwaite"
)

goldberg_genotype_contrasts <- pairs(
  goldberg_genotype_emmeans,
  reverse = TRUE,
  adjust = "none"
)

goldberg_genotype_means <- as.data.frame(
  summary(
    goldberg_genotype_emmeans,
    infer = c(TRUE, TRUE)
  )
)

goldberg_genotype_results <- as.data.frame(
  summary(
    goldberg_genotype_contrasts,
    infer = c(TRUE, TRUE)
  )
) %>%
  mutate(
    p_value_Holm = p.adjust(
      p.value,
      method = "holm"
    )
  )

goldberg_genotype_results

# Calculating the % decrease in the Control to Glioma comparison

library(dplyr)
library(tidyr)
library(emmeans)

group_emmeans <- emmeans(
  goldberg_model,
  ~ Genotype_group | Comparison_group,
  lmer.df = "satterthwaite"
)

group_mean_table <- as.data.frame(group_emmeans)

# Reshape the adjusted means so that Ctrl and Glioma are seperate columns

percent_decrease_results <- group_mean_table %>%
  select(
    Comparison_group,
    Genotype_group,
    emmean
  ) %>%
  pivot_wider(
    names_from = Genotype_group,
    values_from = emmean
  ) %>%
  mutate(
    Difference_Gli_minus_Control = Gli - Control,
    
    Percent_change_Gli_vs_Control =
      100 * (Gli - Control) / Control,
    
    Percent_decrease_in_Gli =
      100 * (Control - Gli) / Control
  )

percent_decrease_results

# Join the percentages to the existing contrast table

control_means <- group_mean_table %>%
  filter(Genotype_group == "Control") %>%
  select(
    Comparison_group,
    Control_adjusted_mean = emmean
  )

contrast_with_percent <- goldberg_genotype_results %>%
  left_join(
    control_means,
    by = "Comparison_group"
  ) %>%
  mutate(
    Percent_decrease_in_Gli =
      -100 * estimate / Control_adjusted_mean
  )

contrast_with_percent

# Goldberg right-minus-left comparisons

goldberg_eye_emmeans <- emmeans(
  goldberg_model,
  ~ Eye | Comparison_group,
  lmer.df = "satterthwaite"
)

goldberg_eye_contrasts <- pairs(
  goldberg_eye_emmeans,
  reverse = TRUE,
  adjust = "none"
)

goldberg_eye_results <- as.data.frame(
  summary(
    goldberg_eye_contrasts,
    infer = c(TRUE, TRUE)
  )
) %>%
  mutate(
    p_value_Holm = p.adjust(
      p.value,
      method = "holm"
    )
  )

goldberg_eye_results

# Gutmann-only 17-week mixed model, with only one sex present

gutmann17_model <- lmer(
  Density_value ~
    Genotype_group +
    Eye +
    (1 | Animal_ID),
  data = gutmann17_model_dat,
  REML = TRUE
)

summary(gutmann17_model)

gutmann17_type3 <- anova(
  gutmann17_model,
  type = 3,
  ddf = "Satterthwaite"
)

gutmann17_genotype_emmeans <- emmeans(
  gutmann17_model,
  ~ Genotype_group,
  lmer.df = "satterthwaite"
)

gutmann17_genotype_means <- as.data.frame(
  summary(
    gutmann17_genotype_emmeans,
    infer = c(TRUE, TRUE)
  )
)

gutmann17_genotype_results <- as.data.frame(
  summary(
    pairs(
      gutmann17_genotype_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

gutmann17_eye_emmeans <- emmeans(
  gutmann17_model,
  ~ Eye,
  lmer.df = "satterthwaite"
)

gutmann17_eye_results <- as.data.frame(
  summary(
    pairs(
      gutmann17_eye_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

# Separate 24-week cohort analysis, one sex present

week24_model <- lmer(
  Density_value ~
    Cohort +
    Eye +
    (1 | Animal_ID),
  data = week24_model_dat,
  REML = TRUE
)

summary(week24_model)

week24_type3 <- anova(
  week24_model,
  type = 3,
  ddf = "Satterthwaite"
)

week24_cohort_emmeans <- emmeans(
  week24_model,
  ~ Cohort,
  lmer.df = "satterthwaite"
)

week24_cohort_means <- as.data.frame(
  summary(
    week24_cohort_emmeans,
    infer = c(TRUE, TRUE)
  )
)

week24_cohort_results <- as.data.frame(
  summary(
    pairs(
      week24_cohort_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

week24_eye_emmeans <- emmeans(
  week24_model,
  ~ Eye,
  lmer.df = "satterthwaite"
)

week24_eye_results <- as.data.frame(
  summary(
    pairs(
      week24_eye_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)


# Optional overall comparisons when interactions are nonsignificant

goldberg_additive_model <- lmer(
  Density_value ~
    Genotype_group +
    Comparison_group +
    Eye +
    (1 | Animal_ID),
  data = goldberg_model_dat,
  REML = TRUE
)

goldberg_overall_genotype_emmeans <- emmeans(
  goldberg_additive_model,
  ~ Genotype_group,
  weights = "equal",
  lmer.df = "satterthwaite"
)

goldberg_overall_genotype_results <- as.data.frame(
  summary(
    pairs(
      goldberg_overall_genotype_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

goldberg_overall_genotype_results

# Convert Type III tables for saving

goldberg_type3_results <- tibble::rownames_to_column(
  as.data.frame(goldberg_type3),
  var = "Model_term"
)

gutmann17_type3_results <- tibble::rownames_to_column(
  as.data.frame(gutmann17_type3),
  var = "Model_term"
)

week24_type3_results <- tibble::rownames_to_column(
  as.data.frame(week24_type3),
  var = "Model_term"
)

# Calculate % change in gutmann mice

library(emmeans)
library(dplyr)

gutmann17_overall_emm <- emmeans(
  gutmann17_model,
  ~ Genotype_group,
  weights = "equal",
  lmer.df = "satterthwaite"
)

gutmann17_overall_emm

emm_table <- as.data.frame(gutmann17_overall_emm)
V <- vcov(gutmann17_overall_emm)

control_i <- which(
  emm_table$Genotype_group == "Control"
)

gli_i <- which(
  emm_table$Genotype_group == "Gli"
)

control_mean <- emm_table$emmean[control_i]
gli_mean <- emm_table$emmean[gli_i]

# Percentage increase in Gli relative to Control
percent_increase <- 100 *
  (gli_mean - control_mean) / control_mean

# Delta-method derivatives
derivative_control <- -100 * gli_mean / control_mean^2
derivative_gli <- 100 / control_mean

gradient <- c(
  derivative_control,
  derivative_gli
)

V_sub <- V[
  c(control_i, gli_i),
  c(control_i, gli_i),
  drop = FALSE
]

percent_variance <- as.numeric(
  t(gradient) %*% V_sub %*% gradient
)

percent_SE <- sqrt(percent_variance)

gutmann17_percent_result <- data.frame(
  Control_adjusted_mean = control_mean,
  Gli_adjusted_mean = gli_mean,
  Difference_Gli_minus_Control =
    gli_mean - control_mean,
  Percent_increase_in_Gli = percent_increase,
  Percent_SE = percent_SE,
  Lower_1_SE = percent_increase - percent_SE,
  Upper_1_SE = percent_increase + percent_SE,
  Lower_95_CI_approx =
    percent_increase - 1.96 * percent_SE,
  Upper_95_CI_approx =
    percent_increase + 1.96 * percent_SE
)

gutmann17_percent_result

# Save the seperate results for Goldberg and Gutmann


write.csv(
  goldberg_counts,
  "Goldberg_animal_counts.csv",
  row.names = FALSE
)

write.csv(
  goldberg_genotype_means,
  "Goldberg_adjusted_genotype_means.csv",
  row.names = FALSE
)

write.csv(
  goldberg_genotype_results,
  "Goldberg_Gli_vs_Control_by_age_gender.csv",
  row.names = FALSE
)

write.csv(
  contrast_with_percent,
  file = "Goldberg_Gli_vs_Control_by_age_gender_with_percent.csv",
  row.names = FALSE
)

write.csv(
  goldberg_eye_results,
  "Goldberg_Right_minus_Left_by_age_gender.csv",
  row.names = FALSE
)

write.csv(
  goldberg_type3_results,
  "Goldberg_mixed_model_Type3_tests.csv",
  row.names = FALSE
)

saveRDS(
  goldberg_model,
  "Goldberg_mixed_model.rds"
)

# Save the Gutmann 17-week results

write.csv(
  gutmann17_counts,
  "Gutmann_17week_animal_counts.csv",
  row.names = FALSE
)

write.csv(
  gutmann17_genotype_means,
  "Gutmann_17week_adjusted_genotype_means.csv",
  row.names = FALSE
)

write.csv(
  gutmann17_genotype_results,
  "Gutmann_17week_Gli_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  gutmann17_eye_results,
  "Gutmann_17week_Right_minus_Left.csv",
  row.names = FALSE
)

write.csv(
  gutmann17_type3_results,
  "Gutmann_17week_mixed_model_Type3_tests.csv",
  row.names = FALSE
)

write.csv(
  gutmann17_percent_result,
  "gutmann17_percent_result.csv",
  row.names = FALSE
)

saveRDS(
  gutmann17_model,
  "Gutmann_17week_mixed_model.rds"
)

# Creating a seperate Goldberg dataset to look at the difference in 
#  all Goldberg riased mice since no group was significant

library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)
library(splines)

ctrl_label <- "NF1 fl/fl, Cre:- Ctrl"
gli_label  <- "NF1 KO/fl, Cre:+, Gli"

# Create an all Goldberg dataset

goldberg_overall_dat <- dat %>%
  mutate(
    ID = trimws(as.character(ID)),
    Raised = trimws(as.character(Raised)),
    Gender = tolower(trimws(as.character(Gender))),
    Genotype_original = trimws(as.character(Genotype)),
    
    Density_value = parse_number(
      as.character(`Density`)
    ),
    
    Age_weeks = parse_number(
      as.character(`Age at Death (weeks)`)
    ),
    
    Eye = case_when(
      toupper(trimws(as.character(Eye))) %in% c("L", "LEFT") ~ "Left",
      toupper(trimws(as.character(Eye))) %in% c("R", "RIGHT") ~ "Right",
      TRUE ~ NA_character_
    ),
    
    Genotype_group = case_when(
      Genotype_original == ctrl_label ~ "Control",
      Genotype_original == gli_label  ~ "Gli",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    Raised == "Goldberg",
    Genotype_group %in% c("Control", "Gli"),
    !is.na(ID),
    ID != "",
    !is.na(Density_value),
    !is.na(Age_weeks),
    !is.na(Gender),
    !is.na(Eye)
  ) %>%
  mutate(
    Genotype_group = factor(
      Genotype_group,
      levels = c("Control", "Gli")
    ),
    Gender = factor(Gender),
    Eye = factor(
      Eye,
      levels = c("Left", "Right")
    ),
    Animal_ID = factor(ID)
  ) %>%
  group_by(
    Animal_ID,
    Genotype_group,
    Gender,
    Age_weeks,
    Eye
  ) %>%
  summarise(
    Density_value = mean(Density_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  droplevels()

# Verify the ages are present

goldberg_overall_dat %>%
  distinct(Animal_ID, Age_weeks, Gender, Genotype_group) %>%
  arrange(Age_weeks) %>%
  print(row.names = FALSE)

# specificlly look at mice that are 8 and 20 weeks of age

goldberg_overall_dat %>%
  distinct(Animal_ID, Age_weeks, Gender, Genotype_group) %>%
  filter(
    between(Age_weeks, 7.5, 8.5) |
      between(Age_weeks, 19.5, 20.5)
  ) %>%
  arrange(Age_weeks) %>%
  print(row.names = FALSE)

# Fit the all-age overall mixed model

goldberg_overall_model <- lmer(
  Density_value ~
    Genotype_group +
    Age_weeks +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = TRUE
)

summary(goldberg_overall_model)

# Obtain the overall Glioma versus Control comparison

emm_options(lmer.df = "satterthwaite")

goldberg_overall_emmeans <- emmeans(
  goldberg_overall_model,
  ~ Genotype_group,
  lmer.df = "satterthwaite"
)

goldberg_overall_means <- as.data.frame(
  summary(
    goldberg_overall_emmeans,
    infer = c(TRUE, TRUE)
  )
)

goldberg_overall_result <- as.data.frame(
  summary(
    pairs(
      goldberg_overall_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

goldberg_overall_means
goldberg_overall_result


# Check the overlap between ages

age_overlap_check <- goldberg_overall_dat %>%
  distinct(Animal_ID, Genotype_group, Age_weeks) %>%
  group_by(Genotype_group) %>%
  summarise(
    N_animals = n_distinct(Animal_ID),
    Minimum_age = min(Age_weeks),
    Maximum_age = max(Age_weeks),
    Median_age = median(Age_weeks),
    .groups = "drop"
  )

age_overlap_check

# Save the all-age results

write.csv(
  goldberg_overall_means,
  "Goldberg_all_ages_adjusted_genotype_means.csv",
  row.names = FALSE
)

write.csv(
  goldberg_overall_result,
  "Goldberg_all_ages_Gli_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  age_overlap_check,
  "Goldberg_all_ages_genotype_age_overlap.csv",
  row.names = FALSE
)

saveRDS(
  goldberg_overall_model,
  "Goldberg_all_ages_mixed_model.rds"
)


# Fit the all-age additive model

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

emm_options(lmer.df = "satterthwaite")

goldberg_overall_dat <- goldberg_overall_dat %>%
  mutate(
    Age_centered = Age_weeks -
      mean(Age_weeks, na.rm = TRUE),
    
    Genotype_group = relevel(
      factor(Genotype_group),
      ref = "Control"
    ),
    
    Gender = relevel(
      factor(Gender),
      ref = "female"
    ),
    
    Eye = factor(
      Eye,
      levels = c("Left", "Right")
    )
  ) %>%
  droplevels()

goldberg_all_age_model <- lmer(
  Density_value ~
    Genotype_group +
    Age_centered +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = TRUE
)

summary(goldberg_all_age_model)

# Obtain the overall tests

all_age_type3 <- anova(
  goldberg_all_age_model,
  type = 3,
  ddf = "Satterthwaite"
)

all_age_type3

# Estimate the age effect

age_trend <- emtrends(
  goldberg_all_age_model,
  ~ 1,
  var = "Age_centered",
  lmer.df = "satterthwaite"
)

age_trend_result <- as.data.frame(
  summary(
    age_trend,
    infer = c(TRUE, TRUE)
  )
)

age_trend_result

# Express the change per 10 weeks

age_trend_per_10_weeks <- age_trend_result %>%
  mutate(
    Estimate_per_10_weeks = Age_centered.trend * 10,
    Lower_CI_per_10_weeks = lower.CL * 10,
    Upper_CI_per_10_weeks = upper.CL * 10
  )

age_trend_per_10_weeks

# Estimate the adjusted gender difference

gender_emmeans <- emmeans(
  goldberg_all_age_model,
  ~ Gender,
  lmer.df = "satterthwaite"
)

gender_contrast <- pairs(
  gender_emmeans,
  reverse = TRUE
)

gender_mean_results <- as.data.frame(
  summary(
    gender_emmeans,
    infer = c(TRUE, TRUE)
  )
)

gender_difference_result <- as.data.frame(
  summary(
    gender_contrast,
    infer = c(TRUE, TRUE)
  )
)

gender_mean_results
gender_difference_result

# Estimate the adjusted eye difference

eye_emmeans <- emmeans(
  goldberg_all_age_model,
  ~ Eye,
  lmer.df = "satterthwaite"
)

eye_contrast <- pairs(
  eye_emmeans,
  reverse = TRUE
)

eye_mean_results <- as.data.frame(
  summary(
    eye_emmeans,
    infer = c(TRUE, TRUE)
  )
)

eye_difference_result <- as.data.frame(
  summary(
    eye_contrast,
    infer = c(TRUE, TRUE)
  )
)

eye_mean_results
eye_difference_result

# Reproduce the adjusted genotype difference

genotype_emmeans <- emmeans(
  goldberg_all_age_model,
  ~ Genotype_group,
  lmer.df = "satterthwaite"
)

genotype_difference_result <- as.data.frame(
  summary(
    pairs(
      genotype_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

genotype_difference_result

# Test whether genotype effects depend on age, gender, or eye

additive_model_ML <- lmer(
  Density_value ~
    Genotype_group +
    Age_centered +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = FALSE
)

interaction_model_ML <- lmer(
  Density_value ~
    Genotype_group * Age_centered +
    Genotype_group * Gender +
    Genotype_group * Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = FALSE
)

anova(
  additive_model_ML,
  interaction_model_ML
)

# Genotype by age interaction

genotype_age_model_ML <- update(
  additive_model_ML,
  . ~ . + Genotype_group:Age_centered
)

genotype_by_age_test <- anova(
  additive_model_ML,
  genotype_age_model_ML
)

genotype_by_age_test

# Gentype by gender interaction

genotype_gender_model_ML <- update(
  additive_model_ML,
  . ~ . + Genotype_group:Gender
)

genotype_by_gender_test <- anova(
  additive_model_ML,
  genotype_gender_model_ML
)

genotype_by_gender_test

# Genotype by eye interaction

genotype_eye_model_ML <- update(
  additive_model_ML,
  . ~ . + Genotype_group:Eye
)

genotype_by_eye_test <- anova(
  additive_model_ML,
  genotype_eye_model_ML
)

genotype_by_eye_test

# Check whether the age effect is nonlinear

library(splines)

linear_age_model_ML <- lmer(
  Density_value ~
    Genotype_group +
    Age_centered +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = FALSE
)

nonlinear_age_model_ML <- lmer(
  Density_value ~
    Genotype_group +
    ns(Age_weeks, df = 3) +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = FALSE
)

nonlinear_age_test <- anova(
  linear_age_model_ML,
  nonlinear_age_model_ML
)

nonlinear_age_test

# Save the results

type3_results <- tibble::rownames_to_column(
  as.data.frame(all_age_type3),
  var = "Model_term"
)

write.csv(
  type3_results,
  "Goldberg_all_ages_Type3_tests.csv",
  row.names = FALSE
)

write.csv(
  genotype_difference_result,
  "Goldberg_all_ages_Gli_vs_Control.csv",
  row.names = FALSE
)

write.csv(
  age_trend_result,
  "Goldberg_all_ages_age_trend.csv",
  row.names = FALSE
)

write.csv(
  gender_mean_results,
  "Goldberg_all_ages_gender_adjusted_means.csv",
  row.names = FALSE
)

write.csv(
  gender_difference_result,
  "Goldberg_all_ages_male_minus_female.csv",
  row.names = FALSE
)

write.csv(
  eye_mean_results,
  "Goldberg_all_ages_eye_adjusted_means.csv",
  row.names = FALSE
)

write.csv(
  eye_difference_result,
  "Goldberg_all_ages_right_minus_left.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(genotype_by_age_test),
  "Goldberg_all_ages_genotype_by_age_test.csv"
)

write.csv(
  as.data.frame(genotype_by_gender_test),
  "Goldberg_all_ages_genotype_by_gender_test.csv"
)

write.csv(
  as.data.frame(genotype_by_eye_test),
  "Goldberg_all_ages_genotype_by_eye_test.csv"
)

write.csv(
  as.data.frame(nonlinear_age_test),
  "Goldberg_all_ages_nonlinear_age_test.csv"
)

saveRDS(
  goldberg_all_age_model,
  "Goldberg_all_ages_additive_mixed_model.rds"
)

# Getting the slope for age change over time

library(emmeans)

age_trend <- emtrends(
  goldberg_all_age_model,
  ~ 1,
  var = "Age_centered",
  lmer.df = "satterthwaite"
)

age_trend_result <- as.data.frame(
  summary(
    age_trend,
    infer = c(TRUE, TRUE)
  )
)

age_trend_result

# Save the age result

write.csv(
  as.data.frame(age_trend_result),
  "age_trend_result.csv"
)


# Display the slope for each genotype

library(lme4)
library(lmerTest)
library(emmeans)

genotype_age_model_REML <- lmer(
  Density_value ~
    Genotype_group * Age_centered +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = TRUE
)

# Estimate the slope for each genotype

age_slopes_by_genotype <- emtrends(
  genotype_age_model_REML,
  ~ Genotype_group,
  var = "Age_centered",
  lmer.df = "satterthwaite"
)

age_slope_results <- as.data.frame(
  summary(
    age_slopes_by_genotype,
    infer = c(TRUE, TRUE)
  )
)

age_slope_results

# Directly compare the slopes

age_slope_difference <- as.data.frame(
  summary(
    pairs(
      age_slopes_by_genotype,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

age_slope_difference

# Create the date to plot density over time for the male/female

library(dplyr)
library(ggplot2)
library(scales)

age_plot_dat <- goldberg_overall_dat %>%
  group_by(
    Animal_ID,
    Age_weeks,
    Gender,
    Genotype_group
  ) %>%
  summarise(
    Mean_density = mean(Density_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Plot_group = interaction(
      Gender,
      Genotype_group,
      sep = " ",
      drop = TRUE
    ),
    
    Plot_group = factor(
      Plot_group,
      levels = c(
        "female Control",
        "female Gli",
        "male Control",
        "male Gli"
      ),
      labels = c(
        "Female Control",
        "Female Glioma",
        "Male Control",
        "Male Glioma"
      )
    )
  ) %>%
  filter(
    !is.na(Age_weeks),
    !is.na(Mean_density),
    !is.na(Plot_group)
  )

# Set axis

max_x <- ceiling(
  max(age_plot_dat$Age_weeks, na.rm = TRUE) / 10
) * 10

major_x <- seq(0, max_x, by = 10)
minor_x <- setdiff(
  seq(0, max_x, by = 2),
  major_x
)

max_y <- 4000

major_y <- seq(0, max_y, by = 500)
minor_y <- setdiff(
  seq(0, max_y, by = 100),
  major_y
)

# Plot the 4 groups age trends

age_trend_four_group_plot <- ggplot(
  age_plot_dat,
  aes(
    x = Age_weeks,
    y = Mean_density,
    color = Plot_group,
    fill = Plot_group,
    shape = Plot_group,
    linetype = Plot_group
  )
) +
  geom_point(
    size = 2.6,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1,
    alpha = 0.10
  ) +
  scale_color_manual(
    values = c(
      "Female Control" = "#0072B2",
      "Female Glioma"  = "#D55E00",
      "Male Control"   = "#009E73",
      "Male Glioma"    = "#CC79A7"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Female Control" = "#0072B2",
      "Female Glioma"  = "#D55E00",
      "Male Control"   = "#009E73",
      "Male Glioma"    = "#CC79A7"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Female Control" = 16,
      "Female Glioma"  = 17,
      "Male Control"   = 15,
      "Male Glioma"    = 18
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Female Control" = "solid",
      "Female Glioma"  = "solid",
      "Male Control"   = "solid",
      "Male Glioma"    = "solid"
    )
  ) +
  scale_x_continuous(
    limits = c(0, max_x),
    breaks = major_x,
    minor_breaks = minor_x,
    expand = expansion(mult = c(0, 0.02)),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_continuous(
    breaks = seq(0, 4000, by = 500),
    minor_breaks = seq(0, 4000, by = 100),
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0)),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  coord_cartesian(
    ylim = c(0, 4000)
  ) +
  labs(
    x = "Age at death (weeks)",
    y = expression("Mean cell density (cells/mm"^2*")"),
    color = "Group",
    fill = "Group",
    shape = "Group",
    linetype = "Group",
    title = "Cell density trends by age, gender, and genotype",
    subtitle = paste(
      "Each point represents the mean of the left and right eyes",
      "for one mouse"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

age_trend_four_group_plot


# Save as a vector image pdf

ggsave(
  filename = "Goldberg_age_trends_by_gender_and_genotype.pdf",
  plot = age_trend_four_group_plot,
  device = cairo_pdf,
  width = 9,
  height = 6.5,
  units = "in",
  bg = "white"
)


# Compare the Age x Sex interaction while adjusting for genotype and paired eyes

library(lme4)
library(lmerTest)

# Common age slope for males and females
age_gender_additive_ML <- lmer(
  Density_value ~
    Genotype_group +
    Age_centered +
    Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = FALSE
)

# Different age slopes allowed for males and females
age_gender_interaction_ML <- lmer(
  Density_value ~
    Genotype_group +
    Age_centered * Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = FALSE
)

age_by_gender_test <- anova(
  age_gender_additive_ML,
  age_gender_interaction_ML
)

age_by_gender_test


# Estimate slopes by sex

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

age_gender_interaction_REML <- lmer(
  Density_value ~
    Genotype_group +
    Age_centered * Gender +
    Eye +
    (1 | Animal_ID),
  data = goldberg_overall_dat,
  REML = TRUE
)

age_slopes_by_gender <- emtrends(
  age_gender_interaction_REML,
  ~ Gender,
  var = "Age_centered",
  lmer.df = "satterthwaite"
)

age_slopes_by_gender_results <- as.data.frame(
  summary(
    age_slopes_by_gender,
    infer = c(TRUE, TRUE)
  )
)

age_slopes_by_gender_results


# Directly compare the male and female slope

age_slope_gender_difference <- as.data.frame(
  summary(
    pairs(
      age_slopes_by_gender,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

age_slope_gender_difference

# Save the slope data

write.csv(
  age_slopes_by_gender_results,
  file = "Goldberg_age_slopes_by_gender.csv",
  row.names = FALSE
)

write.csv(
  age_slope_gender_difference,
  file = "Goldberg_male_vs_female_age_slope_difference.csv",
  row.names = FALSE
)

# Calculate the % reduction in each Goldberg raised group

library(dplyr)
library(tidyr)
library(emmeans)

group_emmeans <- emmeans(
  goldberg_model,
  ~ Genotype_group | Comparison_group,
  lmer.df = "satterthwaite"
)

group_mean_table <- as.data.frame(group_emmeans)

# Prepare the data for Gutmann 24 week and Charles river

library(dplyr)
library(readr)
library(stringr)

# Remove any incomplete version
if (exists("week24_dat")) rm(week24_dat)

gli_label <- "NF1 KO/fl, Cre:+, Gli"

# Confirm required columns exist
required_columns <- c(
  "ID",
  "Raised",
  "Genotype",
  "Eye",
  "Age at Death (weeks)",
  "Density"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Create the 24 week dataset

library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)

gli_label <- "NF1 KO/fl, Cre:+, Gli"

week24_dat <- dat %>%
  mutate(
    ID_original = trimws(as.character(ID)),
    Raised_original = trimws(as.character(Raised)),
    Genotype_original = trimws(as.character(Genotype)),
    
    Age_weeks = parse_number(
      as.character(`Age at Death (weeks)`)
    ),
    
    Density_value = as.numeric(Density),
    
    Eye_standard = case_when(
      toupper(trimws(as.character(Eye))) %in% c("L", "LEFT") ~ "Left",
      toupper(trimws(as.character(Eye))) %in% c("R", "RIGHT") ~ "Right",
      TRUE ~ NA_character_
    ),
    
    Cohort = case_when(
      Raised_original == "Charles River" &
        Genotype_original == "C57Blk/6" ~
        "CharlesRiver_Control",
      
      Raised_original == "Gutmann raised" &
        Genotype_original == gli_label ~
        "Gutmann_Nf1_OPG",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    between(Age_weeks, 23, 25),
    !is.na(Cohort),
    !is.na(ID_original),
    ID_original != "",
    !is.na(Density_value),
    !is.na(Eye_standard)
  ) %>%
  mutate(
    Cohort = factor(
      Cohort,
      levels = c(
        "CharlesRiver_Control",
        "Gutmann_Nf1_OPG"
      )
    ),
    Eye = factor(
      Eye_standard,
      levels = c("Left", "Right")
    ),
    Animal_ID = interaction(
      Raised_original,
      ID_original,
      sep = "_",
      drop = TRUE
    )
  ) %>%
  group_by(
    Animal_ID,
    Cohort,
    Age_weeks,
    Eye
  ) %>%
  summarise(
    Density_value = mean(Density_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  droplevels()

# Verify the resulting data

print(as.data.frame(week24_dat), row.names = FALSE)

week24_dat %>%
  distinct(Animal_ID, Cohort, Age_weeks) %>%
  count(Cohort, name = "N_mice")

# Fit the fixed model

emm_options(lmer.df = "satterthwaite")

week24_model <- lmer(
  Density_value ~ Cohort + Eye + (1 | Animal_ID),
  data = week24_dat,
  REML = TRUE
)

summary(week24_model)
isSingular(week24_model, tol = 1e-4)


# Obtain the adjusted means and the cohort contrast

week24_cohort_emmeans <- emmeans(
  week24_model,
  ~ Cohort,
  lmer.df = "satterthwaite"
)

week24_cohort_means <- as.data.frame(
  summary(
    week24_cohort_emmeans,
    infer = c(TRUE, TRUE)
  )
)

week24_cohort_result <- as.data.frame(
  summary(
    pairs(
      week24_cohort_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

week24_cohort_means
week24_cohort_result


# Save the 24 week results

write.csv(
  week24_cohort_means,
  "Gutmann_CharlesRiver_24week_adjusted_means.csv",
  row.names = FALSE
)

write.csv(
  week24_cohort_result,
  "GutmannGli_vs_CharlesRiverControl_24week.csv",
  row.names = FALSE
)

saveRDS(
  week24_model,
  "Gutmann_CharlesRiver_24week_mixed_model.rds"
)
