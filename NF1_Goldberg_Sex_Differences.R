# Sex based differences in the Goldberg data

library(dplyr)
library(readr)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)

ctrl_label <- "NF1 fl/fl, Cre:- Ctrl"
gli_label  <- "NF1 KO/fl, Cre:+, Gli"

goldberg_density_dat <- dat_6DC %>%
  transmute(
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
    
    Density_value = parse_number(
      as.character(`Density (cells/mm^2)`)
    ),
    
    Sex = case_when(
      str_to_lower(trimws(as.character(Sex))) %in%
        c("female", "f") ~ "female",
      
      str_to_lower(trimws(as.character(Sex))) %in%
        c("male", "m") ~ "male",
      
      TRUE ~ NA_character_
    ),
    
    # Corrected: uppercase values compared with uppercase labels
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
  filter(
    str_to_lower(Raised_original) == "goldberg",
    Genotype_group %in% c("Control", "Gli"),
    !is.na(ID_original),
    ID_original != "",
    !is.na(Density_value),
    !is.na(Age_weeks),
    !is.na(Sex),
    !is.na(Eye)
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
    
    Animal_ID = factor(ID_original)
  ) %>%
  group_by(
    Animal_ID,
    Genotype_group,
    Sex,
    Age_weeks,
    Eye
  ) %>%
  summarise(
    Density_value = mean(
      Density_value,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Age_centered =
      Age_weeks - mean(Age_weeks, na.rm = TRUE)
  ) %>%
  droplevels()

# Test whether the genotype effect differs by sex

density_genotype_sex_additive_ML <- lmer(
  Density_value ~
    Genotype_group +
    Sex +
    Age_centered +
    Eye +
    (1 | Animal_ID),
  data = goldberg_density_dat,
  REML = FALSE
)

density_genotype_sex_interaction_ML <- lmer(
  Density_value ~
    Genotype_group * Sex +
    Age_centered +
    Eye +
    (1 | Animal_ID),
  data = goldberg_density_dat,
  REML = FALSE
)

density_genotype_by_sex_test <- anova(
  density_genotype_sex_additive_ML,
  density_genotype_sex_interaction_ML
)

density_genotype_by_sex_test


# Sex specific emmeans

density_genotype_sex_model <- lmer(
  Density_value ~
    Genotype_group * Sex +
    Age_centered +
    Eye +
    (1 | Animal_ID),
  data = goldberg_density_dat,
  REML = TRUE
)

density_genotype_sex_emmeans <- emmeans(
  density_genotype_sex_model,
  ~ Genotype_group | Sex,
  lmer.df = "satterthwaite"
)

density_genotype_sex_means <- as.data.frame(
  summary(
    density_genotype_sex_emmeans,
    infer = c(TRUE, TRUE)
  )
)

density_genotype_sex_results <- as.data.frame(
  summary(
    pairs(
      density_genotype_sex_emmeans,
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

density_genotype_sex_means
density_genotype_sex_results

# Directly compare the sex-specific genotype effects

density_genotype_sex_emmeans_all <- emmeans(
  density_genotype_sex_model,
  ~ Genotype_group * Sex,
  lmer.df = "satterthwaite"
)

density_genotype_sex_emmeans_all

density_genotype_effect_sex_difference <-
  contrast(
    density_genotype_sex_emmeans_all,
    method = list(
      "Male genotype effect - Female genotype effect" =
        c(1, -1, -1, 1)
    )
  )

density_genotype_effect_sex_difference_result <-
  as.data.frame(
    summary(
      density_genotype_effect_sex_difference,
      infer = c(TRUE, TRUE)
    )
  )

density_genotype_effect_sex_difference_result

# Primary overall genotype analysis

density_additive_model <- lmer(
  Density_value ~
    Genotype_group +
    Sex +
    Age_centered +
    Eye +
    (1 | Animal_ID),
  data = goldberg_density_dat,
  REML = TRUE
)

density_additive_type3 <- anova(
  density_additive_model,
  type = 3,
  ddf = "Satterthwaite"
)

density_additive_type3

# Estimate the overall genotype difference

density_overall_genotype_emmeans <- emmeans(
  density_additive_model,
  ~ Genotype_group,
  lmer.df = "satterthwaite"
)

density_overall_genotype_result <- as.data.frame(
  summary(
    pairs(
      density_overall_genotype_emmeans,
      reverse = TRUE
    ),
    infer = c(TRUE, TRUE)
  )
)

density_overall_genotype_result

# % decrease in female and male mice

# Obtain all four age adjusted means into one emmeans object

library(emmeans)
library(dplyr)

sex_emm <- emmeans(
  density_genotype_sex_model,
  ~ Genotype_group * Sex,
  lmer.df = "satterthwaite"
)

sex_emm_table <- as.data.frame(sex_emm)
sex_emm_table

# Calculate the percentage decrease and its standard error

# Covariance matrix of the adjusted means

V <- vcov(sex_emm)

percent_results <- lapply(
  levels(droplevels(sex_emm_table$Sex)),
  function(sex_name) {
    
    control_index <- which(
      sex_emm_table$Sex == sex_name &
        sex_emm_table$Genotype_group == "Control"
    )
    
    gli_index <- which(
      sex_emm_table$Sex == sex_name &
        sex_emm_table$Genotype_group == "Gli"
    )
    
    control_mean <- sex_emm_table$emmean[control_index]
    gli_mean <- sex_emm_table$emmean[gli_index]
    
    # Percentage decrease = 100 × (Control - Gli) / Control
    percent_decrease <- 100 *
      (control_mean - gli_mean) / control_mean
    
    # Delta-method derivatives
    derivative_control <- 100 * gli_mean / control_mean^2
    derivative_gli <- -100 / control_mean
    
    # Relevant 2 × 2 covariance matrix
    covariance_submatrix <- V[
      c(control_index, gli_index),
      c(control_index, gli_index),
      drop = FALSE
    ]
    
    gradient <- c(
      derivative_control,
      derivative_gli
    )
    
    percent_variance <- as.numeric(
      t(gradient) %*%
        covariance_submatrix %*%
        gradient
    )
    
    percent_SE <- sqrt(percent_variance)
    
    data.frame(
      Sex = sex_name,
      Control_adjusted_mean = control_mean,
      Gli_adjusted_mean = gli_mean,
      Percent_decrease_in_Gli = percent_decrease,
      Percent_SE = percent_SE,
      Lower_1_SE = percent_decrease - percent_SE,
      Upper_1_SE = percent_decrease + percent_SE,
      Lower_95_CI = percent_decrease - 1.96 * percent_SE,
      Upper_95_CI = percent_decrease + 1.96 * percent_SE
    )
  }
) %>%
  bind_rows()

percent_results

# Save data

write.csv(
  percent_results,
  "6DC_Goldberg_overall_percent_decrease_by_sex.csv",
  row.names = FALSE
)

density_genotype_by_sex_test_results <-
  tibble::rownames_to_column(
    as.data.frame(density_genotype_by_sex_test),
    var = "Model"
  )

write.csv(
  density_genotype_by_sex_test_results,
  "6DC_Goldberg_reported_density_genotype_by_sex_interaction_test.csv",
  row.names = FALSE
)

write.csv(
  density_genotype_sex_means,
  "6DC_Goldberg_reported_density_adjusted_means_by_sex.csv",
  row.names = FALSE
)

write.csv(
  density_genotype_sex_results,
  "6DC_Goldberg_reported_density_Gli_vs_Control_by_sex.csv",
  row.names = FALSE
)

write.csv(
  density_genotype_effect_sex_difference_result,
  "6DC_Goldberg_reported_density_genotype_effect_difference_by_sex.csv",
  row.names = FALSE
)

density_additive_type3_results <-
  tibble::rownames_to_column(
    as.data.frame(density_additive_type3),
    var = "Model_term"
  )

write.csv(
  density_additive_type3_results,
  "6DC_Goldberg_reported_density_additive_model_Type3.csv",
  row.names = FALSE
)

write.csv(
  density_overall_genotype_result,
  "6DC_Goldberg_reported_density_overall_Gli_vs_Control.csv",
  row.names = FALSE
)

density_overall_genotype_means <- as.data.frame(
  summary(
    density_overall_genotype_emmeans,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  density_overall_genotype_means,
  "6DC_Goldberg_reported_density_overall_adjusted_genotype_means.csv",
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
