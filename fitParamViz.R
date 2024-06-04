# looking at the different parameter values for different people

library(tidyverse)
source("pairwise_comparisons.R")

fpvals <- read_csv("data/fitted_parameter_values.csv") |>
  select(!`...1`)

upc_exported <- read_csv("data/upc_data.csv") |>
  select(ID, `Age (years)`, Race, Ethnicity, group) |>
  unique()

vals_w_demo <- fpvals |>
  left_join(y = upc_exported, by = "ID")


vals_w_demo |>
  ggplot(aes(x = Race, y = value, color = Race)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  #ggpubr::stat_compare_means() + 
  facet_grid(model_type~parameter) +
  labs(title = "Fitted Parameter Values by Race") +
  ggpubr::theme_pubr()

# looks same (in line with the treatment not really working)
vals_w_demo |>
  ggplot(aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  #ggpubr::stat_compare_means() + 
  facet_grid(model_type~parameter) +
  ggpubr::theme_pubr()

vals_w_demo |>
  ggplot(aes(x = Ethnicity, y = value, color = Ethnicity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  #ggpubr::stat_compare_means() + 
  facet_grid(model_type~parameter) +
  ggpubr::theme_pubr()

vals_w_demo |>
  mutate(r_e = str_c(Race, Ethnicity, sep = "_")) |>
  ggplot(aes(x = r_e, y = value, color = r_e)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  #ggpubr::stat_compare_means() + 
  facet_grid(model_type~parameter) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle = 90))

vals_w_demo |>
  ggplot(aes(x = `Age (years)`, y = value, color = parameter)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  ggpubr::stat_cor(method = "kendall", cor.coef.name = "tau", label.y = 0.75) +
  facet_grid(model_type~parameter) +
  labs(title = "Fitted Parameter Values by Age", subtitle = "Kendall's Tau used for correlations") +
  ggpubr::theme_pubr()


vals_w_demo |>
  ggplot(aes(x = `Age (years)`, y = value, color = parameter)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  ggpubr::stat_cor(method = "kendall", cor.coef.name = "tau", label.y = 0.75) +
  facet_grid(model_type*Race~parameter) +
  ggpubr::theme_pubr()


