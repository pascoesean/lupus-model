# all of this data comes from this [https://doi.org/10.1002/art.38790] study; and can be downloaded
# from the ITN trial share network here [https://www.itntrialshare.org/project/Studies/ITN034AIPUBLIC/Study%20Data/begin.view]
# see README in the `raw_from_itn` directory for more info

library(tidyverse)


efficacy_one <- read_csv("data/raw_from_itn/ADEFF3_2024-05-22_10-50-24.csv")
efficacy_two <- read_csv("data/raw_from_itn/ADVISIT1_2024-05-22_10-51-40.csv")
lab_results <- read_csv("data/raw_from_itn/ADLB1_2024-05-22_10-52-48.csv")

platelet_counts <- lab_results |>
  filter(`Laboratory Test Name` == "Platelet Count")

neutro_lymph_ratio <- lab_results |>
  filter(`Laboratory Test Name` == "Neutrophils" | `Laboratory Test Name` == "Lymphocytes")

ef1_filtered <- efficacy_one |>
  select(`Participant ID`, Visit, `Treatment Group (Char)`, `Completed Study`, `Age (years)`, Race, Ethnicity, `Sex (Char)`, `Visit Number`, `Visit Study Day`, `Urinary Protein-to-Creatinine Ratio-spot`, `Renal Flare`, eGFR)

ef2_filtered <- efficacy_two |>
  select(`Participant ID`, Visit, `CH50 (units/mL)`, `IgA (mg/dL)`, `IgG (mg/dL)`)

ef2_filtered <- efficacy_two |>
  select(`Participant ID`, Visit, `C3 (mg/dL)`, `C4 (mg/dL)`, `CH50 (units/mL)`)

# need to keep doing...

merged <- ef1_filtered |>
  left_join(ef2_filtered, by = c("Participant ID", "Visit"))

above_fifteen <- merged |>
  group_by(`Participant ID`)|>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  filter(n > 15) |>
  select(`Participant ID`)
  #drop_na()

four_or_more_visits <- merged |>
  group_by(`Participant ID`)|>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  filter(n > 3) |>
  select(`Participant ID`)

merged_fifteen <- merged |>
  filter(`Participant ID` %in% above_fifteen$`Participant ID`)


merged_fifteen |>
  filter(Race == "Black" | Race == "White") |>
  ggplot(aes(x = `Visit Study Day`, y = `Urinary Protein-to-Creatinine Ratio-spot`)) +
  geom_line(aes(color = `Participant ID`)) +
  facet_wrap(~ Race, nrow = 1) +
  ggpubr::theme_pubr()

merged_fifteen |>
  filter(Race == "Black" | Race == "White") |>
  ggplot(aes(x = `Visit Study Day`, y = `C4 (mg/dL)`)) +
  geom_line(aes(color = `Participant ID`)) +
  facet_wrap(~ Race, nrow = 1)


merged |>
  filter(`Participant ID` == "ACCESS_966565") |>
  #select(`Visit Study Day`, eGFR) |>
  View()

merged_upc <- merged |>
  select(`Participant ID`, `Visit Study Day`, `Urinary Protein-to-Creatinine Ratio-spot`, `Age (years)`, Race, Ethnicity) |>
  drop_na()

above_fifteen <- merged_upc |>
  group_by(`Participant ID`)|>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  filter(n >= 14) |>
  select(`Participant ID`)

merged_upc |>
  filter(`Participant ID` %in% above_fifteen[[1]]) |>
  filter(Race == "Black" | Race == "White") |>
  ggplot(aes(x = `Visit Study Day`, y = `Urinary Protein-to-Creatinine Ratio-spot`)) +
  geom_line(aes(color = `Participant ID`)) +
  facet_wrap(~ Race, nrow = 1) +
  ggpubr::theme_pubr()

upc_for_export <- merged_upc |>
  filter(`Participant ID` %in% above_fifteen[[1]]) |>
  filter(Race == "Black" | Race == "White") |>
  dplyr::rename("ID" = "Participant ID", "day" = "Visit Study Day", "uPCR" = "Urinary Protein-to-Creatinine Ratio-spot")

readr::write_csv(upc_for_export, file = "data/upc_data.csv")
