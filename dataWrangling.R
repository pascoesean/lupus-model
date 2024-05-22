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
  select(`Participant ID`, Visit, `Treatment Group (Char)`, `Completed Study`, `Age (years)`, Race, Ethnicity, `Sex (Char)`, `Visit Number`, `Visit Study Day`, `Urinary Protein-to-Creatinine Ratio-spot`, `Blood Albumin (g/dL)`, `Renal Flare`, eGFR)

ef2_filtered <- efficacy_two |>
  select(`Participant ID`, Visit, `C3 (mg/dL)`, `C4 (mg/dL)`, `CH50 (units/mL)`, `ANA Result`, `ANA Titer, if Positive`, `IgA (mg/dL)`, `IgG (mg/dL)`, `IgM (mg/dL)`)

# need to keep doing...