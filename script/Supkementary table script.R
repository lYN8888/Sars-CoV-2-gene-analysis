library(tidyverse)
install.packages("openxlsx")
library(openxlsx)


library(tidyverse)
library(openxlsx)
library(stringr)
library(readr)
library(lubridate)

# Define Target Output Directory
target_dir <- "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data"
output_filepath <- file.path(target_dir, "Supplementary_Table_Unified_Lineages.xlsx")

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

# ==============================================================================
# 1. STRICT UNIFIED LINEAGE CLASSIFICATION FUNCTION
# ==============================================================================
assign_lineage_group_strict <- function(lineage) {
  case_when(
    # JN.1+
    grepl("^JN\\.|^JN1|^JN\\.1|^JN11|^JN14|^KP\\.|^KB\\.|^LB\\.|^JF\\.|^JS\\.", lineage) ~ "JN.1+",
    
    # XBB
    grepl("^XBB|^XBL\\.|^XBU\\.|^EG\\.|^GJ\\.|^GE\\.|^CH\\.|^BM\\.|^BQ\\.|^BL\\.|^BF\\.|^FL\\.|^FD\\.|^HK\\.|^LF\\.|^HH\\.|^FE\\.|^EU\\.|^GM\\.|^BR\\.|^BY\\.|^GR\\.|^HN\\.", lineage) ~ "XBB",
    
    # BA.4/BA.5
    grepl("^BA\\.4|^BA\\.5|^BE\\.|^BF\\.|^LZ|^FP\\.", lineage) ~ "BA.4/BA.5",
    
    # BA.2.75+ (Includes BA.2.9+, BA.2.10+, BA.1, BA.2, BN, BS, EJ)
    grepl("^B\\.1\\.1\\.529|^BA\\.1|^BA\\.2|^EJ\\.|^BN\\.|^BS\\.", lineage) ~ "BA.2.75+",
    
    # Delta
    grepl("^AY\\.|^GY\\.|^B\\.1\\.617|^B\\.1\\.617\\.2", lineage) ~ "Delta",
    
    # Alpha
    grepl("^B\\.1\\.1\\.7", lineage) ~ "Alpha",
    
    # Recombinants
    grepl("^X", lineage) ~ "Recombinants",
    
    TRUE ~ "Other"
  )
}

# ==============================================================================
# 2. HOUSEHOLD COHORT PROCESSING
# ==============================================================================
lineage_data_hh <- read_csv("~/Covid_study/data/lineage_data(samples 2022-2024).csv") %>%
  mutate(
    collection_date = as.Date(collection_date, format = "%m/%d/%Y"),
    Unified_Lineage_Group = assign_lineage_group_strict(lineage_group)
  )

df_hh_supp <- lineage_data_hh %>%
  filter(!is.na(collection_date), LINEAGE != "Unassigned") %>%
  select(
    Sample_ID = pid,
    Collection_Date = collection_date,
    Assigned_Pangolin_Lineage = LINEAGE,
    Unified_Lineage_Group
  ) %>%
  mutate(Cohort = "Chiang Mai Household Cohort")

# ==============================================================================
# 3. GISAID COHORT PROCESSING
# ==============================================================================
df_gisaid_raw <- readr::read_delim(
  "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv",
  delim = ";",
  show_col_types = FALSE
)

df_gisaid_supp <- df_gisaid_raw %>%
  separate(seqName, into = c("virus_name", "accession_id", "Collection_date"), sep = "\\|") %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    Assigned_Pangolin_Lineage = as.character(Nextclade_pango),
    Unified_Lineage_Group = assign_lineage_group_strict(Assigned_Pangolin_Lineage)
  ) %>%
  filter(!is.na(Collection_date)) %>%
  select(
    Sample_ID = accession_id,
    Collection_Date = Collection_date,
    Assigned_Pangolin_Lineage,
    Unified_Lineage_Group
  ) %>%
  mutate(Cohort = "GISAID Cohort") %>%
  distinct()

# ==============================================================================
# 4. EXPORT TO EXCEL
# ==============================================================================
master_supp_table <- bind_rows(df_hh_supp, df_gisaid_supp) %>%
  select(Cohort, Sample_ID, Collection_Date, Assigned_Pangolin_Lineage, Unified_Lineage_Group)

lineage_order <- c("BA.2.75+", "BA.4/BA.5", "XBB", "JN.1+", "Recombinants", "Alpha", "Delta", "Other")

summary_table <- master_supp_table %>%
  group_by(Cohort, Unified_Lineage_Group, Assigned_Pangolin_Lineage) %>%
  summarise(Total_Samples = n(), .groups = "drop") %>%
  mutate(Unified_Lineage_Group = factor(Unified_Lineage_Group, levels = lineage_order)) %>%
  arrange(Cohort, Unified_Lineage_Group, desc(Total_Samples))

wb <- createWorkbook()
addWorksheet(wb, "Lineage Mapping Summary")
addWorksheet(wb, "Full Supplementary Table")

writeData(wb, "Lineage Mapping Summary", summary_table)
writeData(wb, "Full Supplementary Table", master_supp_table)

saveWorkbook(wb, output_filepath, overwrite = TRUE)
cat("Strictly mapped supplementary table exported to:\n", output_filepath, "\n")

