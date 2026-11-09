library(tidyverse)
library(lubridate)
library(readr)

# ------------------------------------------------------------------------------
# 1. Read Data & Extract Date from `seqName`
# ------------------------------------------------------------------------------
df_raw <- readr::read_delim(
  "data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv",
  delim = ";"
)

df_clean <- df_raw %>%
  # Extract YYYY-MM-DD pattern directly from seqName
  mutate(
    extracted_date_str = str_extract(seqName, "\\d{4}-\\d{2}-\\d{2}"),
    date_clean = ymd(extracted_date_str)
  ) %>%
  # Remove records with unparseable or missing dates
  filter(!is.na(date_clean)) %>%
  # Define time periods
  mutate(
    Time_Period = case_when(
      date_clean >= ymd("2022-07-01") & date_clean <= ymd("2023-03-31") ~ "Jul 2022 - Mar 2023",
      date_clean >= ymd("2023-03-01") & date_clean <= ymd("2023-10-31") ~ "Mar 2023 - Oct 2023",
      date_clean >= ymd("2023-10-01") & date_clean <= ymd("2024-05-31") ~ "Oct 2023 - May 2024",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Time_Period)) %>%
  # Set correct chronological order for facets
  mutate(Time_Period = factor(Time_Period, levels = c(
    "Jul 2022 - Mar 2023",
    "Mar 2023 - Oct 2023",
    "Oct 2023 - May 2024"
  )))

# ------------------------------------------------------------------------------
# 2. Total Sequences per Period (for exact % prevalence)
# ------------------------------------------------------------------------------
period_totals <- df_clean %>%
  group_by(Time_Period) %>%
  summarise(total_seqs = n_distinct(seqName), .groups = "drop")

# ------------------------------------------------------------------------------
# 3. Process Alterations (Substitutions + Deletions + Insertions)
# ------------------------------------------------------------------------------
# Ensure aaInsertions column exists even if missing in raw data
if (!"aaInsertions" %in% colnames(df_clean)) {
  df_clean$aaInsertions <- NA_character_
}

gene_mutations <- df_clean %>%
  select(seqName, Time_Period, aaSubstitutions, aaDeletions, aaInsertions) %>%
  pivot_longer(
    cols = c(aaSubstitutions, aaDeletions, aaInsertions),
    names_to = "Type",
    values_to = "Mutation_List"
  ) %>%
  filter(!is.na(Mutation_List) & Mutation_List != "" & Mutation_List != "None") %>%
  separate_rows(Mutation_List, sep = ",") %>%
  mutate(
    raw_mut = trimws(Mutation_List),
    Gene = sub(":.*", "", raw_mut),
    # Add (del) or (ins) tags where applicable for clarity
    Mutation = case_when(
      Type == "aaDeletions"  & !str_detect(raw_mut, "(del|Deletion)") ~ paste0(raw_mut, " (del)"),
      Type == "aaInsertions" & !str_detect(raw_mut, "(ins|Insertion)") ~ paste0(raw_mut, " (ins)"),
      TRUE ~ raw_mut
    )
  ) %>%
  filter(Gene %in% c("S", "N"))

# Tally prevalence per time window
mutation_prevalence <- gene_mutations %>%
  count(Time_Period, Gene, Mutation, sort = TRUE) %>%
  left_join(period_totals, by = "Time_Period") %>%
  mutate(Prevalence = (n / total_seqs) * 100) %>%
  filter(Prevalence >= 5) # Threshold cutoff (>= 5%)

# ------------------------------------------------------------------------------
# 4. Plot Landscape
# ------------------------------------------------------------------------------
ggplot(mutation_prevalence, aes(x = reorder(Mutation, n), y = Prevalence, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_grid(Gene ~ Time_Period, scales = "free", space = "free_y") +
  scale_fill_manual(values = c("S" = "#2B5C8F", "N" = "#D95F02")) +
  labs(
    title = "S and N Gene Alterations (Substitutions, Deletions & Insertions ≥5%)",
    subtitle = "Extracted from GISAID seqName dates (July 2022 – May 2024)",
    x = "Substitution / Deletion / Insertion Locus",
    y = "Prevalence (% of sequences)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8, face = "bold"),
    panel.spacing = unit(0.8, "lines")
  )

print(period_totals)
