library(tidyverse)
library(lubridate)
library(readr)

# ==============================================================================
# 1. Load Household & Lineage Data
# ==============================================================================
lineage_data <-read_csv("Covid_study/data/lineage_data(samples 2022-2024).csv")%>%
  select(pid, collection_date, LINEAGE, lineage_group) %>%
  mutate(collection_date = mdy(collection_date))

d35 <- read_csv("Covid_study/data/household samples/d35_threshold_hhtransmission.csv") %>%
  
  mutate(
    sampleid     = str_replace(sampleid, "^ht", "HT"),
    household_id = str_sub(sampleid, 3, 5)
  )

nextclade <- read_delim("Covid_study/data/household samples/nextclade_all_hhconsensus seq.csv", 
                        delim = ";", 
                        show_col_types = FALSE
) %>%
  mutate(
    sample_id    = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1      = str_detect(seqName, "Day1"),
    is_infected  = str_detect(seqName, "Infected")
  ) %>%
  left_join(lineage_data, by = c("sample_id" = "pid")) %>%
  mutate(
    period = case_when(
      collection_date >= as.Date("2022-07-01") & collection_date < as.Date("2023-03-01") ~ "Jul 2022 - Mar 2023",
      collection_date >= as.Date("2023-03-01") & collection_date < as.Date("2023-10-01") ~ "Mar 2023 - Oct 2023",
      collection_date >= as.Date("2023-10-01") & collection_date <= as.Date("2024-05-31") ~ "Oct 2023 - May 2024",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024"))
  )

# Define household transmission groups
d35_hh_ids <- unique(d35$household_id)

community_infected_hh_ids <- nextclade %>%
  filter(is_infected & !(sample_id %in% d35$sampleid)) %>%
  pull(household_id) %>%
  unique() %>%
  setdiff(d35_hh_ids)

day1_classified_periods <- nextclade %>%
  filter(is_day1, !is.na(period)) %>%
  mutate(
    day1_group = case_when(
      household_id %in% d35_hh_ids ~ "Transmitted",
      household_id %in% community_infected_hh_ids ~ "Community Infected",
      TRUE ~ "Non-Transmitted"
    ),
    group_2class = if_else(day1_group == "Transmitted", "Transmitted", "Non-Transmitted")
  )

# ==============================================================================
# 2. Extract Significant Household Candidate Loci
# ==============================================================================
overall_totals <- day1_classified_periods %>%
  group_by(group_2class) %>%
  summarise(total_n = n_distinct(index), .groups = "drop")

total_trans_n <- overall_totals %>% filter(group_2class == "Transmitted") %>% pull(total_n)
total_non_trans_n <- overall_totals %>% filter(group_2class == "Non-Transmitted") %>% pull(total_n)

hh_sn_alterations <- day1_classified_periods %>%
  select(index, seqName, group_2class, period, aaSubstitutions, aaDeletions, aaInsertions) %>%
  pivot_longer(
    cols = c(aaSubstitutions, aaDeletions, aaInsertions),
    names_to = "alteration_type",
    values_to = "mutation_list"
  ) %>%
  filter(!is.na(mutation_list) & mutation_list != "") %>%
  separate_rows(mutation_list, sep = ",") %>%
  mutate(
    raw_mut = trimws(mutation_list),
    gene = str_extract(raw_mut, "^[^:]+"),
    mutation = case_when(
      alteration_type == "aaDeletions" & !str_detect(raw_mut, "(del|Deletion)") ~ paste0(raw_mut, " (del)"),
      alteration_type == "aaInsertions" & !str_detect(raw_mut, "(ins|Insertion)") ~ paste0(raw_mut, " (ins)"),
      TRUE ~ raw_mut
    )
  ) %>%
  filter(gene %in% c("S", "N"))

fisher_all_mutations <- hh_sn_alterations %>%
  group_by(gene, mutation, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  pivot_wider(names_from = group_2class, values_from = count, values_fill = 0)

# Ensure both required group columns exist
if (!"Transmitted" %in% names(fisher_all_mutations)) fisher_all_mutations$Transmitted <- 0
if (!"Non-Transmitted" %in% names(fisher_all_mutations)) fisher_all_mutations$`Non-Transmitted` <- 0

# Safe Fisher's Exact Test function
safe_fisher_full <- function(a, b, N_trans, N_nontrans) {
  c_val <- max(0, N_trans - a)
  d_val <- max(0, N_nontrans - b)
  
  mat <- matrix(c(a, c_val, b, d_val), nrow = 2, byrow = FALSE)
  
  if (dim(mat)[1] == 2 && dim(mat)[2] == 2 && all(mat >= 0)) {
    res <- fisher.test(mat)
    return(tibble(
      p_value = res$p.value,
      odds_ratio = as.numeric(res$estimate),
      or_lower = res$conf.int[1],
      or_upper = res$conf.int[2]
    ))
  } else {
    return(tibble(
      p_value = NA_real_,
      odds_ratio = NA_real_,
      or_lower = NA_real_,
      or_upper = NA_real_
    ))
  }
}

# Run Fisher's test safely with OR and CI extraction
fisher_all_mutations <- fisher_all_mutations %>%
  rowwise() %>%
  mutate(
    stats = list(safe_fisher_full(Transmitted, `Non-Transmitted`, total_trans_n, total_non_trans_n))
  ) %>%
  unnest(stats) %>%
  filter(!is.na(p_value)) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

target_candidate_mutations <- fisher_all_mutations %>%
  filter(p_value < 0.05) %>%
  pull(mutation)

# Format Household prevalence dataset for downstream binding
hh_period_group_totals <- day1_classified_periods %>%
  group_by(period, group_2class) %>%
  summarise(total_period_n = n_distinct(index), .groups = "drop")

hh_summary <- hh_sn_alterations %>%
  filter(mutation %in% target_candidate_mutations) %>%
  group_by(gene, mutation, period, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  complete(nesting(gene, mutation), period, group_2class, fill = list(count = 0)) %>%
  left_join(hh_period_group_totals, by = c("period", "group_2class")) %>%
  mutate(
    Prevalence = (count / total_period_n) * 100,
    Group = group_2class,
    Time_Period = period
  ) %>%
  select(Time_Period, gene, mutation, Group, Prevalence)

# ==============================================================================
# 3. Load GISAID Background Data & Compute Prevalence
# ==============================================================================
df_gisaid_raw <- readr::read_delim(
  "Covid_study/data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv",
  delim = ";",
  show_col_types = FALSE
)

if (!"aaInsertions" %in% colnames(df_gisaid_raw)) {
  df_gisaid_raw$aaInsertions <- NA_character_
}

df_gisaid_clean <- df_gisaid_raw %>%
  mutate(
    extracted_date_str = str_extract(seqName, "\\d{4}-\\d{2}-\\d{2}"),
    date_clean = ymd(extracted_date_str)
  ) %>%
  filter(!is.na(date_clean)) %>%
  mutate(
    Time_Period = case_when(
      date_clean >= ymd("2022-07-01") & date_clean <= ymd("2023-03-31") ~ "Jul 2022 - Mar 2023",
      date_clean >= ymd("2023-03-01") & date_clean <= ymd("2023-10-31") ~ "Mar 2023 - Oct 2023",
      date_clean >= ymd("2023-10-01") & date_clean <= ymd("2024-05-31") ~ "Oct 2023 - May 2024",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Time_Period))

gisaid_totals <- df_gisaid_clean %>%
  group_by(Time_Period) %>%
  summarise(total_seqs = n_distinct(seqName), .groups = "drop")

# ==============================================================================
# 4. Format Household Summary with Counts & Denominators
# ==============================================================================
hh_summary <- hh_sn_alterations %>%
  filter(mutation %in% target_candidate_mutations) %>%
  group_by(gene, mutation, period, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  complete(nesting(gene, mutation), period, group_2class, fill = list(count = 0)) %>%
  left_join(hh_period_group_totals, by = c("period", "group_2class")) %>%
  mutate(
    Prevalence = (count / total_period_n) * 100,
    Group = group_2class,
    Time_Period = period,
    total_n = total_period_n
  ) %>%
  select(Time_Period, gene, mutation, Group, count, total_n, Prevalence)

# ==============================================================================
# 5. Format GISAID Summary with Counts & Denominators
# ==============================================================================
gisaid_summary <- df_gisaid_clean %>%
  select(seqName, Time_Period, aaSubstitutions, aaDeletions, aaInsertions) %>%
  pivot_longer(
    cols = c(aaSubstitutions, aaDeletions, aaInsertions),
    names_to = "Type",
    values_to = "Mut_List"
  ) %>%
  filter(!is.na(Mut_List) & Mut_List != "" & Mut_List != "None") %>%
  separate_rows(Mut_List, sep = ",") %>%
  mutate(
    raw_mut = trimws(Mut_List),
    gene = sub(":.*", "", raw_mut),
    mutation = case_when(
      Type == "aaDeletions" & !str_detect(raw_mut, "(del|Deletion)") ~ paste0(raw_mut, " (del)"),
      Type == "aaInsertions" & !str_detect(raw_mut, "(ins|Insertion)") ~ paste0(raw_mut, " (ins)"),
      TRUE ~ raw_mut
    )
  ) %>%
  filter(mutation %in% target_candidate_mutations) %>%
  group_by(Time_Period, gene, mutation) %>%
  summarise(count = n_distinct(seqName), .groups = "drop") %>%
  complete(nesting(gene, mutation), Time_Period, fill = list(count = 0)) %>%
  left_join(gisaid_totals, by = "Time_Period") %>%
  mutate(
    Prevalence = (count / total_seqs) * 100,
    Group = "GISAID (Background)",
    Time_Period = factor(Time_Period, levels = levels(hh_summary$Time_Period)),
    total_n = total_seqs
  ) %>%
  select(Time_Period, gene, mutation, Group, count, total_n, Prevalence)

# ==============================================================================
# 6. Combine & Create Formatted Display Labels
# ==============================================================================
combined_df <- bind_rows(hh_summary, gisaid_summary) %>%
  mutate(
    Group = factor(Group, levels = c("GISAID (Background)", "Non-Transmitted", "Transmitted")),
    Time_Period = factor(Time_Period, levels = c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024")),
    bar_label = sprintf("%.1f%% (%d/%d)", Prevalence, count, total_n)
  ) %>%
  left_join(
    fisher_all_mutations %>% select(mutation, p_value, adj_p_value, odds_ratio, or_lower, or_upper),
    by = "mutation"
  ) %>%
  mutate(
    p_str  = if_else(p_value < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_value))),
    q_str  = if_else(adj_p_value < 0.001, "q < 0.001", paste0("q = ", sprintf("%.3f", adj_p_value))),
    or_str = sprintf("OR: %.2f (95%% CI: %.2f-%.2f)", odds_ratio, or_lower, or_upper),
    # Combine into facet label (or use or_str for table display)
    facet_y_label = paste0(mutation, "\n", or_str, "\n(", p_str, ", ", q_str, ")")
  )

# ==============================================================================
# 7. Plot with Detailed Sample Counts
# ==============================================================================
fill_colors <- c(
  "GISAID (Background)" = "#808080",
  "Non-Transmitted"      = "#2B5C8F",
  "Transmitted"          = "#D95F02"
)

fig_integrated <- ggplot(combined_df, aes(x = Prevalence, y = reorder(facet_y_label, Prevalence), fill = Group)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = bar_label),
    position = position_dodge(width = 0.85),
    hjust = -0.1,
    size = 2.4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = fill_colors) +
  # Expand upper x-axis limit to 135 to accommodate longer labels without clipping
  scale_x_continuous(limits = c(0, 135), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
  facet_grid(Time_Period ~ gene, scales = "free_y", space = "free_y") +
  labs(
    title = "Comparison of Sigificant S & N Alteration Prevalence within Household Transmission Index Groups with GISAID as a background in Northern Thailand",
    subtitle = "Format: Prevalence % (Positive Samples / Group Total N)",
    x = "Prevalence within Group (%)",
    y = NULL,
    fill = "Cohort Group:",
    caption = "Note: Due to sample size limitations, unadjusted raw p-values (p < 0.05) from Fisher's exact test are used\nfor exploratory analysis, with false discovery rates (q-values) reported for transparency."
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.4),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(face = "bold", size = 9),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 9),
    panel.spacing = unit(0.6, "lines")
  )

print(fig_integrated)
