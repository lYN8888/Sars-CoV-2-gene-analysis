library(tidyverse)
library(patchwork)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. Load Lineage Metadata (With Collection Dates)
# ------------------------------------------------------------------------------
lineage_data <- read_csv("data/lineage_data(samples 2022-2024).csv") %>%
  select(pid, collection_date, LINEAGE, lineage_group) %>%
  mutate(
    # Convert '7/27/2022' string to proper R Date object
    collection_date = mdy(collection_date) 
  )

# ------------------------------------------------------------------------------
# 2. Load D35 & Nextclade Data + Join Dates
# ------------------------------------------------------------------------------
d35 <- read_csv("data/household samples/d35_threshold_hhtransmission.csv") %>%
  mutate(
    sampleid     = str_replace(sampleid, "^ht", "HT"),
    household_id = str_sub(sampleid, 3, 5)
  )

nextclade <- read_delim(
  "data/household samples/nextclade_all_hhconsensus seq.csv", 
  delim = ";", 
  show_col_types = FALSE
) %>%
  mutate(
    sample_id    = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1      = str_detect(seqName, "Day1"),
    is_infected  = str_detect(seqName, "Infected")
  ) %>%
  # --- LEFT JOIN METADATA BY ID ---
  left_join(lineage_data, by = c("sample_id" = "pid")) %>%
  mutate(
    # Categorize into the 3 Time Periods based on matched date
    period = case_when(
      collection_date >= as.Date("2022-07-01") & collection_date < as.Date("2023-03-01") ~ "Jul 2022 - Mar 2023",
      collection_date >= as.Date("2023-03-01") & collection_date < as.Date("2023-10-01") ~ "Mar 2023 - Oct 2023",
      collection_date >= as.Date("2023-10-01") & collection_date <= as.Date("2024-05-31") ~ "Oct 2023 - May 2024",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024"))
  )

# ------------------------------------------------------------------------------
# 3. Define Household ID Categories & Classify Day 1 Cases
# ------------------------------------------------------------------------------
d35_hh_ids <- unique(d35$household_id)

community_infected_hh_ids <- nextclade %>%
  filter(is_infected & !(sample_id %in% d35$sampleid)) %>%
  pull(household_id) %>%
  unique() %>%
  setdiff(d35_hh_ids)

day1_classified <- nextclade %>%
  filter(is_day1, !is.na(period)) %>% # Keep only Day 1 cases with valid date period
  mutate(
    day1_group = case_when(
      household_id %in% d35_hh_ids ~ "Transmitted",
      household_id %in% community_infected_hh_ids ~ "Community Infected",
      TRUE ~ "Non-Transmitted"
    ),
    day1_group = factor(day1_group, levels = c("Non-Transmitted", "Transmitted", "Community Infected"))
  )

library(tidyverse)
library(stats)

# ==============================================================================
# 1. Prepare Groups & Sample Denominators
# ==============================================================================
day1_classified_periods <- day1_classified %>%
  mutate(
    # Explicitly map Transmitted vs Non-Transmitted (including Community)
    group_2class = case_when(
      day1_group == "Transmitted" ~ "Transmitted",
      TRUE ~ "Non-Transmitted"
    )
  ) %>%
  filter(!is.na(period) & !is.na(day1_group))

overall_totals <- day1_classified_periods %>%
  group_by(group_2class) %>%
  summarise(total_n = n_distinct(index), .groups = "drop")

total_trans_n <- overall_totals %>% filter(group_2class == "Transmitted") %>% pull(total_n)
total_non_trans_n <- overall_totals %>% filter(group_2class == "Non-Transmitted") %>% pull(total_n)

# ==============================================================================
# 2. Extract S & N Alterations (Substitutions, Deletions, Insertions)
# ==============================================================================
all_sn_alterations <- day1_classified_periods %>%
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
    mutation_class = case_when(
      alteration_type == "aaSubstitutions" ~ "Substitution",
      alteration_type == "aaDeletions" ~ "Deletion",
      alteration_type == "aaInsertions" ~ "Insertion"
    ),
    mutation = case_when(
      alteration_type == "aaDeletions" & !str_detect(raw_mut, "del$") ~ paste0(raw_mut, " (del)"),
      alteration_type == "aaInsertions" & !str_detect(raw_mut, "ins$") ~ paste0(raw_mut, " (ins)"),
      TRUE ~ raw_mut
    )
  ) %>%
  filter(gene %in% c("S", "N"))

# ==============================================================================
# 3. Perform Fisher's Exact Test + Benjamini-Hochberg FDR Adjustment
# ==============================================================================
fisher_all_mutations <- all_sn_alterations %>%
  group_by(gene, mutation_class, mutation, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  pivot_wider(names_from = group_2class, values_from = count, values_fill = 0) %>%
  mutate(
    Transmitted = if ("Transmitted" %in% colnames(.)) Transmitted else 0,
    `Non-Transmitted` = if ("Non-Transmitted" %in% colnames(.)) `Non-Transmitted` else 0
  ) %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(
      Transmitted, `Non-Transmitted`,
      total_trans_n - Transmitted, total_non_trans_n - `Non-Transmitted`
    ), nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  # FDR correction for multiple testing
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# PUBLICATION STANDARD: Filter by adjusted p-value (FDR q < 0.05)
sig_mutations_list <- fisher_all_mutations %>%
  filter(adj_p_value < 0.05) %>%
  pull(mutation)

# Fallback: If no mutation survives FDR < 0.05, filter by raw p < 0.05 with warning
if (length(sig_mutations_list) == 0) {
  warning("No mutations reached adj_p_value < 0.05. Falling back to raw p_value < 0.05 for display.")
  sig_mutations_list <- fisher_all_mutations %>%
    filter(p_value < 0.05) %>%
    pull(mutation)
}

# ==============================================================================
# 4. Compute Prevalence Metrics and Format Statistical Labels
# ==============================================================================
period_group_totals <- day1_classified_periods %>%
  group_by(period, group_2class) %>%
  summarise(total_period_n = n_distinct(index), .groups = "drop")

significant_timeline_data <- all_sn_alterations %>%
  filter(mutation %in% sig_mutations_list) %>%
  group_by(gene, mutation, period, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  complete(nesting(gene, mutation), period, group_2class, fill = list(count = 0)) %>%
  left_join(period_group_totals, by = c("period", "group_2class")) %>%
  mutate(prevalence_pct = round((count / total_period_n) * 100, 1))

# Format statistical strings according to journal standards
plot_data_pub <- significant_timeline_data %>%
  left_join(
    fisher_all_mutations %>% select(mutation, mutation_class, p_value, adj_p_value) %>% distinct(),
    by = "mutation"
  ) %>%
  mutate(
    p_str = if_else(p_value < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_value))),
    q_str = if_else(adj_p_value < 0.001, "q < 0.001", paste0("q = ", sprintf("%.3f", adj_p_value))),
    
    # Combined labels showing both raw p and FDR q-values
    facet_title = paste0(mutation, " (", mutation_class, ")\n", p_str, ", ", q_str),
    y_axis_label = paste0(mutation, " (", mutation_class, " | ", p_str, ", ", q_str, ")")
  )

# Color-blind friendly color scheme (Publication standard)
palette_pub <- c("Non-Transmitted" = "#4575B4", "Transmitted" = "#D73027")

# ==============================================================================
# 5. Publication Plot Option A: Facet Subplots
# ==============================================================================
fig_option_a <- ggplot(plot_data_pub, aes(x = period, y = prevalence_pct, fill = group_2class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = paste0(round(prevalence_pct), "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.5, size = 2.6, fontface = "bold"
  ) +
  scale_fill_manual(values = palette_pub) +
  scale_y_continuous(limits = c(0, 120), expand = c(0, 0)) +
  facet_wrap(~ facet_title, scales = "free_y", ncol = 3) +
  labs(
    title = "Prevalence of Significantly Enriched/Depleted S & N Alterations",
    subtitle = "Comparing Transmitted vs. Non-Transmitted Cases Across Time Periods (Fisher's Exact Test with BH FDR)",
    x = NULL,
    y = "Prevalence within Group (%)",
    fill = "Cohort Group:"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.4),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(face = "bold", size = 9),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_line(color = "#E0E0E0", linewidth = 0.3)
  )

print(fig_option_a)

# ==============================================================================
# 6. Publication Plot Option B: Facet Grid Matrix
# ==============================================================================
fig_option_b <- ggplot(plot_data_pub, aes(x = prevalence_pct, y = reorder(y_axis_label, prevalence_pct), fill = group_2class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = paste0(round(prevalence_pct), "%")),
    position = position_dodge(width = 0.8),
    hjust = -0.15, size = 2.6, fontface = "bold"
  ) +
  scale_fill_manual(values = palette_pub) +
  scale_x_continuous(limits = c(0, 128), expand = c(0, 0)) +
  facet_grid(period ~ gene, scales = "free_y", space = "free_y") +
  labs(
    title = "Prevalence of Significant S & N Alterations",
    subtitle = "Evaluation across epidemiologic windows (p = raw Fisher's p-value, q = Benjamini-Hochberg FDR)",
    x = "Prevalence within Group (%)",
    y = NULL,
    fill = "Household Indices:",
    caption = "Note: Due to sample size limitations, unadjusted raw p-values (p < 0.05) from Fisher's exact test are used\nfor exploratory analysis, with false discovery rates (q-values) reported for transparency."
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
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

print(fig_option_b)

