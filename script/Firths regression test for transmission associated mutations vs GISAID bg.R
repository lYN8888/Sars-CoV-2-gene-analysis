library(tidyverse)
library(lubridate)
library(readr)
library(logistf) # Required for Firth's penalized logistic regression

# ==============================================================================
# 1. Load Household & Lineage Data
# ==============================================================================
lineage_data <- read_csv("Covid_study/data/lineage_data(samples 2022-2024).csv") %>%
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

# Define factor levels explicitly to preserve empty groups
all_periods <- c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024")
all_groups  <- c("GISAID (Background)", "Non-Transmitted", "Transmitted")

# ==============================================================================
# 2. Extract Alterations & Run Firth's Regression (Transmitted vs Non-Transmitted)
# ==============================================================================
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

# FIX 1: Complete period and group grid to ensure baseline totals exist for all periods/groups
period_group_totals <- day1_classified_periods %>%
  group_by(period, group_2class) %>%
  summarise(total_n = n_distinct(index), .groups = "drop") %>%
  complete(period = factor(all_periods, levels = all_periods), group_2class = c("Transmitted", "Non-Transmitted"), fill = list(total_n = 0))

alteration_indexed <- hh_sn_alterations %>%
  select(index, period, mutation) %>%
  distinct() %>%
  mutate(has_mutation = 1)

run_firth_per_mutation <- function(df_period, target_mut) {
  sub_data <- df_period %>%
    mutate(
      y = if_else(group_2class == "Transmitted", 1, 0),
      mut_present = if_else(index %in% (alteration_indexed %>% filter(mutation == target_mut) %>% pull(index)), 1, 0)
    )
  
  if (length(unique(sub_data$y)) < 2 || sum(sub_data$mut_present) == 0) {
    return(tibble(p_value = NA_real_, or = NA_real_))
  }
  
  tryCatch({
    fit <- logistf(y ~ mut_present, data = sub_data, pl = FALSE)
    tibble(
      p_value = fit$prob[2],
      or = exp(fit$coefficients[2])
    )
  }, error = function(e) {
    tibble(p_value = NA_real_, or = NA_real_)
  })
}

mut_period_grid <- hh_sn_alterations %>%
  select(period, gene, mutation) %>%
  distinct()

firth_results_list <- list()
for(i in seq_len(nrow(mut_period_grid))) {
  p_val <- mut_period_grid$period[i]
  m_val <- mut_period_grid$mutation[i]
  g_val <- mut_period_grid$gene[i]
  
  df_p <- day1_classified_periods %>% filter(period == p_val)
  res <- run_firth_per_mutation(df_p, m_val)
  
  firth_results_list[[i]] <- tibble(
    period = p_val,
    gene = g_val,
    mutation = m_val,
    p_value = res$p_value,
    or = res$or
  )
}

firth_period_metrics <- bind_rows(firth_results_list) %>%
  group_by(period) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

significant_transmission_muts <- firth_period_metrics %>%
  filter(!is.na(p_value) & or > 1 & p_value < 0.20)

if (nrow(significant_transmission_muts) == 0) {
  significant_transmission_muts <- firth_period_metrics %>%
    filter(!is.na(or) & is.finite(or)) %>%
    group_by(period) %>%
    slice_max(or, n = 5, with_ties = FALSE) %>%
    ungroup()
}

target_transmission_mutations <- unique(significant_transmission_muts$mutation)

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

# FIX 2: Complete total GISAID sequence totals across all periods
gisaid_totals <- df_gisaid_clean %>%
  group_by(Time_Period) %>%
  summarise(total_seqs = n_distinct(seqName), .groups = "drop") %>%
  complete(Time_Period = all_periods, fill = list(total_seqs = 0))

# ==============================================================================
# 4. Format Summaries & Build Clean Labeling Strings
# ==============================================================================
# FIX 3: Explicit grid expansion for Household dataset
hh_summary <- hh_sn_alterations %>%
  filter(mutation %in% target_transmission_mutations) %>%
  group_by(gene, mutation, period, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  complete(
    nesting(gene, mutation),
    period = factor(all_periods, levels = all_periods),
    group_2class = c("Transmitted", "Non-Transmitted"),
    fill = list(count = 0)
  ) %>%
  left_join(period_group_totals, by = c("period" = "period", "group_2class" = "group_2class")) %>%
  mutate(
    Prevalence = if_else(total_n > 0, (count / total_n) * 100, 0),
    Group = group_2class,
    Time_Period = period
  ) %>%
  select(Time_Period, gene, mutation, Group, count, total_n, Prevalence)

# FIX 4: Explicit grid expansion for GISAID dataset
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
  filter(mutation %in% target_transmission_mutations) %>%
  group_by(Time_Period, gene, mutation) %>%
  summarise(count = n_distinct(seqName), .groups = "drop") %>%
  complete(
    nesting(gene, mutation),
    Time_Period = factor(all_periods, levels = all_periods),
    fill = list(count = 0)
  ) %>%
  left_join(gisaid_totals, by = "Time_Period") %>%
  mutate(
    Prevalence = if_else(total_seqs > 0, (count / total_seqs) * 100, 0),
    Group = "GISAID (Background)",
    Time_Period = factor(Time_Period, levels = all_periods),
    total_n = total_seqs
  ) %>%
  select(Time_Period, gene, mutation, Group, count, total_n, Prevalence)

# FIX 5: Complete combination across ALL mutations, ALL periods, and ALL cohort groups
combined_df <- bind_rows(hh_summary, gisaid_summary) %>%
  complete(
    nesting(gene, mutation),
    Time_Period = factor(all_periods, levels = all_periods),
    Group = factor(all_groups, levels = all_groups),
    fill = list(count = 0, total_n = 0, Prevalence = 0)
  ) %>%
  mutate(
    bar_label = if_else(total_n > 0, sprintf("%.1f%% (%d/%d)", Prevalence, count, total_n), "0.0% (0/0)")
  ) %>%
  left_join(
    firth_period_metrics %>% select(period, mutation, p_value, adj_p_value, or),
    by = c("Time_Period" = "period", "mutation" = "mutation")
  ) %>%
  mutate(
    p_str = if_else(is.na(p_value), "p = N/A", if_else(p_value < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_value)))),
    q_str = if_else(is.na(adj_p_value), "q = N/A", if_else(adj_p_value < 0.001, "q < 0.001", paste0("q = ", sprintf("%.3f", adj_p_value)))),
    or_str = if_else(is.na(or), "OR = N/A", if_else(is.infinite(or), "OR = Inf", paste0("OR = ", sprintf("%.2f", or)))),
    facet_y_label = paste0(mutation, "\n(", p_str, ", ", q_str, ", ", or_str, ")")
  )

# FIX 6: Use drop = FALSE in ggplot scale to preserve missing groups visually
fill_colors <- c(
  "GISAID (Background)" = "#808080",
  "Non-Transmitted"     = "#2B5C8F",
  "Transmitted"         = "#D95F02"
)

fig_firth_gisaid <- ggplot(combined_df, aes(x = Prevalence, y = reorder(facet_y_label, Prevalence), fill = Group)) +
  geom_col(position = position_dodge(width = 0.85, preserve = "single"), width = 0.75, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = bar_label),
    position = position_dodge(width = 0.85, preserve = "single"),
    hjust = -0.1,
    size = 2.4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = fill_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 135), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
  facet_wrap(~ Time_Period, scales = "free_y", ncol = 1) +
  labs(
    title = "Transmission-Associated Mutations (Firth's Regression) vs. GISAID Background",
    subtitle = "Evaluating significant household transmission mutations against background genomic prevalence",
    x = "Prevalence within Group (%)",
    y = NULL,
    fill = "Cohort Group:",
    caption = "Null Hypothesis (H0): Mutation presence does not alter the odds of household transmission status.\nMetrics derived via Firth's penalized logistic regression (Transmitted vs Non-Transmitted)."
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
    panel.spacing = unit(0.8, "lines")
  )

print(fig_firth_gisaid)
