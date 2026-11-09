library(tidyverse)
library(lubridate)
library(readr)
library(logistf)

# Define master factors to enforce consistent grids
all_periods <- c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024")
all_groups  <- c("GISAID (Background)", "Non-Transmitted", "Transmitted")

# --- Load Data & Classify ---
lineage_data <- read_csv("Covid_study/data/lineage_data(samples 2022-2024).csv") %>%
  select(pid, collection_date, LINEAGE, lineage_group) %>%
  mutate(collection_date = mdy(collection_date))

d35 <- read_csv("Covid_study/data/household samples/d35_threshold_hhtransmission.csv") %>%
  mutate(sampleid = str_replace(sampleid, "^ht", "HT"), household_id = str_sub(sampleid, 3, 5))

d35_hh_ids <- unique(d35$household_id)

nextclade <- read_delim("Covid_study/data/household samples/nextclade_all_hhconsensus seq.csv", delim = ";", show_col_types = FALSE) %>%
  mutate(
    sample_id = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1 = str_detect(seqName, "Day1"),
    is_infected = str_detect(seqName, "Infected")
  ) %>%
  left_join(lineage_data, by = c("sample_id" = "pid")) %>%
  mutate(
    period = case_when(
      collection_date >= as.Date("2022-07-01") & collection_date < as.Date("2023-03-01") ~ "Jul 2022 - Mar 2023",
      collection_date >= as.Date("2023-03-01") & collection_date < as.Date("2023-10-01") ~ "Mar 2023 - Oct 2023",
      collection_date >= as.Date("2023-10-01") & collection_date <= as.Date("2024-05-31") ~ "Oct 2023 - May 2024",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = all_periods)
  )

community_infected_hh_ids <- nextclade %>%
  filter(is_infected & !(sample_id %in% d35$sampleid)) %>%
  pull(household_id) %>% unique() %>% setdiff(d35_hh_ids)

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

# --- Extract Alterations & Run Firth Regression ---
hh_sn_alterations <- day1_classified_periods %>%
  select(index, seqName, group_2class, period, aaSubstitutions, aaDeletions, aaInsertions) %>%
  pivot_longer(cols = c(aaSubstitutions, aaDeletions, aaInsertions), names_to = "alteration_type", values_to = "mutation_list") %>%
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

alteration_indexed <- hh_sn_alterations %>% select(index, period, mutation) %>% distinct()

run_firth_per_mutation <- function(df_period, target_mut) {
  sub_data <- df_period %>%
    mutate(
      y = if_else(group_2class == "Transmitted", 1, 0),
      mut_present = if_else(index %in% (alteration_indexed %>% filter(mutation == target_mut) %>% pull(index)), 1, 0)
    )
  if (length(unique(sub_data$y)) < 2 || sum(sub_data$mut_present) == 0) return(tibble(p_value = NA_real_, or = NA_real_))
  tryCatch({
    fit <- logistf(y ~ mut_present, data = sub_data, pl = FALSE)
    tibble(p_value = fit$prob[2], or = exp(fit$coefficients[2]))
  }, error = function(e) tibble(p_value = NA_real_, or = NA_real_))
}

mut_period_grid <- hh_sn_alterations %>% select(period, gene, mutation) %>% distinct()

firth_results_list <- vector("list", nrow(mut_period_grid))
for(i in seq_len(nrow(mut_period_grid))) {
  p_val <- mut_period_grid$period[i]; m_val <- mut_period_grid$mutation[i]; g_val <- mut_period_grid$gene[i]
  df_p <- day1_classified_periods %>% filter(period == p_val)
  res <- run_firth_per_mutation(df_p, m_val)
  firth_results_list[[i]] <- tibble(period = p_val, gene = g_val, mutation = m_val, p_value = res$p_value, or = res$or)
}

firth_period_metrics <- bind_rows(firth_results_list) %>%
  group_by(period) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Select top target mutations to evaluate consistently
target_mutations <- firth_period_metrics %>%
  filter(!is.na(p_value) & or > 1 & p_value < 0.20) %>%
  pull(mutation) %>% unique()

if (length(target_mutations) == 0) {
  target_mutations <- firth_period_metrics %>% filter(!is.na(or)) %>% slice_max(or, n = 5) %>% pull(mutation) %>% unique()
}

# --- Load GISAID ---
df_gisaid_clean <- readr::read_delim("Covid_study/data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv", delim = ";", show_col_types = FALSE) %>%
  mutate(date_clean = ymd(str_extract(seqName, "\\d{4}-\\d{2}-\\d{2}"))) %>%
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

# Compute denominators
period_group_totals <- day1_classified_periods %>%
  group_by(period, group_2class) %>%
  summarise(total_n = n_distinct(index), .groups = "drop") %>%
  complete(period = factor(all_periods, levels = all_periods), group_2class = c("Transmitted", "Non-Transmitted"), fill = list(total_n = 0))

gisaid_totals <- df_gisaid_clean %>%
  group_by(Time_Period) %>%
  summarise(total_seqs = n_distinct(seqName), .groups = "drop") %>%
  complete(Time_Period = all_periods, fill = list(total_seqs = 0))

# --- Build Master Dataset with Uniform Grid Expansion ---
hh_summary <- hh_sn_alterations %>%
  filter(mutation %in% target_mutations) %>%
  group_by(gene, mutation, period, group_2class) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  complete(nesting(gene, mutation), period = factor(all_periods, levels = all_periods), group_2class = c("Transmitted", "Non-Transmitted"), fill = list(count = 0)) %>%
  left_join(period_group_totals, by = c("period", "group_2class")) %>%
  mutate(Prevalence = if_else(total_n > 0, (count / total_n) * 100, 0), Group = group_2class, Time_Period = period)

gisaid_summary <- df_gisaid_clean %>%
  select(seqName, Time_Period, aaSubstitutions, aaDeletions, aaInsertions) %>%
  pivot_longer(cols = starts_with("aa"), names_to = "Type", values_to = "Mut_List") %>%
  filter(!is.na(Mut_List) & Mut_List != "" & Mut_List != "None") %>%
  separate_rows(Mut_List, sep = ",") %>%
  mutate(
    raw_mut = trimws(Mut_List), gene = sub(":.*", "", raw_mut),
    mutation = case_when(
      Type == "aaDeletions" & !str_detect(raw_mut, "(del|Deletion)") ~ paste0(raw_mut, " (del)"),
      Type == "aaInsertions" & !str_detect(raw_mut, "(ins|Insertion)") ~ paste0(raw_mut, " (ins)"),
      TRUE ~ raw_mut
    )
  ) %>%
  filter(mutation %in% target_mutations) %>%
  group_by(Time_Period, gene, mutation) %>%
  summarise(count = n_distinct(seqName), .groups = "drop") %>%
  complete(nesting(gene, mutation), Time_Period = factor(all_periods, levels = all_periods), fill = list(count = 0)) %>%
  left_join(gisaid_totals, by = "Time_Period") %>%
  mutate(Prevalence = if_else(total_seqs > 0, (count / total_seqs) * 100, 0), Group = "GISAID (Background)", Time_Period = factor(Time_Period, levels = all_periods), total_n = total_seqs)

# --- Build Master Dataset with p-value, q-value, and OR Strings ---
master_combined <- bind_rows(hh_summary, gisaid_summary) %>%
  left_join(
    firth_period_metrics %>% select(period, mutation, p_value, adj_p_value, or), 
    by = c("Time_Period" = "period", "mutation" = "mutation")
  ) %>%
  mutate(
    sig_asterisk = if_else(!is.na(adj_p_value) & adj_p_value < 0.05, " *", ""),
    bar_label = if_else(total_n > 0, sprintf("%.1f%% (%d/%d)", Prevalence, count, total_n), "0.0% (0/0)"),
    
    # Format statistical strings cleanly
    p_str  = if_else(is.na(p_value), "p = N/A", if_else(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))),
    q_str  = if_else(is.na(adj_p_value), "q = N/A", if_else(adj_p_value < 0.001, "q < 0.001", sprintf("q = %.3f", adj_p_value))),
    or_str = if_else(is.na(or), "OR = N/A", if_else(is.infinite(or), "OR = Inf", sprintf("OR = %.2f", or))),
    
    # Combined y-axis facet label with p, q, and OR
    facet_y_label = paste0(mutation, "\n(", p_str, ", ", q_str, ", ", or_str, ")")
  )

# --- Reusable Plotting Function ---
plot_cohort_comparison <- function(data, selected_groups, fill_colors, title_text, null_hypothesis_text) {
  df_sub <- data %>%
    filter(Group %in% selected_groups) %>%
    mutate(Group = factor(Group, levels = selected_groups))
  
  # Combine Null Hypothesis and statistical metric legend in the caption
  full_caption <- paste0(
    "Null Hypothesis (H0): ", null_hypothesis_text, "\n",
    "Metrics: p = Firth's regression p-value | q = BH adjusted p-value (FDR) | OR = Odds Ratio (Transmitted vs Non-Transmitted)"
  )
  
  ggplot(df_sub, aes(x = Prevalence, y = reorder(facet_y_label, Prevalence), fill = Group)) +
    geom_col(position = position_dodge(width = 0.85, preserve = "single"), width = 0.75, color = "black", linewidth = 0.2) +
    geom_text(
      aes(label = bar_label), 
      position = position_dodge(width = 0.85, preserve = "single"), 
      hjust = -0.1, 
      size = 2.2, 
      fontface = "bold"
    ) +
    scale_fill_manual(values = fill_colors, drop = FALSE) +
    scale_x_continuous(limits = c(0, 145), expand = c(0, 0), breaks = c(0, 25, 50, 75, 100)) +
    facet_wrap(~ Time_Period, scales = "free_y", ncol = 1) +
    labs(
      title = title_text, 
      x = "Prevalence within Group (%)", 
      y = NULL, 
      fill = "Cohort Group:",
      caption = full_caption
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 7.5, color = "#333333", hjust = 0, lineheight = 1.2, face = "italic"),
      strip.text = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.4),
      axis.text.y = element_text(size = 7.5, face = "bold"),
      legend.position = "top"
    )
}

# ==============================================================================
# Generate the 3 Figures with Specific Null Hypotheses
# ==============================================================================

# --- Plot 1: Transmitted vs Non-Transmitted ---
fig1 <- plot_cohort_comparison(
  master_combined,
  selected_groups = c("Non-Transmitted", "Transmitted"),
  fill_colors = c("Non-Transmitted" = "#2B5C8F", "Transmitted" = "#D95F02"),
  title_text = "Household Cohort: Transmitted vs. Non-Transmitted",
  null_hypothesis_text = "Mutation presence does not alter the odds of household transmission status (OR = 1)."
)

# --- Plot 2: Transmitted vs GISAID ---
fig2 <- plot_cohort_comparison(
  master_combined,
  selected_groups = c("GISAID (Background)", "Transmitted"),
  fill_colors = c("GISAID (Background)" = "#808080", "Transmitted" = "#D95F02"),
  title_text = "Transmitted Cohort vs. GISAID Background",
  null_hypothesis_text = "Mutation prevalence in transmitted household cases equals community background prevalence."
)

# --- Plot 3: Non-Transmitted vs GISAID ---
fig3 <- plot_cohort_comparison(
  master_combined,
  selected_groups = c("GISAID (Background)", "Non-Transmitted"),
  fill_colors = c("GISAID (Background)" = "#808080", "Non-Transmitted" = "#2B5C8F"),
  title_text = "Non-Transmitted Cohort vs. GISAID Background",
  null_hypothesis_text = "Mutation prevalence in non-transmitted household cases equals community background prevalence."
)

# Display the plots
print(fig1)
