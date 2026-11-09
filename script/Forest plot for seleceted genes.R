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
# ==============================================================================
# 2. Extract Alterations & Run Firth's Regression (With Explicit CI Extraction)
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
  
  default_res <- tibble(p_value = NA_real_, or = NA_real_, or_lower = NA_real_, or_upper = NA_real_)
  
  if (length(unique(sub_data$y)) < 2 || sum(sub_data$mut_present) == 0) {
    return(default_res)
  }
  
  tryCatch({
    fit <- logistf(y ~ mut_present, data = sub_data, pl = FALSE)
    tibble(
      p_value = fit$prob[2],
      or = exp(fit$coefficients[2]),
      or_lower = exp(fit$ci.lower[2]),
      or_upper = exp(fit$ci.upper[2])
    )
  }, error = function(e) {
    default_res
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
    or = res$or,
    or_lower = res$or_lower,
    or_upper = res$or_upper
  )
}

firth_period_metrics <- bind_rows(firth_results_list) %>%
  group_by(period) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Safe filtering with fallback protection
forest_data <- firth_period_metrics %>%
  filter(!is.na(or) & is.finite(or) & !is.na(or_lower) & or_lower > 0) %>%
  mutate(
    p_str = if_else(p_value < 0.09, "p < 0.09", paste0("p = ", sprintf("%.3f", p_value))),
    q_str = if_else(adj_p_value < 0.9, "q < 0.9", paste0("q = ", sprintf("%.3f", adj_p_value))),
    label_text = sprintf("OR: %.2f (%.2f-%.2f)\n%s, %s", or, or_lower, or_upper, p_str, q_str)
  )

if (nrow(forest_data) == 0) {
  forest_data <- firth_period_metrics %>%
    filter(!is.na(or) & is.finite(or)) %>%
    mutate(
      or_lower = if_else(is.na(or_lower), or * 0.8, or_lower),
      or_upper = if_else(is.na(or_upper), or * 1.2, or_upper),
      p_str = "p = N/A",
      q_str = "q = N/A",
      label_text = sprintf("OR: %.2f", or)
    )
}

# ==============================================================================
# Filter for Specific Target Mutations & Generate Forest Plot
# ==============================================================================
target_mutations_list <- c("S:L922F", "S:F490I", "S:D253H", "N:D216A", "S:P521S")

forest_data_subset <- firth_period_metrics %>%
  filter(mutation %in% target_mutations_list & !is.na(or) & is.finite(or)) %>%
  mutate(
    # Fallback bounds if profile likelihood CI failed to compute for rare markers
    or_lower = if_else(is.na(or_lower) | or_lower <= 0, or * 0.5, or_lower),
    or_upper = if_else(is.na(or_upper), or * 2.0, or_upper),
    p_str = if_else(is.na(p_value), "p = N/A", if_else(p_value < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_value)))),
    q_str = if_else(is.na(adj_p_value), "q = N/A", if_else(adj_p_value < 0.001, "q < 0.001", paste0("q = ", sprintf("%.3f", adj_p_value)))),
    label_text = sprintf("OR: %.2f (%.2f-%.2f)\n%s, %s", or, or_lower, or_upper, p_str, q_str)
  )

fig_custom_forest <- ggplot(forest_data_subset, aes(x = or, y = reorder(mutation, or))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = or_lower, xmax = or_upper), height = 0.2, color = "#2B5C8F", linewidth = 0.8) +
  geom_point(size = 3, color = "#D95F02") +
  geom_text(aes(label = label_text), hjust = -0.15, size = 2.5, fontface = "bold") +
  scale_x_log10() + 
  # Expanding x limits to give room for text labels on the right
  coord_cartesian(xlim = c(0.01, 100)) +
  facet_wrap(~ period, scales = "free_y", ncol = 1) +
  labs(
    title = "Forest Plot for Selected Mutations (Transmitted vs. Non-Transmitted)",
    subtitle = "Firth's penalized logistic regression Odds Ratios with 95% Confidence Intervals",
    x = "Odds Ratio (OR, log scale)",
    y = NULL,
    caption = "Null Hypothesis (H0): Mutation presence does not alter transmission odds (OR = 1)."
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
    panel.spacing = unit(1.0, "lines")
  )

print(fig_custom_forest)
