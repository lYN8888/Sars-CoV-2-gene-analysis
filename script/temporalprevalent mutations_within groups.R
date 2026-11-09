library(tidyverse)
library(lubridate)
library(readr)

# ==============================================================================
# 1. Load Household & Lineage Data
# ==============================================================================
lineage_path <- "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/lineage_data(samples 2022-2024).csv"
d35_path     <- "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/household samples/d35_threshold_hhtransmission.csv"
nextclade_path <- "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/household samples/nextclade_all_hhconsensus seq.csv"
gisaid_path  <- "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv"

lineage_data <- read_csv(lineage_path, show_col_types = FALSE) %>%
  select(pid, collection_date, LINEAGE, lineage_group, any_of(c("age", "sex"))) %>%
  mutate(collection_date = mdy(collection_date))

d35 <- read_csv(d35_path, show_col_types = FALSE) %>%
  mutate(
    sampleid     = str_replace(sampleid, "^ht", "HT"),
    household_id = str_sub(sampleid, 3, 5)
  )

nextclade <- read_delim(nextclade_path, delim = ";", show_col_types = FALSE) %>%
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
    group_2class = if_else(day1_group == "Transmitted", "Transmitted", "Non-Transmitted"),
    is_transmitted_bin = if_else(group_2class == "Transmitted", 1, 0)
  )

hh_sn_alterations <- day1_classified_periods %>%
  select(any_of(c("index", "seqName", "group_2class", "is_transmitted_bin", "period", "age", "sex", "aaSubstitutions", "aaDeletions", "aaInsertions"))) %>%
  pivot_longer(
    cols = any_of(c("aaSubstitutions", "aaDeletions", "aaInsertions")),
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

# ==============================================================================
# 2. Period-Specific Statistical Analysis (Transmitted vs Non-Transmitted)
# ==============================================================================
run_robust_stats_for_period <- function(df_period_subset) {
  period_name <- unique(df_period_subset$period)
  
  has_age <- "age" %in% names(df_period_subset) && all(!is.na(df_period_subset$age))
  has_sex <- "sex" %in% names(df_period_subset) && all(!is.na(df_period_subset$sex))
  
  cols_to_keep <- c("index", "is_transmitted_bin", "mutation")
  if (has_age) cols_to_keep <- c(cols_to_keep, "age")
  if (has_sex) cols_to_keep <- c(cols_to_keep, "sex")
  
  mutation_matrix <- df_period_subset %>%
    distinct(across(any_of(cols_to_keep))) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = mutation, values_from = present, values_fill = 0)
  
  meta_cols <- intersect(c("index", "is_transmitted_bin", "age", "sex"), names(mutation_matrix))
  mut_cols <- setdiff(names(mutation_matrix), meta_cols)
  if(length(mut_cols) == 0) return(NULL)
  
  map_dfr(mut_cols, function(mut_name) {
    sub_df <- mutation_matrix %>% select(any_of(c("is_transmitted_bin", "age", "sex")), mut_val = all_of(mut_name))
    
    tab <- table(sub_df$mut_val, sub_df$is_transmitted_bin)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble(Time_Period = period_name, mutation = mut_name, p_value = 1.0, odds_ratio = 1.0, adj_or = NA_real_, regression_p = NA_real_))
    }
    
    fish_res <- tryCatch({ fisher.test(tab) }, error = function(e) { NULL })
    p_val <- if (!is.null(fish_res)) fish_res$p.value else 1.0
    or_val <- if (!is.null(fish_res)) as.numeric(fish_res$estimate) else 1.0
    
    reg_formula <- if (has_age & has_sex) {
      as.formula("is_transmitted_bin ~ mut_val + age + sex")
    } else if (has_age) {
      as.formula("is_transmitted_bin ~ mut_val + age")
    } else {
      as.formula("is_transmitted_bin ~ mut_val")
    }
    
    glm_fit <- tryCatch({ glm(reg_formula, data = sub_df, family = binomial) }, error = function(e) { NULL })
    
    reg_p <- NA_real_
    adjusted_or <- NA_real_
    if (!is.null(glm_fit)) {
      glm_coefs <- summary(glm_fit)$coefficients
      if ("mut_val" %in% rownames(glm_coefs)) {
        reg_p <- glm_coefs["mut_val", "Pr(>|z|)"]
        adjusted_or <- exp(glm_coefs["mut_val", "Estimate"])
      }
    }
    
    tibble(
      Time_Period = period_name,
      mutation = mut_name,
      p_value = p_val,
      odds_ratio = or_val,
      adj_or = adjusted_or,
      regression_p = reg_p
    )
  })
}

fisher_all_mutations <- day1_classified_periods %>%
  filter(!is.na(period)) %>%
  inner_join(hh_sn_alterations %>% select(index, period, mutation), by = c("index", "period")) %>%
  group_by(period) %>%
  group_split() %>%
  map_dfr(run_robust_stats_for_period) %>%
  group_by(Time_Period) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# ==============================================================================
# 3. Period-Specific Prevalence Calculation (Household & GISAID)
# ==============================================================================
period_group_totals <- day1_classified_periods %>%
  filter(!is.na(period)) %>%
  group_by(period, group_2class) %>%
  summarise(total_n_group = n_distinct(index), .groups = "drop")

hh_group_summary <- hh_sn_alterations %>%
  filter(!is.na(period)) %>%
  group_by(period, group_2class, gene, mutation) %>%
  summarise(count = n_distinct(index), .groups = "drop") %>%
  left_join(period_group_totals, by = c("period", "group_2class")) %>%
  mutate(
    Prevalence = (count / total_n_group) * 100,
    Group = group_2class,
    Time_Period = period
  )

df_gisaid_raw <- read_delim(gisaid_path, delim = ";", show_col_types = FALSE)
if (!"aaInsertions" %in% colnames(df_gisaid_raw)) df_gisaid_raw$aaInsertions <- NA_character_

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

gisaid_summary <- df_gisaid_clean %>%
  select(seqName, Time_Period, aaSubstitutions, aaDeletions, aaInsertions) %>%
  pivot_longer(cols = c(aaSubstitutions, aaDeletions, aaInsertions), names_to = "Type", values_to = "Mut_List") %>%
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
  filter(gene %in% c("S", "N")) %>%
  group_by(Time_Period, gene, mutation) %>%
  summarise(count = n_distinct(seqName), .groups = "drop") %>%
  left_join(gisaid_totals, by = "Time_Period") %>%
  mutate(
    Prevalence = (count / total_seqs) * 100,
    Group = "GISAID (Background)",
    total_n_group = total_seqs
  ) %>%
  select(Time_Period, gene, mutation, Group, count, total_n_group, Prevalence)

# ==============================================================================
# 4. Combine & Isolate Truly Period-Specific Top Mutations per Group
# ==============================================================================
combined_all <- bind_rows(
  hh_group_summary %>% select(-any_of("Time_Period")) %>% rename(Time_Period = period), 
  gisaid_summary
) %>%
  mutate(
    Group = factor(Group, levels = c("GISAID (Background)", "Non-Transmitted", "Transmitted")),
    Time_Period = factor(Time_Period, levels = c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024")),
    bar_label = sprintf("%.1f%% (%d/%d)", Prevalence, count, total_n_group)
  ) %>%
  left_join(fisher_all_mutations, by = c("Time_Period", "mutation")) %>%
  mutate(
    stat_label = if_else(
      !is.na(p_value) & p_value < 1.0,
      sprintf("p: %.3f | q: %.3f | OR: %.2f", p_value, adj_p_value, odds_ratio),
      "Background/Baseline"
    ),
    facet_y_label = paste0(mutation)
  )

# Dynamically select the top prevalent mutations unique to each specific epoch and group window
top_mutations_per_epoch <- combined_all %>%
  group_by(Time_Period, Group) %>%
  slice_max(Prevalence, n = 7, with_ties = FALSE) %>%
  ungroup()

# ==============================================================================
# 5. Facet Plot Output with Independent Period Lineage Distributions
# ==============================================================================
fill_colors <- c(
  "GISAID (Background)" = "#808080",
  "Non-Transmitted"      = "#2B5C8F",
  "Transmitted"          = "#D95F02"
)

fig_separated_stats <- ggplot(top_mutations_per_epoch, aes(x = Prevalence, y = reorder(facet_y_label, Prevalence), fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = fill_colors) +
  scale_x_continuous(limits = c(0, 130), expand = c(0, 0)) +
  facet_grid(Time_Period ~ Group, scales = "free", space = "free_y") +
  labs(
    title = "Temporal Prevalent S and N Mutations",
    x = "Prevalence within Group (%)",
    y = NULL,
    fill = "Cohort Group:"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "#EFEFEF", color = "#333333", linewidth = 0.4),
    axis.text.y = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(face = "bold", size = 8),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 8),
    panel.spacing = unit(0.5, "lines")
  )

print(fig_separated_stats)

