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
# 3. Process Alterations (Substitutions + Deletions)
# ------------------------------------------------------------------------------
gene_mutations <- df_clean %>%
  select(seqName, Time_Period, aaSubstitutions, aaDeletions) %>%
  pivot_longer(
    cols = c(aaSubstitutions, aaDeletions),
    names_to = "Type",
    values_to = "Mutation_List"
  ) %>%
  filter(!is.na(Mutation_List) & Mutation_List != "") %>%
  separate_rows(Mutation_List, sep = ",") %>%
  mutate(
    Mutation = trimws(Mutation_List),
    Gene = sub(":.*", "", Mutation)
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
    title = "S and N Gene Alterations (≥5% Prevalence) Across Epidemic Windows",
    subtitle = "Extracted from GISAID seqName dates (July 2022 – May 2024)",
    x = "Mutation / Deletion Locus",
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






library(tidyverse)
library(readr)

# --- Step 1: Load and Clean Datasets ---
spec(nextclade)

d35 <-  read_csv("data/household samples/d35_threshold_hhtransmission.csv") %>%
  mutate(
    sampleid = str_replace(sampleid, "^ht", "HT"),
    household_id = str_sub(sampleid, 3, 5)
  )

nextclade <- read_delim("data/household samples/nextclade_all_hhconsensus seq.csv", delim = ";", show_col_types= FALSE) %>%
  mutate(
    sample_id    = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1      = str_detect(seqName, "Day1"),
    is_infected  = str_detect(seqName, "Infected")
  )

# --- Step 2: Define Household ID Categories ---
d35_hh_ids <- unique(d35$household_id)

community_infected_hh_ids <- nextclade %>%
  filter(is_infected & !(sample_id %in% d35$sampleid)) %>%
  pull(household_id) %>%
  unique() %>%
  setdiff(d35_hh_ids)

# --- Step 3: Classify Day 1 Cases (3-Group Classification) ---
day1_classified <- nextclade %>%
  filter(is_day1) %>%
  mutate(
    day1_group = case_when(
      household_id %in% d35_hh_ids ~ "Transmitted",
      household_id %in% community_infected_hh_ids ~ "Community Infected",
      TRUE ~ "Non-Transmitted"
    )
  )

# --- Step 4: Extract and Separate S and N Gene Mutations ---
day1_sn_mutations <- day1_classified %>%
  # Keep relevant columns
  select(index, seqName, household_id, day1_group, aaSubstitutions) %>%
  # Split comma-separated mutations into individual rows
  separate_rows(aaSubstitutions, sep = ",") %>%
  # Extract Gene Prefix and Mutation Name
  mutate(
    gene = str_extract(aaSubstitutions, "^[^:]+"),             # Extracts 'S', 'N', 'ORF1a', etc.
    mutation = str_remove(aaSubstitutions, "^[^:]+:")          # Extracts the mutation e.g., 'T21R'
  ) %>%
  # Filter ONLY for S and N genes
  filter(gene %in% c("S", "N"))

# View sample of parsed S/N mutations
head(day1_sn_mutations)


# --- Step 5: Summarize Most Frequent S/N Mutations per Group ---
sn_mutation_summary <- day1_sn_mutations %>%
  group_by(day1_group, gene, aaSubstitutions) %>%
  summarise(
    count = n(),
    sample_percentage = round((n() / n_distinct(index)) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(day1_group, gene, desc(count))

print(sn_mutation_summary)


# --- Step 6: Extract Wide Format Matrix (Samples vs S/N Mutations) ---
# Useful if you want a binary matrix (0/1 presence) of mutations across samples
sn_binary_matrix <- day1_sn_mutations %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = c(seqName, day1_group),
    names_from = aaSubstitutions,
    values_from = present,
    values_fill = 0
  )

head(sn_binary_matrix)

# Color Palette for Minimalistic Styling
palette_colors <- c(
  "Non-Transmitted"    = "#A8DADC", # Soft Sage / Light Teal
  "Transmitted"        = "#D8B4F8", # Soft Lavender
  "Community Infected" = "#FFB3A7"  # Soft Coral
)

# --- Plot 1: Overall Group Counts (Left Panel) ---
plot_counts <- day1_classified %>%
  count(day1_group) %>%
  mutate(
    percentage = n / sum(n) * 100,
    pct_label  = paste0(round(percentage, 1), "%")
  )

p1 <- ggplot(plot_counts, aes(x = day1_group, y = n, fill = day1_group)) +
  geom_bar(stat = "identity", width = 0.5, color = NA) +
  geom_text(aes(label = n), vjust = -0.5, fontface = "bold", size = 4.5, color = "#222222") +
  geom_text(aes(label = pct_label, y = n / 2), fontface = "bold", size = 4, color = "#222222") +
  scale_fill_manual(values = palette_colors) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  labs(
    title = "Overall Category Counts",
    x = NULL,
    y = "Number of Day 1 Cases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "#222222"),
    axis.text.x = element_text(face = "bold", color = "#333333"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E0E0E0", linetype = "dashed")
  )

# --- Step 3: Extract Top S & N Mutations ---
day1_sn_mutations <- day1_classified %>%
  select(index, seqName, day1_group, aaSubstitutions) %>%
  separate_rows(aaSubstitutions, sep = ",") %>%
  mutate(
    gene = str_extract(aaSubstitutions, "^[^:]+"),
    mutation = aaSubstitutions
  ) %>%
  filter(gene %in% c("S", "N"))

top_mutations <- day1_sn_mutations %>%
  group_by(day1_group, mutation) %>%
  summarise(
    count = n(),
    group_total = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(percentage = round((count / group_total) * 100, 1)) %>%
  group_by(day1_group) %>%
  slice_max(order_by = percentage,n= 50, with_ties = FALSE) %>%
  ungroup()

# Function to generate sub-plots for top mutations per group
make_mutation_subplot <- function(group_name, group_color) {
  df <- top_mutations %>%
    filter(day1_group == group_name)
  
  ggplot(df, aes(x = percentage, y = reorder(mutation, percentage))) +
    geom_bar(stat = "identity", fill = group_color, width = 0.55) +
    geom_text(aes(label = paste0(round(percentage), "%")), hjust = -0.2, fontface = "bold", size = 3.5, color = "#333333") +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0)) +
    labs(
      title = toupper(group_name),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10, color = "#333333", hjust = 0),
      axis.text.y = element_text(face = "bold", color = "#444444"),
      axis.text.x = element_blank(),
      panel.grid = element_blank()
    )
}
install.packages("patchwork")
library(patchwork)
library(tidyverse)

# 2. Define the palette
palette_colors <- c(
  "Non-Transmitted"    = "#A8DADC",
  "Transmitted"        = "#D8B4F8",
  "Community Infected" = "#FFB3A7"
)

# 3. Re-define the subplots
p_non   <- make_mutation_subplot("Non-Transmitted", palette_colors["Non-Transmitted"])
p_trans <- make_mutation_subplot("Transmitted", palette_colors["Transmitted"])
p_comm  <- make_mutation_subplot("Community Infected", palette_colors["Community Infected"])

# 4. Combine vertically (this will work now)
p_mutations <- (p_non / p_trans / p_comm) + 
  plot_annotation(
    title = "Top S and N Gene Amino Acid Mutations", 
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5, color = "#222222"))
  )

# 5. Combine with the left plot side-by-side
final_plot <- (p1 | p_mutations) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "Spike (S) and Nucleocapsid (N) Gene Mutations across Day 1 Index Groups",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#111111"))
  )

print(final_plot)

library(tidyverse)
library(patchwork)
library(lubridate)

# --- Step 1: Read and Clean Datasets ---
d35 <- read_csv("data/household samples/d35_threshold_hhtransmission.csv") %>%
  mutate(
    sampleid     = str_replace(sampleid, "^ht", "HT"),
    household_id = str_sub(sampleid, 3, 5)
  )

nextclade <- read_delim("data/household samples/nextclade_all_hhconsensus seq.csv", delim = ";", show_col_types = FALSE) %>%
  mutate(
    sample_id    = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1      = str_detect(seqName, "Day1"),
    is_infected  = str_detect(seqName, "Infected"),
    
    # Extract date (assumes format YYYY-MM-DD or standard date pattern in seqName/date column)
    # Adjust regex/column name if collection date is stored in a dedicated column like 'date'
    collection_date = as.Date(str_extract(seqName, "\\d{4}-\\d{2}-\\d{2}")), 
    
    # Define Time Periods
    period = case_when(
      collection_date >= as.Date("2022-07-01") & collection_date < as.Date("2023-03-01") ~ "Jul 2022 - Mar 2023",
      collection_date >= as.Date("2023-03-01") & collection_date < as.Date("2023-10-01") ~ "Mar 2023 - Oct 2023",
      collection_date >= as.Date("2023-10-01") & collection_date <= as.Date("2024-05-31") ~ "Oct 2023 - May 2024",
      TRUE ~ "Other / Missing Date"
    ),
    period = factor(period, levels = c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024"))
  )
library(tidyverse)
library(patchwork)

# --- Read Data and Assign Period by Clade ---
nextclade <- read_delim(
  "data/household samples/nextclade_all_hhconsensus seq.csv", 
  delim = ";", 
  show_col_types = FALSE
) %>%
  mutate(
    sample_id    = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1      = str_detect(seqName, "Day1"),
    is_infected  = str_detect(seqName, "Infected"),
    
    # Classify time period by clade designation
    period = case_when(
      str_detect(clade, "21|22") ~ "Jul 2022 - Mar 2023",
      str_detect(clade, "23|XBB") ~ "Mar 2023 - Oct 2023",
      str_detect(clade, "24|JN")  ~ "Oct 2023 - May 2024",
      TRUE ~ "Jul 2022 - Mar 2023" # Fallback if unclassified
    ),
    period = factor(period, levels = c("Jul 2022 - Mar 2023", "Mar 2023 - Oct 2023", "Oct 2023 - May 2024"))
  )
# --- Step 2: Define Household ID Categories & Classify Day 1 Cases ---
d35_hh_ids <- unique(d35$household_id)

community_infected_hh_ids <- nextclade %>%
  filter(is_infected & !(sample_id %in% d35$sampleid)) %>%
  pull(household_id) %>%
  unique() %>%
  setdiff(d35_hh_ids)

day1_classified <- nextclade %>%
  filter(is_day1, !is.na(period), period != "Other / Missing Date") %>%
  mutate(
    day1_group = case_when(
      household_id %in% d35_hh_ids ~ "Transmitted",
      household_id %in% community_infected_hh_ids ~ "Community Infected",
      TRUE ~ "Non-Transmitted"
    )
  )

# --- Step 3: Combine Substitutions and Deletions for S and N Genes ---
# Gather aaSubstitutions and aaDeletions
day1_sn_mutations <- day1_classified %>%
  select(index, seqName, household_id, day1_group, period, aaSubstitutions, aaDeletions) %>%
  # Pivot columns long to treat substitutions and deletions together
  pivot_longer(
    cols = c(aaSubstitutions, aaDeletions),
    names_to = "mutation_type",
    values_to = "mutations_str"
  ) %>%
  filter(!is.na(mutations_str) & mutations_str != "") %>%
  separate_rows(mutations_str, sep = ",") %>%
  mutate(
    gene = str_extract(mutations_str, "^[^:]+"),
    mutation = mutations_str
  ) %>%
  filter(gene %in% c("S", "N"))

# --- Step 4: Summarize Frequencies per Group and Period ---
sn_summary_by_period <- day1_sn_mutations %>%
  group_by(period, day1_group, gene, mutation) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  # Calculate sample counts per (period x group) for accurate percentage denominator
  left_join(
    day1_classified %>%
      group_by(period, day1_group) %>%
      summarise(total_samples = n_distinct(index), .groups = "drop"),
    by = c("period", "day1_group")
  ) %>%
  mutate(percentage = round((count / total_samples) * 100, 1))

# --- Step 5: Filter Top Mutations per Period and Group for Visualization ---
top_sn_period <- sn_summary_by_period %>%
  group_by(period, day1_group) %>%
  slice_max(order_by = percentage, n = 15, with_ties = FALSE) %>%
  ungroup()

# --- Step 6: Plot Visualization ---
palette_colors <- c(
  "Non-Transmitted"    = "#A8DADC",
  "Transmitted"        = "#D8B4F8",
  "Community Infected" = "#FFB3A7"
)

plot_sn_by_period <- ggplot(top_sn_period, aes(x = percentage, y = reorder(mutation, percentage), fill = day1_group)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(round(percentage), "%")), 
    hjust = -0.15, size = 2.8, fontface = "bold", color = "#333333"
  ) +
  scale_x_continuous(limits = c(0, 115), expand = c(0, 0)) +
  scale_fill_manual(values = palette_colors) +
  facet_grid(period ~ day1_group, scales = "free_y", space = "free_y") +
  labs(
    title = "S and N Gene Mutations (Substitutions & Deletions) by Time Period",
    subtitle = "Top mutations shown per group across three epidemiologic windows",
    x = "Prevalence within Group (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "#555555"),
    strip.text.x     = element_text(face = "bold", size = 11),
    strip.text.y     = element_text(face = "bold", size = 10, angle = 270),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    axis.text.y      = element_text(size = 8.5, color = "#222222"),
    axis.text.x      = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.8, "lines")
  )

print(plot_sn_by_period)

