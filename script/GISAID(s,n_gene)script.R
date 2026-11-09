library(tidyverse)
library(lubridate)
library(readr)
# 1. Read correctly using read.csv with sep = ";"
metaGISAID22_24_nxtclde <- readr::read_delim(
  "data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv",
  delim = ";"
)
library(tidyverse)

# 1. Calculate total number of unique sequences in your dataset
total_sequences <- n_distinct(metaGISAID22_24_nxtclde$seqName)

# 2. Filter mutations that appear in at least 5% of sequences
# (Change 0.05 to 0.10 if you prefer a stricter 10% cutoff)
mutation_counts <- gene_data %>%
  count(Gene, Mutation, sort = TRUE) %>%
  mutate(Prevalence = (n / total_sequences) * 100) %>%
  filter(Prevalence >= 5) # Cutoff: >= 5% prevalence

# 3. Plot the justifiable mutation set
ggplot(mutation_counts, aes(x = reorder(Mutation, n), y = n, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Gene, scales = "free") +
  scale_fill_manual(values = c("S" = "#2B5C8F", "N" = "#D95F02")) +
  labs(
    title = "GISAID data S and N Gene Substitutions Present in ≥ 5% of Sequences from samples collected during 2022-2024 ",
    subtitle = paste("Prevalence ≥ 5% across N =", total_sequences, "sequences"),
    x = "Substitution Locus",
    y = "Total sample count with distinct(≥ 5% of total sequences) substitutions"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 9, face = "bold")
  )

# 2. Extract, separate, and filter S and N gene mutations
gene_data <- metaGISAID22_24_nxtclde %>%
  select(seqName, aaSubstitutions) %>%
  filter(!is.na(aaSubstitutions) & aaSubstitutions != "") %>%
  separate_rows(aaSubstitutions, sep = ",") %>%
  mutate(Mutation = trimws(aaSubstitutions)) %>%
  mutate(Gene = sub(":.*", "", Mutation)) %>%
  filter(Gene %in% c("S", "N"))

# 3. View output
head(gene_data)
mutation_counts <- gene_data %>%
  count(Gene, Mutation, sort = TRUE) %>%
  group_by(Gene) %>%
  slice_max(n, n = 15, with_ties = FALSE)

library(tidyverse)

# 1. Total sequences count
total_sequences <- n_distinct(metaGISAID22_24_nxtclde$seqName)

# 2. Extract BOTH aaSubstitutions AND aaDeletions
gene_data_all <- metaGISAID22_24_nxtclde %>%
  select(seqName, aaSubstitutions, aaDeletions) %>%
  # Pivot both columns into a single long column of mutations
  pivot_longer(
    cols = c(aaSubstitutions, aaDeletions), 
    names_to = "Mutation_Type", 
    values_to = "Mutation_List"
  ) %>%
  filter(!is.na(Mutation_List) & Mutation_List != "") %>%
  separate_rows(Mutation_List, sep = ",") %>%
  mutate(Mutation = trimws(Mutation_List)) %>%
  mutate(Gene = sub(":.*", "", Mutation)) %>%
  filter(Gene %in% c("S", "N"))

# 3. Count frequencies
mutation_counts_all <- gene_data_all %>%
  count(Gene, Mutation, sort = TRUE) %>%
  mutate(Prevalence = (n / total_sequences) * 100) %>%
  filter(Prevalence >= 5)

# 4. Plot full landscape (Substitutions + Deletions)
ggplot(mutation_counts_all, aes(x = reorder(Mutation, n), y = n, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Gene, scales = "free") +
  scale_fill_manual(values = c("S" = "#2B5C8F", "N" = "#D95F02")) +
  labs(
    title = "S and N Gene Amino Acid Alterations (Substitutions & Deletions) of GISAID samples collected during 2022-2024 ",
    subtitle = paste("Prevalence ≥ 5% across N =", total_sequences, "sequences"),
    x = "Mutation Locus",
    y = "Total count of the samples with distinct(>5%) mutations"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 9, face = "bold")
  )
library(tidyverse)

# 1. Total sequences count
total_sequences <- n_distinct(metaGISAID22_24_nxtclde$seqName)

# 2. Extract ONLY aaDeletions
deletions_data <- metaGISAID22_24_nxtclde %>%
  select(seqName, aaDeletions) %>%
  filter(!is.na(aaDeletions) & aaDeletions != "") %>%
  separate_rows(aaDeletions, sep = ",") %>%
  mutate(Mutation = trimws(aaDeletions)) %>%
  mutate(Gene = sub(":.*", "", Mutation)) %>%
  filter(Gene %in% c("S", "N"))

# 3. Count frequencies and apply the >= 5% prevalence threshold
deletion_counts <- deletions_data %>%
  count(Gene, Mutation, sort = TRUE) %>%
  mutate(Prevalence = (n / total_sequences) * 100) %>%
  filter(Prevalence >= 5)

# 4. Plot ONLY deletions
ggplot(deletion_counts, aes(x = reorder(Mutation, n), y = n, fill = Gene)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~Gene, scales = "free") +
  scale_fill_manual(values = c("S" = "#2B5C8F", "N" = "#D95F02")) +
  labs(
    title = "GISAID data S and N Gene Deletions Present in ≥ 5% of Sequences from samples collected during 2022-2024",
    subtitle = paste("Prevalence ≥ 5% across N =", total_sequences, "sequences"),
    x = "Deletion Locus",
    y = "Total count of the samples with distinct(>5%) deletions"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 9, face = "bold")
  )

