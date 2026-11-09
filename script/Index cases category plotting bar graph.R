
library(tidyverse)
library(readr)
d35<- read_csv("~/Covid_study/data/household samples/d35_threshold_hhtransmission.csv")%>%
  mutate(
    sampleid = str_replace(sampleid, "^ht", "HT"),
    household_id = str_sub(sampleid, 3, 5)
  )

nextclade <- read_delim("~/Covid_study/data/household samples/nextclade_all_hhconsensus seq.csv", delim = ";") %>%
  mutate(
    sample_id    = str_extract(seqName, "HT\\d{5}"),
    household_id = str_sub(sample_id, 3, 5),
    is_day1      = str_detect(seqName, "Day1"),
    is_infected  = str_detect(seqName, "Infected")
  )

# --- 2. Identify Household ID Sets ---
# Households present in d35 (transmitted)
d35_hh_ids <- unique(d35$household_id)

# Households with infected cases in Nextclade but NOT present in d35 file
community_infected_hh_ids <- nextclade %>%
  filter(is_infected & !(sample_id %in% d35$sampleid)) %>%
  pull(household_id) %>%
  unique() %>%
  # Exclude households already identified as transmitted in d35
  setdiff(d35_hh_ids)

# --- 3. Filter and Classify Day 1 Cases Only ---
day1_indices_classified <- nextclade %>%
  filter(is_day1) %>%
  mutate(
    day1_category = case_when(
      # Household has infected case(s) in d35
      household_id %in% d35_hh_ids ~ "transmitted_indices",
      
      # Household has infected case(s) outside d35 file
      household_id %in% community_infected_hh_ids ~ "community_infected_indices",
      
      # Household has no infected cases
      TRUE ~ "non_transmitted_indices"
    )
  )

# --- 4. Summary Output ---
# Count of Day 1 cases per category
print("Day 1 Cases Classification Summary:")
table(day1_indices_classified$day1_category)

# Extract row indices for each category (matching Nextclade 'index' column)
transmitted_indices_list <- day1_indices_classified %>% 
  filter(day1_category == "transmitted_indices") %>% 
  pull(index)

community_infected_indices_list <- day1_indices_classified %>% 
  filter(day1_category == "community_infected_indices") %>% 
  pull(index)

non_transmitted_indices_list <- day1_indices_classified %>% 
  filter(day1_category == "non_transmitted_indices") %>% 
  pull(index)

#Plotting

# --- 1. Aggregate Category Counts ---
category_counts <- day1_indices_classified %>%
  count(day1_category) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, " (", round(percentage, 1), "%)"),
    day1_category = factor(
      day1_category,
      levels = c("non_transmitted_indices", "transmitted_indices", "community_infected_indices"),
      labels = c("Non-Transmitted Indices", "Transmitted Indices", "Community Infected Indices")
    )
  )

# --- 2. Create Bar Plot ---
ggplot(category_counts, aes(x = day1_category, y = n, fill = day1_category)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_text(aes(label = label), vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_manual(values = c(
    "Non-Transmitted Indices"    = "#4C72B0",
    "Transmitted Indices"        = "#DD8452",
    "Community Infected Indices" = "#55A868"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Count of the Northern Thailand Household Collection Day 1 Index Cases by Transmission Category",
    x = "Category",
    y = "Number of Day 1 Cases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

