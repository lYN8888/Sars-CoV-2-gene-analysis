library(dplyr)
library(ggplot2)
library(scales)
library(readr)

# ==============================================================================
# 1. Load Household & Lineage Data
# ==============================================================================
lineage_data <- read_csv("data/lineage_data(samples 2022-2024).csv")

monthly_samples <- lineage_data %>%
  # Parse single slashes with Month first (%m), Year second (%Y), and Day last (%d)
  mutate(collection_date = as.Date(collection_date, format = "%m/%d/%Y")) %>%
  
  # The rest of the pipeline handles the plotting aggregation cleanly:
  filter(!is.na(collection_date)) %>%
  mutate(month_year = as.Date(format(collection_date, "%Y-%m-01"))) %>%
  group_by(month_year) %>%
  summarise(sample_size = n_distinct(pid), .groups = "drop")

ggplot(monthly_samples, aes(x = month_year, y = sample_size)) +
  geom_line(color = "#1d3557", linewidth = 1.2) +
  geom_point(color = "#e63946", size = 2.5) +
  theme_minimal(base_size = 14) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Northern thailand SARS Cov-2 Sample collection size distribution over time  ",
    x = "Collection Date",
    y = "Number of Samples"
  )
