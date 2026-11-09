library(dplyr)
library(ggplot2)
library(scales)
library(readr)

# ==============================================================================
# 1. Load Household & Lineage Data
# ==============================================================================
lineage_data <- read_csv("~/Covid_study/data/lineage_data(samples 2022-2024).csv")

monthly_samples <- lineage_data %>%
  # Parse single slashes with Month first (%m), Year second (%Y), and Day last (%d)
  mutate(collection_date = as.Date(collection_date, format = "%m/%d/%Y")) %>%
  filter(!is.na(collection_date)) %>%
  mutate(month_year = as.Date(format(collection_date, "%Y-%m-01"))) %>%
  group_by(month_year) %>%
  summarise(sample_size = n_distinct(pid), .groups = "drop")

# Define period boundary dates for vertical lines
cutoff_dates <- as.Date(c("2023-03-01", "2023-10-01"))

# ==============================================================================
# 2. Plot with Separating Lines and Period Annotations
# ==============================================================================
ggplot(monthly_samples, aes(x = month_year, y = sample_size)) +
  
  # --- Alternating Shaded Background Regions ---
  annotate("rect", xmin = as.Date("2022-07-01"), xmax = as.Date("2023-03-01"),
           ymin = -Inf, ymax = Inf, fill = "#1d3557", alpha = 0.05) +
  annotate("rect", xmin = as.Date("2023-03-01"), xmax = as.Date("2023-10-01"),
           ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.05) +
  annotate("rect", xmin = as.Date("2023-10-01"), xmax = as.Date("2024-05-31"),
           ymin = -Inf, ymax = Inf, fill = "#1d3557", alpha = 0.05) +
  
  # --- Separating Vertical Lines ---
  geom_vline(xintercept = cutoff_dates, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  
  # --- Period Text Labels ---
  annotate("text", x = as.Date("2022-11-01"), y = 20, label = "Jul 2022 - Mar 2023", 
           fontface = "bold", size = 3.5, color = "#5d3557") +
  annotate("text", x = as.Date("2023-06-15"), y = 20, label = "Mar 2023 - Oct 2023", 
           fontface = "bold", size = 3.5, color = "#e63946") +
  annotate("text", x = as.Date("2024-01-15"), y = 20, label = "Oct 2023 - May 2024", 
           fontface = "bold", size = 3.5, color = "#2d1000") +
  
  # --- Line & Point Layers ---
  geom_line(color = "#1d3557", linewidth = 1.2) +
  geom_point(color = "#e63946", size = 2.5) +
  
  # --- Formatting & Axes ---
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y",
               limits = c(as.Date("2022-07-01"), as.Date("2024-05-31")),
               expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 22)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Northern Thailand Household SARS-CoV-2 Sample Collection Size Distribution across three periods",
    x = "Collection Date",
    y = "Number of Samples"
  )

