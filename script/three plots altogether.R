library(tidyverse)
library(scales)
library(readr)
library(lubridate)
library(stringr)
library(patchwork)

# Shared Period Boundaries
period_lines <- as.Date(c("2023-03-01", "2023-10-01"))

# ==============================================================================
# 1. FIXED UNIFIED FACTOR LEVELS & COLOR PALETTE
# ==============================================================================
lineage_order <- c("BA.2.75+", "BA.4/BA.5", "XBB", "JN.1+", "Alpha", "Delta", "Recombinant", "Other")

shared_colors <- c(
  "BA.2.75+"    = "#2171B5",
  "BA.4/BA.5"   = "#084594",
  "XBB"         = "#41AB5D",
  "JN.1+"       = "#E87F24",
  "Alpha"       = "#587F24",
  "Delta"       = "#004B87",
  "Recombinant" = "#9E9AC8",
  "Other"       = "#BAB0AC"
)

# ==============================================================================
# FIGURE A: Household Cohort
# ==============================================================================
lineage_data_hh <- read_csv("~/Covid_study/data/lineage_data(samples 2022-2024).csv") %>%
  mutate(
    collection_date = as.Date(collection_date, format = "%m/%d/%Y"),
    month = as.Date(month, format = "%m/%d/%Y"),
    lineage_group = case_when(
      lineage_group == "BA.2.75+/ BA.2.9+/ BA.2.10+" ~ "BA.2.75+",
      grepl("BA.4/BA.5", lineage_group)             ~ "BA.4/BA.5",
      grepl("XBB", lineage_group)                   ~ "XBB",
      grepl("JN.1", lineage_group)                  ~ "JN.1+",
      TRUE                                          ~ lineage_group
    )
  )

area_hh <- lineage_data_hh %>%
  filter(!is.na(month), lineage_group != "Unassigned") %>%
  count(month, lineage_group) %>%
  complete(month, lineage_group = lineage_order, fill = list(n = 0)) %>%
  group_by(month) %>%
  mutate(
    proportion = n / sum(n),
    lineage_group = factor(lineage_group, levels = lineage_order)
  ) %>%
  ungroup()

p1 <- ggplot(area_hh, aes(x = month, y = proportion, fill = lineage_group, group = lineage_group)) +
  geom_area(position = "fill", color = "white", linewidth = 0.2, alpha = 0.92) +
  geom_vline(xintercept = period_lines, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = as.Date("2022-11-01"), y = 1.03, label = "Period 1:\nJuly 2022 - Mar 2023", size = 2.8, fontface = "bold", vjust = 0) +
  annotate("text", x = as.Date("2023-06-15"), y = 1.03, label = "Period 2:\nMar 2023 - Oct 2023", size = 2.8, fontface = "bold", vjust = 0) +
  annotate("text", x = as.Date("2024-01-15"), y = 1.03, label = "Period 3:\nOct 2023 - May 2024", size = 2.8, fontface = "bold", vjust = 0) +
  scale_fill_manual(values = shared_colors, breaks = lineage_order, limits = lineage_order, name = NULL, drop = FALSE) +
  scale_x_date(limits = c(as.Date("2022-07-01"), as.Date("2024-05-31")), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12)), breaks = c(0, .25, .50, .75, 1)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(color = "grey70", linewidth = 0.3),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(t = 15, r = 10, b = 2, l = 10)
  ) +
  labs(
    title = "A. Percentage of SARS-CoV-2 Lineages (Chiang Mai Household Cohort)",
    y = "Percentage"
  )

# ==============================================================================
# FIGURE B: Sample Size Distribution
# ==============================================================================
monthly_samples <- lineage_data_hh %>%
  filter(!is.na(collection_date)) %>%
  mutate(month_year = as.Date(format(collection_date, "%Y-%m-01"))) %>%
  group_by(month_year) %>%
  summarise(sample_size = n_distinct(pid), .groups = "drop")

p2 <- ggplot(monthly_samples, aes(x = month_year, y = sample_size)) +
  annotate("rect", xmin = as.Date("2022-07-01"), xmax = as.Date("2023-03-01"), ymin = -Inf, ymax = Inf, fill = "#1d3557", alpha = 0.05) +
  annotate("rect", xmin = as.Date("2023-03-01"), xmax = as.Date("2023-10-01"), ymin = -Inf, ymax = Inf, fill = "#e63946", alpha = 0.05) +
  annotate("rect", xmin = as.Date("2023-10-01"), xmax = as.Date("2024-05-31"), ymin = -Inf, ymax = Inf, fill = "#1d3557", alpha = 0.05) +
  geom_vline(xintercept = period_lines, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_line(color = "#1d3557", linewidth = 1) +
  geom_point(color = "#e63946", size = 2) +
  scale_x_date(limits = c(as.Date("2022-07-01"), as.Date("2024-05-31")), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 22)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(t = 5, r = 10, b = 2, l = 10)
  ) +
  labs(
    title = "B. Northern Thailand Household Sample Collection Size Across Three Periods",
    y = "Sample Count"
  )

# ==============================================================================
# FIGURE C: GISAID Cohort
# ==============================================================================
df_gisaid_raw <- readr::read_delim(
  "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv",
  delim = ";",
  show_col_types = FALSE
)

assign_lineage_group_unified <- function(lineage) {
  case_when(
    grepl("^JN\\.|^JN1|^JN\\.1|^JN11|^JN14|^KP\\.|^KB\\.|^LB\\.|^JF\\.|^JS\\.", lineage) ~ "JN.1+",
    grepl("^XBB|^XBL\\.|^XBU\\.|^EG\\.|^GJ\\.|^GE\\.|^CH\\.|^BM\\.|^BQ\\.|^BL\\.|^BF\\.|^FL\\.|^FD\\.|^HK\\.|^LF\\.|^HH\\.|^FE\\.|^EU\\.|^GM\\.|^BR\\.|^BY\\.|^GR\\.|^HN\\.", lineage) ~ "XBB",
    grepl("^BA\\.4|^BA\\.5|^BE\\.|^BF\\.|^LZ|^FP\\.", lineage) ~ "BA.4/BA.5",
    grepl("^B\\.1\\.1\\.529|^BA\\.1|^BA\\.2|^EJ\\.|^BN\\.|^BS\\.", lineage) ~ "BA.2.75+",
    grepl("^AY\\.|^GY\\.|^B\\.1\\.617|^B\\.1\\.617\\.2", lineage) ~ "Delta",
    grepl("^B\\.1\\.1\\.7", lineage) ~ "Alpha",
    grepl("^X", lineage) ~ "Recombinant",
    TRUE ~ "Other"
  )
}

lineage_gisaid_clean <- df_gisaid_raw %>%
  mutate(lineage = as.character(Nextclade_pango)) %>%
  separate(seqName, into = c("virus_name", "accession_id", "Collection_date"), sep = "\\|") %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    month = floor_date(Collection_date, "month"),
    lineage_group = assign_lineage_group_unified(as.character(lineage))
  ) %>%
  select(accession_id, Collection_date, month, lineage, lineage_group) %>%
  distinct()

area_gisaid <- lineage_gisaid_clean %>%
  filter(between(Collection_date, as.Date("2022-07-01"), as.Date("2024-05-01")),
         !is.na(month), lineage_group != "Unassigned") %>%
  count(month, lineage_group) %>%
  complete(month, lineage_group = lineage_order, fill = list(n = 0)) %>%
  group_by(month) %>%
  mutate(
    proportion = n / sum(n),
    lineage_group = factor(lineage_group, levels = lineage_order)
  ) %>%
  ungroup()

p3 <- ggplot(area_gisaid, aes(x = month, y = proportion, fill = lineage_group, group = lineage_group)) +
  geom_area(position = "fill", color = "white", linewidth = 0.2, alpha = 0.92) +
  geom_vline(xintercept = period_lines, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = shared_colors, breaks = lineage_order, limits = lineage_order, name = NULL, drop = FALSE) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y", limits = c(as.Date("2022-07-01"), as.Date("2024-05-31")), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05)), breaks = c(0, .25, .50, .75, 1)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(color = "grey70", linewidth = 0.3),
    axis.text.x = element_text(hjust = 0.5),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
  ) +
  labs(
    title = "C. Percentage of SARS-CoV-2 Lineages (GISAID Cohort)",
    x = "Month and year of specimen collection",
    y = "Percentage"
  )

# ==============================================================================
# COMBINE WITH SINGLE COLLECTED LEGEND
# ==============================================================================
combined_figure <- (p1 / p2 / p3) + 
  plot_layout(heights = c(1, 0.7, 1), guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 8.5)
  )

print(combined_figure)

ggsave("Combined_SARS_CoV2_SingleLegend_Fixed.png", combined_figure, width = 10, height = 12, dpi = 300)
