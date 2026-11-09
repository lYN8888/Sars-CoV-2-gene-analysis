
library(tidyverse)
library(scales)

library(readr)
lineage_data <- read_csv("~/Covid_study/data/lineage_data(samples 2022-2024).csv")

lineage_data <- lineage_data %>%
  mutate(
    collection_date = as.Date(collection_date,
                              format = "%m/%d/%Y"))

lineage_data <- lineage_data %>%
  mutate(
    month = as.Date(month,
    format = "%m/%d/%Y")
  )



## Assigning colors for each lineage
lineage_colors <- c(

  "BA.2.75+" = "#2171B5",
  "BA.4/BA.5"= "#084594",
  "XBB"      = "#41AB5D",
  "JN.1+"    = "#E87F24",
  "Other"    = "#BAB0AC"
)

lineage_area <- lineage_data %>%
  filter(!is.na(month), lineage_group != "Unassigned") %>%
  count(month, lineage_group) %>%
  complete(month, lineage_group, fill = list(n = 0)) %>%
  group_by(month) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  group_by(lineage_group) %>%
  mutate(first_seen = min(month[n > 0])) %>%
  ungroup() %>%
  mutate(lineage_group = fct_reorder(lineage_group, first_seen))
nrow(lineage_area)

period_lines <- as.Date(c("2023-03-01", "2023-10-01"))

p_area <- ggplot(lineage_area, aes(x = month, y = proportion,
                                   fill = lineage_group, group = lineage_group)) +
  geom_area(position = "fill", color = "white", linewidth = 0.2, alpha = 0.92) +
  
  # --- Add Vertical Lines to Separate Periods ---
  geom_vline(xintercept = period_lines, linetype = "dashed", color = "black", linewidth = 0.6) +
  
  # --- Add Period Annotations at the Top ---
  annotate("text", x = as.Date("2022-11-01"), y = 1.03, label = "Period 1:\nJuly 2022 - Mar 2023", 
           size = 3, fontface = "bold", vjust = 0) +
  annotate("text", x = as.Date("2023-06-15"), y = 1.03, label = "Period 2:\nMar 2023 - Oct 2023", 
           size = 3, fontface = "bold", vjust = 0) +
  annotate("text", x = as.Date("2024-01-15"), y = 1.03, label = "Period 3:\nOct 2023 - May 2024", 
           size = 3, fontface = "bold", vjust = 0) +
  
  scale_fill_manual(values = lineage_colors, name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 3)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y", expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12)),
                     breaks = c(0, .25, .50, .75, 1),
                     sec.axis = dup_axis(name = NULL)) +
  coord_cartesian(clip = "off") + # Allows text labels outside the plot boundary
  theme_minimal(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    legend.key.size = unit(0.45, "cm"),
    legend.text     = element_text(size = 9),
    axis.ticks      = element_line(color = "grey70", linewidth = 0.3),
    axis.text.x     = element_text(hjust = 0.5),
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(color = "grey40", size = 10),
    plot.margin     = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  labs(
    title    = "Percentage of SARS-CoV-2 Lineages",
    subtitle = "Chiang Mai Household Cohort, 2022–2024",
    x = "Month and year of specimen collection",
    y = "Percentage of total lineages"
  )

print(p_area)
