library(tidyverse)   # Data manipulation and plotting
library(lubridate)   # Working with dates
library(scales)      # Percentage labels on graphs

library(readr)
df_gisaid_raw <- readr::read_delim(
  "C:/Users/Lenovo/OneDrive/Documents/Covid_study/data/Gisaid_data/Data.csv/metaGISAID22_24_nxtclde_nohighcover.csv",
  delim = ";",
  show_col_types = FALSE
)
lineage_GISAID <- df_gisaid_raw %>%
  mutate(
    lineage = as.character(Nextclade_pango))

all_variants <- lineage_GISAID %>%
  separate(
    seqName,
    into = c("virus_name",
             "accession_id",
             "Collection_date"),
    sep = "\\|"
  )

call_variants <- all_variants %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    month = floor_date(Collection_date, "month")
  )

GISAID_data <- call_variants %>%
  mutate(
    Collection_Date = as.Date(Collection_date, format = "%Y/%m/%d")
  )


GISAIDLINEAGE_data <- GISAID_data %>%
  select(
    accession_id,
    Collection_date,
    month,
    lineage) 

# Create a function to create group of lineages

library(stringr)
library(dplyr)

assign_lineage_group <- function(lineage) {
  case_when(
    # --- ERA 3 (Oct 2023 – May 2024): JN.1 Era ---
    # Fixed: added generalized patterns for JN sub-branches like JN.11, JN.14, etc.
    grepl("^JN\\.|^JN1|^JN\\.1|^JN11|^JN14|^KP\\.|^KB\\.|^LB\\.|^JF\\.|^JS\\.", lineage) ~ "JN.1+ (Oct 2023-April 2024)",
    
    # --- ERA 2 (Mar 2023 – Oct 2023): XBB Era ---
    # Added HN sublineages here as well
    grepl("^XBB|^XBL\\.|^XBU\\.|^EG\\.|^GJ\\.|^GE\\.|^CH\\.|^BM\\.|^BQ\\.|^BL\\.|^BF\\.|^FL\\.|^FD\\.|^HK\\.|^LF\\.|^HH\\.|^FE\\.|^EU\\.|^GM\\.|^BR\\.|^BY\\.|^GR\\.|^HN\\.", lineage) ~ "XBB (Mar 2023-Oct 2023)",
    
    # --- ERA 1 (July 2022 – Mar 2023): BA.4/BA.5 & BA.2 Era ---
    grepl("^BA\\.4|^BA\\.5|^BE\\.|^BF\\.|^LZ|^FP\\.", lineage) ~ "BA.4/BA.5 (Jul 2022-Mar 2023)",
    grepl("^B\\.1\\.1\\.529|^BA\\.1|^BA\\.2|^EJ\\.|^BN\\.|^BS\\.", lineage) ~ "BA.2.75+/ BA.2.9+/ BA.2.10+",
    
    # --- Older Variants ---
    grepl("^AY\\.|^GY\\.|^B\\.1\\.617|^B\\.1\\.617\\.2", lineage) ~ "Delta",
    grepl("^B\\.1\\.1\\.7", lineage) ~ "Alpha",
    
    # --- Recombinants ---
    grepl("^X", lineage) ~ "Other Recombinant (X*)",
    
    TRUE ~ "Other"
  )
}
# Apply the function to our data 

lineage_data <- GISAIDLINEAGE_data %>%
  mutate( lineage_group = assign_lineage_group(as.character(lineage))) %>%
  select(
    accession_id,
    Collection_date,
    month,
    lineage,
    lineage_group) %>%
  
  distinct()



## Assign colors for each lineage group
lineage_colors <- c(
  "BA.2.75+/ BA.2.9+/ BA.2.10+"    = "#6BAED6",
  "Alpha"                          = "#587F24",
  "Delta"                          = "#2171B5",
  "BA.4/BA.5 (Jul 2022-Mar 2023)"  = "#2150AC",
  "XBB (Mar 2023-Oct 2023)"        = "#41AB5D",
  "JN.1+ (Oct 2023-April 2024)"      = "#E87F24",
  "Other Recombinant (X*)"         = "#9E9AC8",
  "Other"                          = "#BDBDBD"
)

lineage_area <- lineage_data %>%
  filter((between(Collection_date,
                  as.Date("2022-07-01"),
                  as.Date("2024-05-01"))),
         !is.na(month), lineage_group != "Unassigned") %>%
  count(month, lineage_group) %>%
  complete(month, lineage_group, fill = list(n = 0)) %>%
  group_by(month) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  group_by(lineage_group) %>%
  mutate(first_seen = min(month[n > 0])) %>%
  ungroup() %>%
  mutate(lineage_group = fct_reorder(lineage_group, first_seen))
nrow(lineage_data)

period_lines <- as.Date(c("2023-03-01", "2023-10-01"))

p_area <- ggplot(lineage_area, aes(x = month, y = proportion,
                                   fill = lineage_group, group = lineage_group)) +
  geom_area(position = "fill", color = "white", linewidth = 0.2, alpha = 0.92) +
  
  # --- Add Vertical Lines to Separate Periods ---
  geom_vline(xintercept = period_lines, linetype = "dashed", color = "black", linewidth = 0.6) +
  
  # --- Add Period Labels at the Top ---
  annotate("text", x = as.Date("2022-11-01"), y = 1.03, label = "Period 1:\nJuly 2022 - Mar 2023", 
           size = 3, fontface = "bold", vjust = 0) +
  annotate("text", x = as.Date("2023-06-15"), y = 1.03, label = "Period 2:\nMar 2023 - Oct 2023", 
           size = 3, fontface = "bold", vjust = 0) +
  annotate("text", x = as.Date("2024-01-15"), y = 1.03, label = "Period 3:\nOct 2023 - April 2024", 
           size = 3, fontface = "bold", vjust = 0) +
  
  scale_fill_manual(values = lineage_colors, name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 3)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y", expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12)),
                     breaks = c(0, .25, .50, .75, 1),
                     sec.axis = dup_axis(name = NULL)) +
  coord_cartesian(clip = "off") + # Allows period text labels to sit neatly above the plot area
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
    title    = "Percentage of hcov-19 Lineages",
    subtitle = "GISAID Cohort, 2022–2024",
    x = "Month and year of specimen collection",
    y = "Percentage of total lineages"
  )

print(p_area)
