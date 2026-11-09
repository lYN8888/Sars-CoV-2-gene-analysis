

library(tidyverse)   # Data manipulation and plotting
library(lubridate)   # Working with dates
library(scales)      # Percentage labels on graphs

library(readr)
cnx <- read_csv("data/Gisaid_data/Data.csv/gisaid_hcov-19_CMI lineage report.csv")
cei <- read_csv("data/Gisaid_data/Data.csv/gisaid_hcov-19_CRI lineage report.csv")
lp <- read_csv("data/Gisaid_data/Data.csv/gisaid_hcov-19_Lampang lineage report.csv")
mhs <- read_csv("data/Gisaid_data/Data.csv/gisaid_hcov-19_MaeHS lineage report.csv")
lph <- read_csv("data/Gisaid_data/Data.csv/gisaid_hcov-19_LAMPHUN lineage report.csv")
tak <- read_csv("data/Gisaid_data/Data.csv/gisaid_hcov-19_Tak lineage report.csv")
NCBI_virus <- read_csv("data/NCBI_virus_data/(NCBI_virus) seq.csv")

lineage_NCBI <- NCBI_virus %>%
  mutate(
    lineage = as.character(Pangolin))

all_variants <- bind_rows(
  cnx,
  cei,
  lp,
  mhs,
  lph,
  tak
)

all_variants <- all_variants %>%
  separate(
    taxon,
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

NCBI_data <- lineage_NCBI %>%
  mutate(
    Collection_Date = as.Date(Collection_Date, format = "%m/%d/%Y")
  )

NCBI_variants <- NCBI_data %>%
  mutate(
    accession_id = as.character(Accession))
NCBI_variants <- NCBI_data %>%
  transmute(
    accession_id = as.character(Accession),
    Collection_date = Collection_Date,
    month = floor_date(Collection_Date, "month"),
    lineage
  )


Meta_data <- bind_rows(call_variants, NCBI_variants) %>%
  select(
    accession_id,
    Collection_date,
    month,
    lineage) 

# Create a function to create group of lineages
assign_lineage_group <- function(lineage) {
  case_when(
    # BA.2 like Omicron and BA.1/BA.5 like Omicron
    grepl("^BA\\.2|^BA\\.2\\.10\\.1|^BA\\.2\\.1|^BA\\.2\\.9|^BA\\.2\\.9\\.5|^BN\\.1\\.1\\.1|^BN\\.1\\.2|^BN\\.1\\.3|^BN\\.1\\.2\\.3|^BN\\.1\\.3\\.6|^BS\\.1\\.1|^BA\\.2\\.75\\.2|^BL\\.1|^EJ\\.2|^BM\\.4\\.1\\.1|^CH\\.1\\.1|^BA\\.1|^BA\\.1\\.1|^BA\\.5\\.2", lineage) ~ "BA.2 like",
    # B.1.1.7-like Alpha 
    grepl("^B\\.1\\.1\\.7", lineage) ~ "Alpha",
    # B.1.617.2 like Delta
    grepl("^AY\\.85|^AY\\.79|^AY\\.30|^AY\\.118|^AY\\.31|^AY\\.92|^B\\.1\\.617|^B\\.1\\.617\\.2", lineage)  ~ "Delta",
    # XBB and all Pangolin aliases
    grepl("^XBB\\.1\\.16|^XBB\\.2\\.3\\.4|^XBB\\.1\\.16\\.1",lineage)                  ~ "XBB",
    # JN.1+ / BA.2.86 branch (before BA.2 catch)
    grepl("^JN\\.|^JF\\.|^KB\\.|^KP\\.|^LB\\.|^MC\\.|^NB\\.|^KS\\.|^MV\\.", lineage) ~ "JN.1+",
    TRUE                                                             ~ "Other"
  )
}


# Apply the function to our data 

lineage_data <- Meta_data %>%
  mutate( lineage_group = assign_lineage_group(as.character(lineage))) %>%
  select(
    accession_id,
    Collection_date,
    month,
    lineage,
    lineage_group) %>%
  
  distinct()



## Assing colors for each lineage
lineage_colors <- c(
  "BA.2 like" = "#6BAED6",
  "Alpha"  = "#587F24",
  "Delta" = "#2171B5",
  "XBB"        = "#41AB5D",
  "JN.1+"    = "#E87F24",
  "Other"      = "#BAB0AC"
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

p_area <- ggplot(lineage_area, aes(x = month, y = proportion,
                                   fill = lineage_group, group = lineage_group)) +
  geom_area(position = "fill", color = "white", linewidth = 0.2, alpha = 0.92) +
  scale_fill_manual(values = lineage_colors, name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 3)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y", expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0),
                     breaks = c(0, .25, .50, .75, 1),
                     sec.axis = dup_axis(name = NULL)) +
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
    plot.subtitle   = element_text(color = "grey40", size = 10)
  ) +
  labs(
    title    = "Percentage of hcov-19 Lineages",
    subtitle = "GISAID+NCBI_virus Cohort, 2022–2024",
    x = "Month and year of specimen collection",
    y = "Percentage of total lineages"
  )
p_area
  


