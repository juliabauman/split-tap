#!/usr/bin/env Rscript
# Filter cells based on knee threshold and create filtered datasets
# Snakemake inputs: guide_cdna, guides, cdna, umi_summary
# Snakemake outputs: filtered_guide_cdna, filtered_cells, umi_histogram, stats

library(dplyr)
library(ggplot2)

# Read inputs from Snakemake
guide_cdna_file <- snakemake@input[['guide_cdna']]
guides_file <- snakemake@input[['guides']]
cdna_file <- snakemake@input[['cdna']]
umi_summary_file <- snakemake@input[['umi_summary']]

filtered_guide_cdna_file <- snakemake@output[['filtered_guide_cdna']]
filtered_cells_file <- snakemake@output[['filtered_cells']]
umi_histogram_file <- snakemake@output[['umi_histogram']]
stats_file <- snakemake@output[['stats']]

cell_thresh <- snakemake@params[['cell_thresh']]
max_umi_count <- snakemake@params[['max_umi_count']]

# Read UMI summary
cell_umi_counts <- read.table(umi_summary_file, header=TRUE, sep='\t')

# Get top N cells by UMI count
cells_above_knee <- cell_umi_counts %>%
  slice_head(n = cell_thresh) %>%
  pull(cell)

# Read guide data to filter cells with guides
guides <- read.table(guides_file, header=TRUE, sep='\t')
cells_with_guides <- unique(guides$cell)

# Filter to cells above knee AND with guide UMIs
cells_above_knee_wGuide <- cells_above_knee[cells_above_knee %in% cells_with_guides]

# Read guide_cDNA combined file
guide_cdna <- read.table(guide_cdna_file, header=TRUE, sep='\t')

# Calculate total UMIs per cell for histogram
kept_input <- guide_cdna %>%
  filter(cell %in% cells_above_knee_wGuide) %>%
  group_by(cell) %>%
  summarise(total_umis = sum(count, na.rm = TRUE))

# Filter out high UMI cells (potential multiplets)
highUMI_cells <- kept_input %>%
  filter(total_umis > max_umi_count) %>%
  pull(cell)

cells_final <- cells_above_knee_wGuide[!(cells_above_knee_wGuide %in% highUMI_cells)]

# Create UMI histogram
umi_p_cell_kept <- ggplot(kept_input, aes(x = total_umis)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  xlim(0, 300) +
  labs(
    title = paste0("Histogram of Total UMIs (Top ", cell_thresh, " Cells)"),
    x = "Total UMIs",
    y = "Number of Cells"
  ) +
  theme_minimal()

ggsave(umi_histogram_file, plot = umi_p_cell_kept, width = 8, height = 6, dpi = 300)

# Filter the guide_cDNA file
filtered_df <- guide_cdna %>%
  filter(cell %in% cells_final)

write.table(filtered_df, filtered_guide_cdna_file, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

# Save filtered cell list
writeLines(cells_final, filtered_cells_file)

# Write statistics
cdna_data <- read.csv(cdna_file, sep='\t')
subs_cdna <- cdna_data %>% filter(cell %in% cells_final)
umi_per_cell <- subs_cdna %>%
  group_by(cell) %>%
  summarise(total_count = sum(count)) %>%
  pull(total_count)

guides_filt <- guides %>% filter(cell %in% cells_final)
crop_per_cell <- guides_filt %>%
  filter(count > 1) %>%
  group_by(cell) %>%
  summarise(total_count = sum(count)) %>%
  pull(total_count)

# Write stats
stats_text <- paste0(
  "Total cells: ", nrow(cell_umi_counts), "\n",
  "Cells above knee threshold (", cell_thresh, "): ", length(cells_above_knee), "\n",
  "Cells with guide UMIs: ", length(cells_with_guides), "\n",
  "Cells above knee with guides: ", length(cells_above_knee_wGuide), "\n",
  "High UMI cells filtered: ", length(highUMI_cells), "\n",
  "Final cells retained: ", length(cells_final), "\n",
  "Average cDNA UMIs/cell: ", sprintf("%.2f", mean(umi_per_cell)), "\n",
  "Average CROP UMIs/cell: ", sprintf("%.2f", mean(crop_per_cell)), "\n"
)

writeLines(stats_text, stats_file)
