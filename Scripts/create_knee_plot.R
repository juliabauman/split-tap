#!/usr/bin/env Rscript
# Generate knee plot to determine cell filtering threshold
# Snakemake inputs: combined_cDNA_counts.tsv
# Snakemake outputs: plot, summary

library(dplyr)
library(ggplot2)

# Read inputs from Snakemake
cdna_file <- snakemake@input[[1]]
plot_file <- snakemake@output[['plot']]
summary_file <- snakemake@output[['summary']]
test_thresh <- snakemake@params[['test_thresh']]

# Read cDNA data
cdna_data <- read.csv(cdna_file, sep='\t')

# Calculate total UMI counts per cell
cell_umi_counts <- cdna_data %>%
  group_by(cell) %>%
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count))

# Create knee plot
knee_plot <- ggplot(cell_umi_counts, aes(x = 1:nrow(cell_umi_counts), y = total_count)) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Knee Plot of Summed Counts",
    x = "Ranked Cells (log scale)",
    y = "Total UMIs (log scale)"
  ) +
  theme_minimal() +
  geom_vline(xintercept = test_thresh, linetype = "dashed", color = "red") +
  annotate("text", x = test_thresh, y = max(cell_umi_counts$total_count),
           label = paste0(test_thresh, " cells"), color = "red", vjust = -0.5)

ggsave(plot_file, plot = knee_plot, width = 8, height = 6, dpi = 300)

# Save summary
write.table(cell_umi_counts, summary_file, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
