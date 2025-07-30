
####  3rd January 2025 ####
### Quality control of SNP raw and filtered data 


# Load necessary libraries
library(ggplot2)
library(gridExtra)


setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/QTL_mapping_exp/RADseq_analysis_2023/")



# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Assuming the data is loaded from a CSV file
data <- read.table("gstacks_distribs.csv", header = TRUE, sep = ",")

# Ensure 'primary_kept' and 'kept_frac' are numeric
data$primary_kept <- as.numeric(data$primary_kept)
data$kept_frac <- as.numeric(data$kept_frac)

# Adjust sample names to be ordered factors for better axis appearance
data$sample <- factor(data$sample, levels = data$sample)

# Panel A: Number of primary kept alignments per sample
p1 <- ggplot(data, aes(x = primary_kept, y = reorder(sample, primary_kept))) +
  geom_point() +
  geom_vline(xintercept = mean(data$primary_kept, na.rm = TRUE), linetype = "dashed", color = "red") +
  xlab("Kept alignments (M)") +
  ylab("") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

# Panel B: Boxplot of primary kept alignments
p2 <- ggplot(data, aes(y = primary_kept)) +
  geom_boxplot(fill = "#0099FF", color = "black") +  # Blue boxplot
  xlab("") +
  ylab("Kept alignments (M)") +
  theme_minimal() +
  theme(axis.text.y = element_blank())

# Panel C: Fraction of alignments kept per sample
p3 <- ggplot(data, aes(x = kept_frac, y = reorder(sample, kept_frac))) +
  geom_point() +
  geom_vline(xintercept = mean(data$kept_frac, na.rm = TRUE), linetype = "dashed", color = "red") +
  xlab("Fraction of Alignments Kept") +
  ylab("") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

# Panel D: Boxplot of fraction of alignments kept
p4 <- ggplot(data, aes(y = kept_frac)) +
  geom_boxplot(fill = "#FFCC00", color = "black") +  # Yellow boxplot
  xlab("") +
  ylab("Fraction of Alignments Kept") +
  theme_minimal() +
  theme(axis.text.y = element_blank())

# Arrange plots in a grid, adjusting spacing and layout
grid.arrange(p1, p2, p3, p4, ncol = 2, heights = c(2, 1), widths = c(2, 1))


############################# 
# Load necessary libraries
library(ggplot2)
library(gridExtra)

# loaded from a CSV file
data <- read.table("gstacks_distribs.csv", header = TRUE, sep = ",")

# Ensure 'primary_kept' and 'kept_frac' are numeric
data$primary_kept <- as.numeric(data$primary_kept)
data$kept_frac <- as.numeric(data$kept_frac)

# Adjust sample names to be ordered factors for better axis appearance
data$sample <- factor(data$sample, levels = data$sample)

# Set y-axis limits for scatter plots (Kept Alignments)
y_limits_kept_alignments <- c(0, max(data$primary_kept, na.rm = TRUE))

# Panel A: Number of primary kept alignments per sample (with y limits)
p1 <- ggplot(data, aes(x = primary_kept, y = reorder(sample, primary_kept))) +
  geom_point() +
  geom_vline(xintercept = mean(data$primary_kept, na.rm = TRUE), linetype = "dashed", color = "red") +
  xlab("Kept alignments (M) (log10)") +
  ylab("") +
  theme_minimal() +
  scale_x_log10() +
  #scale_x_continuous(limits = c(0, max(data$primary_kept))) +  # Set x limits for better visualization
  theme(axis.text.y = element_text(size = 6))

# Panel B: Boxplot of primary kept alignments with axes swapped
p2 <- ggplot(data, aes(y = primary_kept)) +
  geom_boxplot(fill = "#0099FF", color = "black") +  # Blue boxplot
  ylab("Kept alignments (M) (log10)") +
  xlab("") +
  theme_minimal() +
  coord_flip() +  # Swap x and y axis
  scale_x_log10()
  #scale_y_continuous(limits = c(0, max(data$primary_kept)))  # Set y-axis limits

# Panel C: Fraction of alignments kept per sample (with y limits)
p3 <- ggplot(data, aes(x = kept_frac, y = reorder(sample, kept_frac))) +
  geom_point() +
  geom_vline(xintercept = mean(data$kept_frac, na.rm = TRUE), linetype = "dashed", color = "red") +
  xlab("Fraction of Alignments Kept") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +  # Set x limits for fraction (0 to 1)
  theme(axis.text.y = element_text(size = 6))

# Panel D: Boxplot of fraction of alignments kept with axes swapped
p4 <- ggplot(data, aes(y = kept_frac)) +
  geom_boxplot(fill = "#FFCC00", color = "black") +  # Yellow boxplot
  ylab("Fraction of Alignments Kept") +
  xlab("") +
  theme_minimal() +
  coord_flip() +  # Swap x and y axis
  scale_y_continuous(limits = c(0, 1))  # Set y limits from 0 to 1 for fraction

# Arrange plots in a grid, adjusting spacing and layout
grid.arrange(p1, p2, p3, p4, ncol = 2, heights = c(2, 1), widths = c(2, 1))


# Arrange plots in a grid, adjusting spacing and layout
final_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, heights = c(2, 1), widths = c(2, 1))

# Saving plots to PDF, SVG, and PNG formats

# Save as PDF
ggsave("gstacks_distribs_plot.pdf", plot = final_plot, width = 12, height = 8)

# Save as SVG
ggsave("gstacks_distribs_plot.svg", plot = final_plot, width = 12, height = 8)

# Save as PNG
ggsave("gstacks_distribs_plot.png", plot = final_plot, width = 12, height = 8, dpi = 300)


#######################################

# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Assuming the data is loaded from a CSV file
data <- read.table("gstacks_distribs.csv", header = TRUE, sep = ",")

# Ensure 'primary_kept' and 'kept_frac' are numeric
data$primary_kept <- as.numeric(data$primary_kept)
data$kept_frac <- as.numeric(data$kept_frac)

# Adjust sample names to be ordered factors for better axis appearance
data$sample <- factor(data$sample, levels = data$sample)

# Define a custom theme to add borders to all plots
custom_theme_with_border <- theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add black borders to each plot
  )

# Panel A: Number of primary kept alignments per sample (log scale)
p1 <- ggplot(data, aes(x = primary_kept, y = reorder(sample, primary_kept))) +
  geom_point() +
  geom_vline(xintercept = mean(data$primary_kept, na.rm = TRUE), linetype = "dashed", color = "red") +
  xlab("Kept alignments (M, log scale)") +
  ylab("") +
  scale_x_log10() +  # Apply log scale to x-axis
  custom_theme_with_border +
  theme(axis.text.y = element_text(size = 6))

# Panel B: Boxplot of primary kept alignments with axes swapped (log scale)
p2 <- ggplot(data, aes(y = primary_kept)) +
  geom_boxplot(fill = "#0099FF", color = "black") +  # Blue boxplot
  ylab("Kept alignments (M, log scale)") +
  xlab("") +
  coord_flip() +  # Swap x and y axis
  scale_y_log10() +  # Apply log scale to y-axis
  custom_theme_with_border

# Panel C: Fraction of alignments kept per sample (normal scale)
p3 <- ggplot(data, aes(x = kept_frac, y = reorder(sample, kept_frac))) +
  geom_point() +
  geom_vline(xintercept = mean(data$kept_frac, na.rm = TRUE), linetype = "dashed", color = "red") +
  xlab("Fraction of Alignments Kept") +
  ylab("") +
  scale_x_continuous(limits = c(0, 1)) +  # Set x limits for fraction (0 to 1)
  custom_theme_with_border +
  theme(axis.text.y = element_text(size = 6))

# Panel D: Boxplot of fraction of alignments kept with axes swapped (normal scale)
p4 <- ggplot(data, aes(y = kept_frac)) +
  geom_boxplot(fill = "#FFCC00", color = "black") +  # Yellow boxplot
  ylab("Fraction of Alignments Kept") +
  xlab("") +
  coord_flip() +  # Swap x and y axis
  scale_y_continuous(limits = c(0, 1)) +  # Set y limits from 0 to 1 for fraction
  custom_theme_with_border

# Arrange plots in a grid, adjusting spacing and layout
grid.arrange(p1, p2, p3, p4, ncol = 2, heights = c(2, 1), widths = c(2, 1))

# Arrange plots in a grid, adjusting spacing and layout
final_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, heights = c(2, 1), widths = c(2, 1))

# Saving plots to PDF, SVG, and PNG formats

# Save as PDF
ggsave("gstacks_distribs_plot.pdf", plot = final_plot, width = 12, height = 8)

# Save as SVG
ggsave("gstacks_distribs_plot.svg", plot = final_plot, width = 12, height = 8)

# Save as PNG
ggsave("gstacks_distribs_plot.png", plot = final_plot, width = 12, height = 8, dpi = 300)




######################################## Analysis after SNPs calling ##############################
################## RAW VCF density plot ######################

# Install vcfR package if not already installed
install.packages("vcfR")

setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/QTL_mapping_exp/RADseq_analysis_2023/population/popul_vcf_final/vcf_processing/")

# Load necessary libraries
library(vcfR)
library(readr)
library(ggplot2)
library(tidyverse)
library(gplots)
library(tidyverse)
library(ggplot2)

#### raw data 

var_depth <- read_delim("raw_data.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() + xlim (0, 100) 
summary(var_depth$mean_depth)

var_miss <- read_delim("raw_data.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "lightgreen", colour = "black", alpha = 0.3)+theme_light() 
summary(var_miss$fmiss)


# Reading the frequency data (you can also suppress the column type message)
var_freq <- read_delim("raw_data.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf)
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

#var_freq <- read_delim("raw_data.frq", delim = "\t",
#                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

d <- ggplot(var_freq, aes(maf)) + geom_density(fill = "lightgreen", colour = "black", alpha = 0.3)+theme_light()  
summary(var_freq$maf)


ind_depth <- read_delim("raw_data.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)


e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() 
summary(ind_depth$depth) 


ind_miss  <- read_delim("raw_data.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() 


ind_het <- read_delim("raw_data.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
summary(ind_het$f)
g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() 

library(gridExtra)
grid.arrange(b, c, d, e, f, g)

# Arrange the plots into a grid
plot_grid <- arrangeGrob(b, c, d, e, f, g, ncol = 2)  # Customize ncol/nrow based on your desired layout

# Save the grid in different formats
ggsave("grid_raw_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_raw_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_raw_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)



##### Raw VCF with summmary stats on the plots #####


library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

#### raw data 

# Load the mean depth data
var_depth <- read_delim("raw_data.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# Summarize mean depth data
depth_summary <- summary(var_depth$mean_depth)
depth_annotation <- paste0(
  "Min: ", round(depth_summary[1], 2), 
  "\n1st Qu.: ", round(depth_summary[2], 2), 
  "\nMedian: ", round(depth_summary[3], 2),
  "\nMean: ", round(depth_summary[4], 2),
  "\n3rd Qu.: ", round(depth_summary[5], 2),
  "\nMax: ", round(depth_summary[6], 2)
)

# Plot the mean depth
b <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  xlim(0, 100) + 
  annotate("text", x = 50, y = 0.05, label = depth_annotation, size = 3, hjust = 0) + xlab("variant_depth")

# Load missing data per site
var_miss <- read_delim("raw_data.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize missing data
miss_summary <- summary(var_miss$fmiss)
miss_annotation <- paste0(
  "Min: ", round(miss_summary[1], 4), 
  "\n1st Qu.: ", round(miss_summary[2], 4), 
  "\nMedian: ", round(miss_summary[3], 4),
  "\nMean: ", round(miss_summary[4], 4),
  "\n3rd Qu.: ", round(miss_summary[5], 4),
  "\nMax: ", round(miss_summary[6], 4)
)

# Plot missing data
c <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.3, y = 1.26, label = miss_annotation, size = 3, hjust = 0) + xlab("variant_missing")

# Load frequency data
var_freq <- read_delim("raw_data.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf)
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

# Summarize MAF
maf_summary <- summary(var_freq$maf)
maf_annotation <- paste0(
  "Min: ", round(maf_summary[1], 6), 
  "\n1st Qu.: ", round(maf_summary[2], 6), 
  "\nMedian: ", round(maf_summary[3], 6),
  "\nMean: ", round(maf_summary[4], 6),
  "\n3rd Qu.: ", round(maf_summary[5], 6),
  "\nMax: ", round(maf_summary[6], 6)
)

# Plot MAF
d <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.3, y = 3.5, label = maf_annotation, size = 3, hjust = 0)

# Load individual depth data
ind_depth <- read_delim("raw_data.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# Summarize individual depth
ind_depth_summary <- summary(ind_depth$depth)
ind_depth_annotation <- paste0(
  "Min: ", round(ind_depth_summary[1], 3), 
  "\n1st Qu.: ", round(ind_depth_summary[2], 3), 
  "\nMedian: ", round(ind_depth_summary[3], 3),
  "\nMean: ", round(ind_depth_summary[4], 3),
  "\n3rd Qu.: ", round(ind_depth_summary[5], 3),
  "\nMax: ", round(ind_depth_summary[6], 3)
)

# Plot individual depth
e <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 50, y = 35, label = ind_depth_annotation, size = 3, hjust = 0) + xlab("indiv_depth")

# Load individual missing data
ind_miss  <- read_delim("raw_data.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize individual missing data
ind_miss_summary <- summary(ind_miss$fmiss)
ind_miss_annotation <- paste0(
  "Min: ", round(ind_miss_summary[1], 4), 
  "\n1st Qu.: ", round(ind_miss_summary[2], 4), 
  "\nMedian: ", round(ind_miss_summary[3], 4),
  "\nMean: ", round(ind_miss_summary[4], 4),
  "\n3rd Qu.: ", round(ind_miss_summary[5], 4),
  "\nMax: ", round(ind_miss_summary[6], 4)
)

# Plot individual missing data
f <- ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.3, y = 40, label = ind_miss_annotation, size = 3, hjust = 0) + xlab("indiv_missing")

# Load individual heterozygosity
ind_het <- read_delim("raw_data.het", delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

# Summarize heterozygosity
het_summary <- summary(ind_het$f)
het_annotation <- paste0(
  "Min: ", round(het_summary[1], 5), 
  "\n1st Qu.: ", round(het_summary[2], 5), 
  "\nMedian: ", round(het_summary[3], 5),
  "\nMean: ", round(het_summary[4], 5),
  "\n3rd Qu.: ", round(het_summary[5], 5),
  "\nMax: ", round(het_summary[6], 5)
)

# Plot heterozygosity
g <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.45, y = 20, label = het_annotation, size = 3, hjust = 0) + xlab("het")

# Arrange all plots
grid.arrange(b, c, d, e, f, g)

plot_grid <- grid.arrange(b, c, d, e, f, g)


# Save the grid in different formats
ggsave("grid_raw_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_raw_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_raw_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)






#############################################################################
################## Filtered VCF density plot ######################

setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/QTL_mapping_exp/RADseq_analysis_2023/population/popul_vcf_final/vcf_processing/filt_data/")

# Load necessary libraries
library(vcfR)
library(readr)
library(ggplot2)
library(tidyverse)
library(gplots)
library(tidyverse)
library(ggplot2)

#### filtered data 

var_depth <- read_delim("filt_data.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() + xlim (0, 100) 
summary(var_depth$mean_depth)

var_miss <- read_delim("filt_data.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "lightgreen", colour = "black", alpha = 0.3)+theme_light() 
summary(var_miss$fmiss)


# Reading the frequency data (you can also suppress the column type message)
var_freq <- read_delim("filt_data.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf)
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

#var_freq <- read_delim("filt_data.frq", delim = "\t",
#                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

d <- ggplot(var_freq, aes(maf)) + geom_density(fill = "lightgreen", colour = "black", alpha = 0.3)+theme_light()  
summary(var_freq$maf)


ind_depth <- read_delim("filt_data.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)


e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() 
summary(ind_depth$depth) 


ind_miss  <- read_delim("filt_data.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() 


ind_het <- read_delim("filt_data.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
summary(ind_het$f)
g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + theme_light() 

library(gridExtra)
grid.arrange(b, c, d, e, f, g)

# Arrange the plots into a grid
plot_grid <- arrangeGrob(b, c, d, e, f, g, ncol = 2)  # Customize ncol/nrow based on your desired layout

# Save the grid in different formats
ggsave("grid_filt_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)


##### Filtered VCF with summmary stats on the plots #####


library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

#### filtered data 

# Load the mean depth data for filtered data
var_depth <- read_delim("filt_data.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# Summarize mean depth data
depth_summary <- summary(var_depth$mean_depth)
depth_annotation <- paste0(
  "Min: ", round(depth_summary[1], 2), 
  "\n1st Qu.: ", round(depth_summary[2], 2), 
  "\nMedian: ", round(depth_summary[3], 2),
  "\nMean: ", round(depth_summary[4], 2),
  "\n3rd Qu.: ", round(depth_summary[5], 2),
  "\nMax: ", round(depth_summary[6], 2)
)

# Plot the mean depth for filtered data
b <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  xlim(0, 100) + 
  annotate("text", x = 60, y = 0.02, label = depth_annotation, size = 3, hjust = 0)

# Load missing data per site for filtered data
var_miss <- read_delim("filt_data.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize missing data for filtered data
miss_summary <- summary(var_miss$fmiss)
miss_annotation <- paste0(
  "Min: ", round(miss_summary[1], 4), 
  "\n1st Qu.: ", round(miss_summary[2], 4), 
  "\nMedian: ", round(miss_summary[3], 4),
  "\nMean: ", round(miss_summary[4], 4),
  "\n3rd Qu.: ", round(miss_summary[5], 4),
  "\nMax: ", round(miss_summary[6], 4)
)

# Plot missing data for filtered data
c <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.23, y = 6.0, label = miss_annotation, size = 3, hjust = 0)

# Load frequency data for filtered data
var_freq <- read_delim("filt_data.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf) for filtered data
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

# Summarize MAF for filtered data
maf_summary <- summary(var_freq$maf)
maf_annotation <- paste0(
  "Min: ", round(maf_summary[1], 6), 
  "\n1st Qu.: ", round(maf_summary[2], 6), 
  "\nMedian: ", round(maf_summary[3], 6),
  "\nMean: ", round(maf_summary[4], 6),
  "\n3rd Qu.: ", round(maf_summary[5], 6),
  "\nMax: ", round(maf_summary[6], 6)
)

# Plot MAF for filtered data
d <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.15, y = 3.0, label = maf_annotation, size = 3, hjust = 0)

# Load individual depth data for filtered data
ind_depth <- read_delim("filt_data.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# Summarize individual depth for filtered data
ind_depth_summary <- summary(ind_depth$depth)
ind_depth_annotation <- paste0(
  "Min: ", round(ind_depth_summary[1], 3), 
  "\n1st Qu.: ", round(ind_depth_summary[2], 3), 
  "\nMedian: ", round(ind_depth_summary[3], 3),
  "\nMean: ", round(ind_depth_summary[4], 3),
  "\n3rd Qu.: ", round(ind_depth_summary[5], 3),
  "\nMax: ", round(ind_depth_summary[6], 3)
)

# Plot individual depth for filtered data
e <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 110, y = 26, label = ind_depth_annotation, size = 3, hjust = 0)

# Load individual missing data for filtered data
ind_miss  <- read_delim("filt_data.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize individual missing data for filtered data
ind_miss_summary <- summary(ind_miss$fmiss)
ind_miss_annotation <- paste0(
  "Min: ", round(ind_miss_summary[1], 4), 
  "\n1st Qu.: ", round(ind_miss_summary[2], 4), 
  "\nMedian: ", round(ind_miss_summary[3], 4),
  "\nMean: ", round(ind_miss_summary[4], 4),
  "\n3rd Qu.: ", round(ind_miss_summary[5], 4),
  "\nMax: ", round(ind_miss_summary[6], 4)
)

# Plot individual missing data for filtered data
f <- ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.5, y = 50, label = ind_miss_annotation, size = 3, hjust = 0)

# Load individual heterozygosity for filtered data
ind_het <- read_delim("filt_data.het", delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

# Summarize heterozygosity for filtered data
het_summary <- summary(ind_het$f)
het_annotation <- paste0(
  "Min: ", round(het_summary[1], 5), 
  "\n1st Qu.: ", round(het_summary[2], 5), 
  "\nMedian: ", round(het_summary[3], 5),
  "\nMean: ", round(het_summary[4], 5),
  "\n3rd Qu.: ", round(het_summary[5], 5),
  "\nMax: ", round(het_summary[6], 5)
)

# Plot heterozygosity for filtered data
g <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = -0.8, y = 13, label = het_annotation, size = 3, hjust = 0)

# Arrange all plots
grid.arrange(b, c, d, e, f, g)

# Arrange the plots into a grid
plot_grid <- arrangeGrob(b, c, d, e, f, g, ncol = 2)  # Customize ncol/nrow based on your desired layout

# Save the grid in different formats
ggsave("grid_filt_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)



########################################################################################################################################
######################## PCA ###########################################################################################################

# Load the libraries
library(data.table)
library(ggplot2)

# Load PCA data from the .eigenvec file (adjust the path as needed)
pca_data <- fread("pca_results.eigenvec", header = FALSE)

# Assign column names: First two columns for IDs, remaining for PCs
colnames(pca_data) <- c("FID", "IID", paste0("PC", 1:(ncol(pca_data) - 2)))

# Check the first few rows to ensure data is loaded correctly
head(pca_data)



# Plot the first two principal components (PC1 and PC2)
ggplot(pca_data, aes(x = PC1, y = PC2, label = IID)) +
  geom_point(size = 2, color = "blue") +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  geom_text(vjust = 1.5, size = 3)  # Optional: add labels for each point



# Load the .fam file to get population or group information
fam_data <- fread("pca_results.fam", header = FALSE)

# Merge the PCA data with FAM file
# Assuming that population info is in column 6 of the FAM file
pca_data <- merge(pca_data, fam_data[, .(V1, V2, V6)], by.x = c("FID", "IID"), by.y = c("V1", "V2"))
colnames(pca_data)[ncol(pca_data)] <- "Population"

# Convert Population column to a factor
pca_data$Population <- as.factor(pca_data$Population)

# Optional: If -9 represents missing or invalid values, replace them
# Replace -9 with NA or another category (e.g., "Unknown")
pca_data$Population[pca_data$Population == -9] <- NA  # Or another value, like "Unknown"


# Now create the plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA Plot Colored by Population", x = "PC1", y = "PC2") +
  scale_color_manual(values = c("red", "blue", "green", "purple"))  # Customize the colors as needed



############################################
####### without_missing_individuals_SNP_call #############


###### RAW SNPS ##########

setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/QTL_mapping_exp/RADseq_analysis_2023/population/highly_paired_SNP_call/pop_vcf/")

library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

#### raw data 

# Load the mean depth data
var_depth <- read_delim("raw_data.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# Summarize mean depth data
depth_summary <- summary(var_depth$mean_depth)
depth_annotation <- paste0(
  "Min: ", round(depth_summary[1], 2), 
  "\n1st Qu.: ", round(depth_summary[2], 2), 
  "\nMedian: ", round(depth_summary[3], 2),
  "\nMean: ", round(depth_summary[4], 2),
  "\n3rd Qu.: ", round(depth_summary[5], 2),
  "\nMax: ", round(depth_summary[6], 2)
)

# Plot the mean depth
b <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  xlim(0, 100) + 
  annotate("text", x = 50, y = 0.05, label = depth_annotation, size = 3, hjust = 0)

# Load missing data per site
var_miss <- read_delim("raw_data.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize missing data
miss_summary <- summary(var_miss$fmiss)
miss_annotation <- paste0(
  "Min: ", round(miss_summary[1], 4), 
  "\n1st Qu.: ", round(miss_summary[2], 4), 
  "\nMedian: ", round(miss_summary[3], 4),
  "\nMean: ", round(miss_summary[4], 4),
  "\n3rd Qu.: ", round(miss_summary[5], 4),
  "\nMax: ", round(miss_summary[6], 4)
)

# Plot missing data
c <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.3, y = 1.26, label = miss_annotation, size = 3, hjust = 0)

# Load frequency data
var_freq <- read_delim("raw_data.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf)
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

# Summarize MAF
maf_summary <- summary(var_freq$maf)
maf_annotation <- paste0(
  "Min: ", round(maf_summary[1], 6), 
  "\n1st Qu.: ", round(maf_summary[2], 6), 
  "\nMedian: ", round(maf_summary[3], 6),
  "\nMean: ", round(maf_summary[4], 6),
  "\n3rd Qu.: ", round(maf_summary[5], 6),
  "\nMax: ", round(maf_summary[6], 6)
)

# Plot MAF
d <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.3, y = 3.5, label = maf_annotation, size = 3, hjust = 0)

# Load individual depth data
ind_depth <- read_delim("raw_data.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# Summarize individual depth
ind_depth_summary <- summary(ind_depth$depth)
ind_depth_annotation <- paste0(
  "Min: ", round(ind_depth_summary[1], 3), 
  "\n1st Qu.: ", round(ind_depth_summary[2], 3), 
  "\nMedian: ", round(ind_depth_summary[3], 3),
  "\nMean: ", round(ind_depth_summary[4], 3),
  "\n3rd Qu.: ", round(ind_depth_summary[5], 3),
  "\nMax: ", round(ind_depth_summary[6], 3)
)

# Plot individual depth
e <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 50, y = 35, label = ind_depth_annotation, size = 3, hjust = 0)

# Load individual missing data
ind_miss  <- read_delim("raw_data.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize individual missing data
ind_miss_summary <- summary(ind_miss$fmiss)
ind_miss_annotation <- paste0(
  "Min: ", round(ind_miss_summary[1], 4), 
  "\n1st Qu.: ", round(ind_miss_summary[2], 4), 
  "\nMedian: ", round(ind_miss_summary[3], 4),
  "\nMean: ", round(ind_miss_summary[4], 4),
  "\n3rd Qu.: ", round(ind_miss_summary[5], 4),
  "\nMax: ", round(ind_miss_summary[6], 4)
)

# Plot individual missing data
f <- ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.3, y = 40, label = ind_miss_annotation, size = 3, hjust = 0)

# Load individual heterozygosity
ind_het <- read_delim("raw_data.het", delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

# Summarize heterozygosity
het_summary <- summary(ind_het$f)
het_annotation <- paste0(
  "Min: ", round(het_summary[1], 5), 
  "\n1st Qu.: ", round(het_summary[2], 5), 
  "\nMedian: ", round(het_summary[3], 5),
  "\nMean: ", round(het_summary[4], 5),
  "\n3rd Qu.: ", round(het_summary[5], 5),
  "\nMax: ", round(het_summary[6], 5)
)

# Plot heterozygosity
g <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.45, y = 20, label = het_annotation, size = 3, hjust = 0)

# Arrange all plots
grid.arrange(b, c, d, e, f, g)


# Save the grid in different formats
ggsave("grid_raw_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_raw_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_raw_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)



##########################
###### filt SNPS #########
setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/QTL_mapping_exp/RADseq_analysis_2023/population/popul_vcf_final/")

library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)


#### filtered data 

# Load the mean depth data for filtered data
var_depth <- read_delim("filt_based_max_miss05.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# Summarize mean depth data
depth_summary <- summary(var_depth$mean_depth)
depth_annotation <- paste0(
  "Min: ", round(depth_summary[1], 2), 
  "\n1st Qu.: ", round(depth_summary[2], 2), 
  "\nMedian: ", round(depth_summary[3], 2),
  "\nMean: ", round(depth_summary[4], 2),
  "\n3rd Qu.: ", round(depth_summary[5], 2),
  "\nMax: ", round(depth_summary[6], 2)
)

# Plot the mean depth for filtered data
b <- ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  xlim(0, 100) + 
  annotate("text", x = 60, y = 0.02, label = depth_annotation, size = 3, hjust = 0) + xlab("variant_depth")

# Load missing data per site for filtered data
var_miss <- read_delim("filt_based_max_miss05.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize missing data for filtered data
miss_summary <- summary(var_miss$fmiss)
miss_annotation <- paste0(
  "Min: ", round(miss_summary[1], 4), 
  "\n1st Qu.: ", round(miss_summary[2], 4), 
  "\nMedian: ", round(miss_summary[3], 4),
  "\nMean: ", round(miss_summary[4], 4),
  "\n3rd Qu.: ", round(miss_summary[5], 4),
  "\nMax: ", round(miss_summary[6], 4)
)

# Plot missing data for filtered data
c <- ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.23, y = 6.0, label = miss_annotation, size = 3, hjust = 0) + xlab("variant_missing")

# Load frequency data for filtered data
var_freq <- read_delim("filt_based_max_miss05.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf) for filtered data
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

# Summarize MAF for filtered data
maf_summary <- summary(var_freq$maf)
maf_annotation <- paste0(
  "Min: ", round(maf_summary[1], 6), 
  "\n1st Qu.: ", round(maf_summary[2], 6), 
  "\nMedian: ", round(maf_summary[3], 6),
  "\nMean: ", round(maf_summary[4], 6),
  "\n3rd Qu.: ", round(maf_summary[5], 6),
  "\nMax: ", round(maf_summary[6], 6)
)

# Plot MAF for filtered data
d <- ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.15, y = 3.0, label = maf_annotation, size = 3, hjust = 0)

# Load individual depth data for filtered data
ind_depth <- read_delim("filt_based_max_miss05.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# Summarize individual depth for filtered data
ind_depth_summary <- summary(ind_depth$depth)
ind_depth_annotation <- paste0(
  "Min: ", round(ind_depth_summary[1], 3), 
  "\n1st Qu.: ", round(ind_depth_summary[2], 3), 
  "\nMedian: ", round(ind_depth_summary[3], 3),
  "\nMean: ", round(ind_depth_summary[4], 3),
  "\n3rd Qu.: ", round(ind_depth_summary[5], 3),
  "\nMax: ", round(ind_depth_summary[6], 3)
)

# Plot individual depth for filtered data
e <- ggplot(ind_depth, aes(depth)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 110, y = 26, label = ind_depth_annotation, size = 3, hjust = 0) + xlab("indiv_depth")

# Load individual missing data for filtered data
ind_miss  <- read_delim("filt_based_max_miss05.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize individual missing data for filtered data
ind_miss_summary <- summary(ind_miss$fmiss)
ind_miss_annotation <- paste0(
  "Min: ", round(ind_miss_summary[1], 4), 
  "\n1st Qu.: ", round(ind_miss_summary[2], 4), 
  "\nMedian: ", round(ind_miss_summary[3], 4),
  "\nMean: ", round(ind_miss_summary[4], 4),
  "\n3rd Qu.: ", round(ind_miss_summary[5], 4),
  "\nMax: ", round(ind_miss_summary[6], 4)
)

# Plot individual missing data for filtered data
f <- ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = 0.5, y = 50, label = ind_miss_annotation, size = 3, hjust = 0) + xlab("indiv_missing")

# Load individual heterozygosity for filtered data
ind_het <- read_delim("filt_based_max_miss05.het", delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

# Summarize heterozygosity for filtered data
het_summary <- summary(ind_het$f)
het_annotation <- paste0(
  "Min: ", round(het_summary[1], 5), 
  "\n1st Qu.: ", round(het_summary[2], 5), 
  "\nMedian: ", round(het_summary[3], 5),
  "\nMean: ", round(het_summary[4], 5),
  "\n3rd Qu.: ", round(het_summary[5], 5),
  "\nMax: ", round(het_summary[6], 5)
)

# Plot heterozygosity for filtered data
g <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  annotate("text", x = -0.8, y = 13, label = het_annotation, size = 3, hjust = 0)

# Arrange all plots
grid.arrange(b, c, d, e, f, g)

# Arrange the plots into a grid
plot_grid <- arrangeGrob(b, c, d, e, f, g, ncol = 2)  # Customize ncol/nrow based on your desired layout

# Save the grid in different formats
ggsave("grid_filt_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)





############################################################################################################################################################
###### those SNPs obtained fter filtering and comparing with the parent sites #########

setwd("/Users/Shared/Files From d.localized/PhD-Uni-Koeln_2021-2024/PhD_work/seeds data/QTL_mapping_exp/RADseq_analysis_2023/population/popul_vcf_final/missing_removed/filt_data_2/filt_parents_sites/new/")


library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)


#### filtered data 

# Load the mean depth data for filtered data
var_depth <- read_delim("filtered_F3_SNPs_checked_parents_vcf_stats.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# Summarize mean depth data
depth_summary <- summary(var_depth$mean_depth)
depth_annotation <- paste0(
  "Min: ", round(depth_summary[1], 2), 
  "\n1st Qu.: ", round(depth_summary[2], 2), 
  "\nMedian: ", round(depth_summary[3], 2),
  "\nMean: ", round(depth_summary[4], 2),
  "\n3rd Qu.: ", round(depth_summary[5], 2),
  "\nMax: ", round(depth_summary[6], 2)
)

# Plot the mean depth for filtered data
b <- ggplot(var_depth, aes(mean_depth)) + geom_histogram(fill = "blue", colour = "blue", alpha = 0.3) + 
  #geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + 
  xlim(0, 100) +
  theme(
    # Increase axis text sizes
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    
    # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) + xlab("variant_mean_depth")
  #annotate("text", x = 60, y = 0.02, label = depth_annotation, size = 3, hjust = 0) + xlab("variant_depth")

ggsave("Variant_mean_depth.pdf", plot = b, device = "pdf", width = 14, height = 10)


# Load missing data per site for filtered data
var_miss <- read_delim("filtered_F3_SNPs_checked_parents_vcf_stats.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize missing data for filtered data
miss_summary <- summary(var_miss$fmiss)
miss_annotation <- paste0(
  "Min: ", round(miss_summary[1], 4), 
  "\n1st Qu.: ", round(miss_summary[2], 4), 
  "\nMedian: ", round(miss_summary[3], 4),
  "\nMean: ", round(miss_summary[4], 4),
  "\n3rd Qu.: ", round(miss_summary[5], 4),
  "\nMax: ", round(miss_summary[6], 4)
)

# Plot missing data for filtered data
c <- ggplot(var_miss, aes(fmiss)) + geom_histogram(fill = "blue", colour = "blue", alpha = 0.3) +
  #geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + theme(
    # Increase axis text sizes
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    
    # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )  + xlab("variant_missing")
 # annotate("text", x = 0.23, y = 6.0, label = miss_annotation, size = 3, hjust = 0) + xlab("variant_missing")

ggsave("Variant_missing.pdf", plot = c, device = "pdf", width = 14, height = 10)


# Load frequency data for filtered data
var_freq <- read_delim("filtered_F3_SNPs_checked_parents_vcf_stats.frq", delim = "\t", 
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq_ref", "freq_alt"), 
                       skip = 1, show_col_types = FALSE)

# Calculate Minor Allele Frequency (maf) for filtered data
var_freq$maf <- apply(var_freq[, c("freq_ref", "freq_alt")], 1, min)

# Summarize MAF for filtered data
maf_summary <- summary(var_freq$maf)
maf_annotation <- paste0(
  "Min: ", round(maf_summary[1], 6), 
  "\n1st Qu.: ", round(maf_summary[2], 6), 
  "\nMedian: ", round(maf_summary[3], 6),
  "\nMean: ", round(maf_summary[4], 6),
  "\n3rd Qu.: ", round(maf_summary[5], 6),
  "\nMax: ", round(maf_summary[6], 6)
)

# Plot MAF for filtered data
d <- ggplot(var_freq, aes(maf)) + geom_histogram(fill = "blue", colour = "blue", alpha = 0.3) +
  #geom_density(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + theme(
    # Increase axis text sizes
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    
    # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) + ylab("Site")
  #annotate("text", x = 0.15, y = 3.0, label = maf_annotation, size = 3, hjust = 0)

ggsave("Minor_allele_frequency.pdf", plot = d, device = "pdf", width = 14, height = 10)



# Load individual depth data for filtered data
ind_depth <- read_delim("filtered_F3_SNPs_checked_parents_vcf_stats.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# Summarize individual depth for filtered data
ind_depth_summary <- summary(ind_depth$depth)
ind_depth_annotation <- paste0(
  "Min: ", round(ind_depth_summary[1], 3), 
  "\n1st Qu.: ", round(ind_depth_summary[2], 3), 
  "\nMedian: ", round(ind_depth_summary[3], 3),
  "\nMean: ", round(ind_depth_summary[4], 3),
  "\n3rd Qu.: ", round(ind_depth_summary[5], 3),
  "\nMax: ", round(ind_depth_summary[6], 3)
)

# Plot individual depth for filtered data
e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "blue", colour = "blue", alpha = 0.3) +
  #geom_histogram(fill = "lightgreen", colour = "black", alpha = 0.3) + 
  theme_light() + theme(
    # Increase axis text sizes
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    
    # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) + xlab("indiv_depth")
  #annotate("text", x = 110, y = 26, label = ind_depth_annotation, size = 3, hjust = 0) + xlab("indiv_depth")

ggsave("Individual_depth.pdf", plot = e, device = "pdf", width = 14, height = 10)


# Load individual missing data for filtered data
ind_miss  <- read_delim("filtered_F3_SNPs_checked_parents_vcf_stats.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize individual missing data for filtered data
ind_miss_summary <- summary(ind_miss$fmiss)
ind_miss_annotation <- paste0(
  "Min: ", round(ind_miss_summary[1], 4), 
  "\n1st Qu.: ", round(ind_miss_summary[2], 4), 
  "\nMedian: ", round(ind_miss_summary[3], 4),
  "\nMean: ", round(ind_miss_summary[4], 4),
  "\n3rd Qu.: ", round(ind_miss_summary[5], 4),
  "\nMax: ", round(ind_miss_summary[6], 4)
)

# Plot individual missing data for filtered data
f <- ggplot(ind_miss, aes(fmiss)) + 
  geom_histogram(fill = "blue", colour = "blue", alpha = 0.3) + 
  theme_light() + theme(
    # Increase axis text sizes
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    
    # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) + xlab("indiv_missing")
  #annotate("text", x = 0.5, y = 50, label = ind_miss_annotation, size = 3, hjust = 0) + xlab("indiv_missing")

ggsave("Individual_missing.pdf", plot = f, device = "pdf", width = 14, height = 10)


# Load individual heterozygosity for filtered data
ind_het <- read_delim("filtered_F3_SNPs_checked_parents_vcf_stats.het", delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1)

# Summarize heterozygosity for filtered data
het_summary <- summary(ind_het$f)
het_annotation <- paste0(
  "Min: ", round(het_summary[1], 5), 
  "\n1st Qu.: ", round(het_summary[2], 5), 
  "\nMedian: ", round(het_summary[3], 5),
  "\nMean: ", round(het_summary[4], 5),
  "\n3rd Qu.: ", round(het_summary[5], 5),
  "\nMax: ", round(het_summary[6], 5)
)

# Plot heterozygosity for filtered data
g <- ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "blue", colour = "blue", alpha = 0.3) + 
  theme_light() + theme(
    # Increase axis text sizes
    axis.text.x = element_text(size = 14),   # X-axis text size
    axis.text.y = element_text(size = 14),   # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    
    # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) + xlab("Heterozgousity(f)")
  #annotate("text", x = -0.8, y = 13, label = het_annotation, size = 3, hjust = 0)

ggsave("Heterozygousity_f.pdf", plot = g, device = "pdf", width = 14, height = 10)



# Arrange all plots
grid.arrange(b, c, d, e, f, g)

# Arrange the plots into a grid
plot_grid <- arrangeGrob(b, c, d, e, f, g, ncol = 2)  # Customize ncol/nrow based on your desired layout

# Save the grid in different formats
ggsave("grid_filt_vcf_stats.pdf", plot = plot_grid, device = "pdf", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.svg", plot = plot_grid, device = "svg", width = 14, height = 10)
ggsave("grid_filt_vcf_stats.png", plot = plot_grid, device = "png", width = 14, height = 10, dpi = 300)




