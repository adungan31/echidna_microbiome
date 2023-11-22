#set working directory
setwd("C:/Users/dungana/Dropbox/ARC_DP_MinimalMicrobiome/Projects/BacteriaViability/Data/Metabarcoding/R/input_files")

#packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(data.table)


#Import data
my_data <- read.csv("observed_features.csv", header = TRUE, row.names = "sample.id")

#move non_numeric data to a separate data table and remove from my_data
#For my dataset, the last 17 columns were non_numeric
non_numeric <- my_data[, (ncol(my_data) - 16):ncol(my_data)]
my_data <- my_data[, -((ncol(my_data)-16):ncol(my_data))]

# Create a new data frame to store the results with the correct data types
result <- data.frame(group_of_10 = character(), mean_value = numeric(), SE = numeric(), stringsAsFactors = FALSE)

# Calculate the mean and SE for each group of 10 columns using lapply
group_means <- lapply(seq(from = 1, to = ncol(my_data), by = 10), function(i) {
  group_name <- gsub(".*[.]([0-9]+)[_]?iter[.][0-9]+", "\\1", names(my_data)[i])
  group_mean <- apply(my_data[, i:(i+9)], 1, mean, na.rm = TRUE)
  group_SE <- apply(my_data[, i:(i+9)], 1, function(x) sd(x, na.rm = TRUE) / sqrt(10))
  data.frame(group_of_10 = group_name, mean_value = group_mean, SE = group_SE)
})

# Combine the results into a data frame
result <- do.call(rbind, group_means)


# Add current row names as a new column
result <- rownames_to_column(result, var = "sample")
#keep only the first 12 characters
result$sample <- substr(result$sample, 1, 12)

# Create a new column called "Species" in result
result$Species <- NA

# Loop through the row names of the non_numeric table and match with the sample column in the result table
for (i in rownames(non_numeric)) {
  sample_name <- substr(i, 1, 12)
  rows_to_update <- which(result$sample == sample_name)
  result[rows_to_update, "Species"] <- non_numeric[i, "Species"]
}


result$group_of_10 <- as.numeric(result$group_of_10)
data <- na.omit(result)
data <- data.frame(data)
data$SE <- NULL

# Calculate mean and SE for each Species and group_of_10

result_summary <- data %>%
  filter(!Species %in% c("DNA Extract Blank", "PCR Blank", "Blank")) %>%
  group_by(Species, group_of_10) %>%
  summarize(across(mean_value, list(mean_se)), .groups = "drop")

result_summary <- result_summary %>%
  mutate(SE = (mean_value_1$ymax - mean_value_1$ymin) / 2)

result_summary <- result_summary %>%
  mutate(mean_value = mean_value_1 %>% pull(y),
         ymin = mean_value_1 %>% pull(ymin),
         ymax = mean_value_1 %>% pull(ymax)) %>%
  mutate(SE = (ymax - ymin) / 2)

ggplot(result_summary, aes(x = group_of_10, y = mean_value, color = Species)) +
  geom_errorbar(aes(ymin = mean_value - SE, ymax = mean_value + SE), width = 0.2) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")) +
  labs(x = "Read Depth", y = "Observed ASVs") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "lightgray"),
    axis.line = element_line(color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  scale_x_continuous(limits = c(0, 10000))

