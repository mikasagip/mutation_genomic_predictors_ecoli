library(ggplot2)
library(dplyr)
#import in relevant libraries

df <- read.csv("all_properties.csv", stringsAsFactors = FALSE)
#read in csv

y_var   <- "mutation_density"
y_label <- "Mutation Density"
#assigns mutation density as y variable, also the proper label for the axis

var_list <- list(
  list(x_var = "relative_distance_oscillating_cosine",   x_label = "Relative Distance from Mid-replichore (Curved)"),
  list(x_var = "protein_abundance_lb",                   x_label = "Protein Abundance (LB media)")
)
#makes list of x variable associated with their proper label for later cycling through

bar_colours <- c(
  "#1e6091", "#19537c"
)
#list of colours to cycle through

bar_outline <- "black"
bar_alpha   <- 0.85
n_bins       <- 20
label_digits <- 3
#setting up some setting for plotting later on

make_length_weighted_bins <- function(data, x_var, n_bins) {
  
  data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(avg_length), !is.na(.data[[y_var]])) %>%
    arrange(.data[[x_var]])
  #sorts by x variable after removing NAs
  
  total_length  <- sum(data$avg_length)
  target_length <- total_length / n_bins
  #each bin should contain this much total sequence
  
  target_length
  total_length
  #prints the total and target length
  
  bin_id      <- integer(nrow(data))
  #stores bin assignment for each row
  current_bin <- 1L
  cumulative  <- 0
  
  for (k in seq_len(nrow(data))) {
    #loops through genes
    bin_id[k]  <- current_bin
    #assigns gene to current bin
    cumulative <- cumulative + data$avg_length[k]
    #adds gene length to running total for current bin
    if (cumulative >= target_length && current_bin < n_bins) {
      #if bin is full, move to the next bin
      current_bin <- current_bin + 1L
      cumulative  <- 0
      #rests the bin info
    }
  }
  
  data$bin <- bin_id
  #attaches bin label 
  
  bin_edges <- data %>%
    group_by(bin) %>%
    summarise(x_min = min(.data[[x_var]]),
              x_max = max(.data[[x_var]]),
              .groups = "drop")
  #for each bin finds minimun and maximum values
  
  breaks_vals <- c(bin_edges$x_min, tail(bin_edges$x_max, 1))
  #creates bin edge labels for lefts of bins + right of final bin 
  
  bin_summary <- data %>%
    group_by(bin) %>%
    summarise(
      mean_y = mean(.data[[y_var]]),
      se_y   = sd(.data[[y_var]]) / sqrt(n()),
      .groups = "drop"
    )
  #computes summary statistics for each bin to get the mean and standard error
  
  list(summary = bin_summary, breaks = breaks_vals)
}
#returns the summary statistics and bin broundaries

all_summaries <- lapply(var_list, function(v) {
  make_length_weighted_bins(df, v$x_var, n_bins)$summary
})
#applies the binning function to all the different variables

global_ymax <- max(sapply(all_summaries, function(s) {
  max(s$mean_y + s$se_y, na.rm = TRUE)
}))
#finds the global maximums for both mean and se from all variables

global_ymax <- global_ymax * 1.05
#adds 5% headroom above the tallest bar + error bar

make_plot <- function(x_var, x_label, bar_fill, bin_summary, breaks) {
  
  break_labels <- round(breaks, label_digits)
  #rounds the break labels to the pre-determined significangt figures
  
  ggplot(bin_summary, aes(x = bin, y = mean_y)) +
    geom_bar(
      stat   = "identity",
      fill   = bar_fill,
      colour = bar_outline,
      alpha  = bar_alpha,
      width  = 1
    ) +
    geom_errorbar(
      aes(ymin = mean_y - se_y, ymax = mean_y + se_y),
      #applies the error bars
      width     = 0.3,
      colour    = bar_outline,
      linewidth = 0.5
    ) +
    scale_x_continuous(
      breaks = c(seq_len(n_bins) - 0.5, n_bins + 0.5),
      labels = break_labels,
      #applies the break labels
      #aligns the bars with the labels
      expand = expansion(add = 0.5)
    ) +
    scale_y_continuous(
      limits = c(0, global_ymax),
      #applies the global y limit so all main bar graphs have the same y axis
      expand = expansion(mult = c(0, 0))
      #expand = 0 because limits already include headroom
    ) +
    labs(
      x = x_label,
      y = y_label
      #applies the proper labels 
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(colour = "grey50", size = 7),
      axis.title   = element_text(face = "bold"),
      axis.line    = element_line(colour = "grey30"),
      axis.text.x  = element_text(angle = 45, hjust = 1)
      #applies font sizes etc.
    )
}
#makes function for plotting

output_dir <- "figure_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)
#if output folder doesnt exit, create it

for (i in seq_along(var_list)) {
  #loop through x variables
  v         <- var_list[[i]]
  fill      <- bar_colours[[i]]
  binned    <- all_summaries[[i]]
  #retrieves information
  breaks    <- make_length_weighted_bins(df, v$x_var, n_bins)$breaks
  #re-use pre-computed summary but still need the break positions for axis labels
  
  p <- make_plot(v$x_var, v$x_label, fill, binned, breaks)
  #makes the plot
  
  fname <- paste0(output_dir, "/mutation_density_by_", v$x_var, ".png")
  ggsave(fname, plot = p, width = 8, height = 5, dpi = 300, bg = "white")
  #saves the plot
  message("Saved: ", fname)
}