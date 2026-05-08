library(ggplot2)
library(dplyr)
library(patchwork)
#import in libraries 

df <- read.csv("all_properties.csv", stringsAsFactors = FALSE)
#import in the data from csv file

y_var   <- "mutation_density"
y_label <- "Mutation Density"
#assign y variable and label

var_list <- list(
  list(x_var = "gc_content",                             x_label = "GC Content"),
  list(x_var = "cpg_excess",                             x_label = "CpG Excess"),
  list(x_var = "avg_length",                             x_label = "Gene Length"),
  list(x_var = "transcript",                             x_label = "Transcript Abundance"),
  list(x_var = "protein_abundance_pax",                  x_label = "Protein Abundance (PaxDB)"),
  list(x_var = "relative_distance_from_oric",            x_label = "Relative Distance From oriC"),
  list(x_var = "relative_distance_oscillating_triangle", x_label = "Relative Distance from Mid-replichore (Linear)")
)
#makes list of all x variables to loop through

bar_colours <- c(
  "#d9ed92", "#b5e48c", "#99d98c", "#76c893", "#52b69a",
  "#34a0a4", "#168aad" 
)
#makes list of colours for the different graphs

bar_outline  <- "black"
bar_alpha    <- 0.85
n_bins       <- 20
label_digits <- 3
#some settings for the graphs

make_length_weighted_bins <- function(data, x_var, n_bins) {
  
  data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(avg_length), !is.na(.data[[y_var]])) %>%
    arrange(.data[[x_var]])
  #clean and sort data
  
  total_length  <- sum(data$avg_length)
  target_length <- total_length / n_bins
  #define bins size based on total sequence length
  
  bin_id      <- integer(nrow(data))
  current_bin <- 1L
  cumulative  <- 0
  
  for (k in seq_len(nrow(data))) {
    bin_id[k]  <- current_bin
    cumulative <- cumulative + data$avg_length[k]
    if (cumulative >= target_length && current_bin < n_bins) {
      current_bin <- current_bin + 1L
      cumulative  <- 0
    }
  }
  #cycles through genes
  #acculmulated gene lengths in bin 
  #moves to next bin once its full
  
  data$bin <- bin_id
  #store bin ids
  
  bin_edges <- data %>%
    group_by(bin) %>%
    summarise(x_min = min(.data[[x_var]]),
              x_max = max(.data[[x_var]]),
              .groups = "drop")
  #extracts bin boundaries 
  
  breaks_vals <- c(bin_edges$x_min, tail(bin_edges$x_max, 1))
  #finds the bin edges 
  
  bin_summary <- data %>%
    group_by(bin) %>%
    summarise(
      mean_y = mean(.data[[y_var]]),
      se_y   = sd(.data[[y_var]]) / sqrt(n()),
      .groups = "drop"
    )
  #compute summary statistics for each bin 
  
  list(summary = bin_summary, breaks = breaks_vals)
  #outputs results of summary 
}

all_summaries <- lapply(var_list, function(v) {
  make_length_weighted_bins(df, v$x_var, n_bins)$summary
})
#computes the summaries for all bins for all the x variables 
#stores the summaries

global_ymax <- max(sapply(all_summaries, function(s) {
  max(s$mean_y + s$se_y, na.rm = TRUE)
}))
global_ymax <- global_ymax * 1.05
#computes the global y max value 
#then adds 5% padding
#this is to make panels comparable in composite later

make_plot <- function(x_var, x_label, bar_fill, bin_summary, breaks) {
  #builds one plot per x variable
  
  break_labels <- round(breaks, label_digits)
  #creates the break labels
  
  #plotting bar graph
  ggplot(bin_summary, aes(x = bin, y = mean_y)) +
    geom_bar(
      stat   = "identity",
      fill   = bar_fill,
      colour = bar_outline,
      alpha  = bar_alpha,
      width  = 1
      #base of bar graph
    ) +
    geom_errorbar(
      aes(ymin = mean_y - se_y, ymax = mean_y + se_y),
      width     = 0.3,
      colour    = bar_outline,
      linewidth = 0.5
      #adds error bars (+- 1 SE)
    ) +
    scale_x_continuous(
      breaks = c(seq_len(n_bins) - 0.5, n_bins + 0.5),
      labels = break_labels,
      expand = expansion(add = 0.5)
      #applies ticks at breaks 
      #adds the lebels 
    ) +
    scale_y_continuous(
      limits = c(0, global_ymax),
      expand = expansion(mult = c(0, 0))
      #makes y scale fit with the global y max
    ) +
    labs(x = x_label, y = y_label
         #applies the nice labels tpo the axes
    ) +
    theme_classic(base_size = 10) +
    theme(
      #applies theme
      axis.title  = element_text(face = "bold", size = 14),
      axis.text   = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.tag    = element_text(face = "bold", size = 13),
      axis.line   = element_line(colour = "grey30")
    )
}

plot_list <- vector("list", length(var_list))
#creates empty list at length of all x variables

for (i in seq_along(var_list)) {
  #loop throgh the x variables
  v      <- var_list[[i]]
  fill   <- bar_colours[[i]]
  binned <- all_summaries[[i]]
  breaks <- make_length_weighted_bins(df, v$x_var, n_bins)$breaks
  #extract all needed info
  
  plot_list[[i]] <- make_plot(v$x_var, v$x_label, fill, binned, breaks) +
    labs(tag = paste0(letters[i], ")"))
  #makes all the plots 
  #adds a label for each panel
}

composite <- wrap_plots(plot_list, ncol = 2, nrow = 4) +
  plot_layout(guides = "keep")
#arranges the composite

output_dir <- "figure_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)
#sets up the output directory if it doesnt already exist

#saves the composite
ggsave(
  paste0(output_dir, "/composite_mutation_density.png"),
  plot   = composite,
  width  = 10,
  height = 16,
  dpi    = 300,
  bg     = "white"
)
message("Saved composite: ", output_dir, "/composite_mutation_density.png")