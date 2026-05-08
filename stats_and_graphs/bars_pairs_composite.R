library(ggplot2)
library(dplyr)
library(patchwork)
#import in libtraries

df <- read.csv("all_properties.csv", stringsAsFactors = FALSE)
#import in data

pair_list <- list(
  list(
    x_var   = "transcript",
    y_var   = "protein_abundance_pax",
    x_label = "Transcript Abundance",
    y_label = "Protein Abundance (PaxDB)",
    outname = "transcript_vs_pax"
  ),
  list(
    x_var   = "avg_length",
    y_var   = "gc_content",
    x_label = "Gene Length",
    y_label = "GC Content",
    outname = "length_vs_gc"
  ),
  list(
    x_var   = "transcript",
    y_var   = "binary_pec",
    x_label = "Transcript Abundance",
    y_label = "Essentiality (PEC)",
    outname = "transcript_vs_pec"
  ),
  list(
    x_var   = "cpg_excess",
    y_var   = "protein_abundance_pax",
    x_label = "CpG Excess",
    y_label = "Protein Abundance (PaxDB)",
    outname = "cpg_vs_pax"
  ),
  list(
    x_var   = "avg_length",
    y_var   = "transcript",
    x_label = "Gene Length",
    y_label = "Transcript Abundance",
    outname = "length_vs_transcript"
  )
)
#create embedded lists for each variable pair with labels

n_bins       <- 20
label_digits <- 3
#settings for plotting

bar_colours <- c("#f7b267", "#f79d65", "#f4845f", "#f27059", "#d95f4a")
bar_outline  <- "black"
bar_alpha    <- 0.85
#colours for graphs

make_length_weighted_bins <- function(data, x_var, y_var, n_bins) {
  
  data <- data %>%
    filter(!is.na(.data[[x_var]]),
           !is.na(.data[[y_var]]),
           !is.na(avg_length)) %>%
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
    summarise(
      x_min = min(.data[[x_var]]),
      x_max = max(.data[[x_var]]),
      .groups = "drop"
    )
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

make_plot <- function(x_label, y_label, bar_fill, bin_summary, breaks, ymax) {
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
      limits = c(0, ymax),
      expand = expansion(mult = c(0, 0))
      #applies the y axis tpo fir the ymax of the current data
    ) +
    labs(x = x_label, y = y_label
         #applies the nice labels tpo the axes
    ) +
    theme_classic(base_size = 10) +      # <-- base size for composite panels
    theme(
      #applies theme
      axis.title  = element_text(face = "bold", size = 14),
      axis.text   = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.tag    = element_text(face = "bold", size = 13),
      axis.line   = element_line(colour = "grey30")
    )
}

plot_list <- vector("list", length(pair_list))
#creates empty list at length of all x variables

for (i in seq_along(pair_list)) {
  #loop throgh the x variables
  p_info <- pair_list[[i]]
  fill   <- bar_colours[[i]]
  #extract all needed info
  
  binned <- make_length_weighted_bins(df, p_info$x_var, p_info$y_var, n_bins)
  #calculates the length weighted binning
  ymax   <- max(binned$summary$mean_y + binned$summary$se_y, na.rm = TRUE) * 1.05
  #finds the y max value and adds 5% padding at the top
  
  plot_list[[i]] <- make_plot(
    x_label     = p_info$x_label,
    y_label     = p_info$y_label,
    bar_fill    = fill,
    bin_summary = binned$summary,
    breaks      = binned$breaks,
    ymax        = ymax
  ) +
    labs(tag = paste0(letters[i], ")"))
  #makes all the plots 
  #adds a label for each panel
}

composite <- wrap_plots(plot_list, ncol = 2, nrow = 3) +
  plot_layout(guides = "keep")
#arranges the composite

output_dir <- "figure_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)
#sets up the output directory if it doesn't already exist

#saves the composite
ggsave(
  paste0(output_dir, "/composite_pairs.png"),
  plot   = composite,
  width  = 10,
  height = 12,      # <-- 3 rows so shorter than the 4-row version
  dpi    = 300,
  bg     = "white"
)
message("Saved composite: ", output_dir, "/composite_pairs.png")