library(ggplot2)
library(dplyr)
#import in libraries

df <- read.csv("all_properties.csv", stringsAsFactors = FALSE)
#import in data

y_var   <- "mutation_density"
y_label <- "Mutation Density"
#defines y variable and label

x_var = "middle"
x_label = "Middle Coordinate"
#defines x variable and label

bar_colour  <- "#1a759f"
smooth_colour <- "darkred"
#colours for bars and smooth line

bar_outline  <- "black"
bar_alpha    <- 0.85
n_bins       <- 20
label_digits <- 3
#settings for plotting

make_length_weighted_bins <- function(data, x_var, n_bins, origin = NULL) {
  data <- data %>%
    filter(!is.na(.data[[x_var]]), !is.na(avg_length), !is.na(.data[[y_var]]))
  #make sure data is clean
  
  genome_size <- max(data[[x_var]], na.rm = TRUE)
  #calculate total sequence amount 
  
  if (!is.null(origin)) {
    data[[x_var]] <- (data[[x_var]] - origin + genome_size / 2) %% genome_size
  }
  #rotates the data around oric so that oric will be centred in the graph
  
  data <- data %>% arrange(.data[[x_var]])
  #ensures bins follow genome order
  
  total_length  <- sum(data$avg_length)
  target_length <- total_length / n_bins
  #sets up target sequence amount for each bin
  
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
  #cycles through the bins
  #when bin is full, moves onto the next one
  
  data$bin <- bin_id
  #stores the bins ids
  
  bin_edges <- data %>%
    group_by(bin) %>%
    summarise(x_min = min(.data[[x_var]]),
              x_max = max(.data[[x_var]]),
              .groups = "drop")
  #gets the bin edges
  
  breaks_rotated <- c(bin_edges$x_min, tail(bin_edges$x_max, 1))
  breaks_original <- (breaks_rotated + origin - genome_size / 2) %% genome_size
  #converted the rotated labels back to the original coordinates
  
  bin_summary <- data %>%
    group_by(bin) %>%
    summarise(
      mean_y = mean(.data[[y_var]]),
      se_y   = sd(.data[[y_var]]) / sqrt(n()),
      .groups = "drop"
    )
  #computes summary stats 
  
  list(summary = bin_summary, breaks = breaks_original, data = data)
  #saves the summary stats
}

oric_pos <- 3925860
#defines position of oric

binned <- make_length_weighted_bins(df, x_var, n_bins, origin = oric_pos)
#uses the function to divide bins by cumulative sequence length

ymax <- max(binned$summary$mean_y + binned$summary$se_y, na.rm = TRUE)
ymax <- ymax * 1.05
#calculates the y max value and adds 5% padding 

make_plot <- function(x_var, x_label, bar_fill, bin_summary, breaks, raw_data) {
  #start building figure
  
  break_labels <- round(breaks, label_digits)
  #finds the break labels
  
  oric_x      <- (n_bins + 1) / 2
  ter_x_left  <- 0.5  
  ter_x_right <- n_bins + 0.5
  #finds positions for oric and ter labels 
  
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
      #applies  the error bars
    ) +
    geom_smooth(
      data    = bin_summary,
      aes(x = bin, y = mean_y),
      method  = "loess",
      span    = 0.3,
      se      = FALSE,
      colour  = smooth_colour,
      linewidth = 1,
      inherit.aes = FALSE
      #LOESS curve fitted on bin means
    ) +
    annotate("text", x = oric_x, y = ymax, label = "oriC",
             fontface = "bold", vjust = 1.2, size = 3.5
             #add oric label
    ) +
    annotate("segment", x = oric_x, xend = oric_x,
             y = ymax * 0.95, yend = 0,
             colour = "black", linewidth = 0.5, linetype = "dashed"
             #add line fpr oric label
    ) +
    annotate("text", x = ter_x_left, y = ymax, label = "ter",
             fontface = "bold", vjust = 1.2, size = 3.5
             #add ter label
    ) +
    annotate("segment", x = ter_x_left, xend = ter_x_left,
             y = ymax * 0.95, yend = 0,
             colour = "black", linewidth = 0.5, linetype = "dashed"
             #add line for ter label
    ) +
    annotate("text", x = ter_x_right, y = ymax, label = "ter",
             fontface = "bold", vjust = 1.2, size = 3.5
             #add ter label
    ) +
    annotate("segment", x = ter_x_right, xend = ter_x_right,
             y = ymax * 0.95, yend = 0,
             colour = "black", linewidth = 0.5, linetype = "dashed"
             #add line for ter label
    ) +
    scale_x_continuous(
      breaks = c(seq_len(n_bins) - 0.5, n_bins + 0.5),
      labels = break_labels,
      expand = expansion(add = 0.5)
      #sets up the x axis with the bin edge labels
    ) +
    scale_y_continuous(
      limits = c(0, ymax),
      expand = expansion(mult = c(0, 0)),
      oob = scales::squish
      #applies the ymax to the y axis
      #also ensures there are no clipping issues 
    ) +
    labs(x = x_label, y = y_label
         #add nice axis labels
    ) +
    theme_classic(base_size = 11) +
    theme(
      #apply theme
      plot.title   = element_text(face = "bold", hjust = 0.5),
      plot.caption = element_text(colour = "grey50", size = 7),
      axis.title   = element_text(face = "bold"),
      axis.line    = element_line(colour = "grey30"),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    )
}

output_dir <- "figure_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)
#create output doirectory if it doesnt already exist

p <- make_plot(
  x_var,
  x_label,
  bar_colour,
  binned$summary,
  binned$breaks,
  binned$data
)
#actually create the plot

fname <- paste0(output_dir, "/mutation_density_by_middle.png")
#saves name for file

#saves the figure
ggsave(
  fname,
  plot = p,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)