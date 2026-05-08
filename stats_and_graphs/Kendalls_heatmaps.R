library(ggplot2)
library(reshape2)
library(ggpubr)
#import in relevant libraries

file_list <- list(
  Model1 = list(
    tau          = "model1_kendall_tau.csv",
    pval         = "model1_kendall_pvalues.csv",
    partial_tau  = "model1_kendall_partial_tau.csv",
    partial_pval = "model1_kendall_partial_pvalues.csv"
  ),
  Model2 = list(
    tau          = "model2_kendall_tau.csv",
    pval         = "model2_kendall_pvalues.csv",
    partial_tau  = "model2_kendall_partial_tau.csv",
    partial_pval = "model2_kendall_partial_pvalues.csv"
  ),
  Model3 = list(
    tau          = "model3_kendall_tau.csv",
    pval         = "model3_kendall_pvalues.csv",
    partial_tau  = "model3_kendall_partial_tau.csv",
    partial_pval = "model3_kendall_partial_pvalues.csv"
  ),
  Model4 = list(
    tau          = "model4_kendall_tau.csv",
    pval         = "model4_kendall_pvalues.csv",
    partial_tau  = "model4_kendall_partial_tau.csv",
    partial_pval = "model4_kendall_partial_pvalues.csv"
  )
)
#nested list of file names for later importing

model_names <- names(file_list)
#exctracts model names for looping

output_dir <- "heatmap_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)
#assigns the output folder, if it doesnt exist, creates it

var_rename <- c(
  "avg_length"                              = "Length",
  "binary_asymmetry"                        = "Asymmetry",
  "relative_distance_from_oric"             = "Distance From oriC",
  "transcript"                              = "Transcript Abundance",
  "gc_content"                              = "GC Content",
  "cpg_excess"                              = "CpG Excess",
  "relative_distance_oscillating_triangle"  = "Distance from Mid-replichore (Linear)",
  "binary_pec"                              = "Essentiality (PEC)",
  "protein_abundance_pax"                   = "Protein Abundance (PaxDB)",
  "mutation_density"                        = "Mutation Density",
  "relative_distance_oscillating_cosine"    = "Distance from Mid-replichore (Curved)",
  "binary_gerdes"                           = "Essentiality (GERDES)",
  "protein_abundance_lb"                    = "Protein Abundance (LB media)"
)
#maps variable names to better labels 

ALPHA_ONE_STAR   <- 0.05     # *
ALPHA_TWO_STAR   <- 0.01     # **
ALPHA_THREE_STAR <- 0.005    # ***
#assigning significance thresholds

read_matrix <- function(filepath) {
  raw    <- read.csv(filepath, check.names = FALSE)
  rnames <- as.character(raw[[1]])
  #first coloumn is row names
  mat    <- as.matrix(raw[, -1, drop = FALSE])
  #removes the first column and maintains a numeric matrix
  rownames(mat) <- rnames
  #assigns the row names to the rows in the matrix
  if (!identical(colnames(mat), rnames)) colnames(mat) <- rnames
  #ensures the matrix is square and symmertical
  storage.mode(mat) <- "numeric"
  #ensures matrix is stored as numeric
  mat
}
#CSV to matrix function


lower_tri_na <- function(mat) {
  mat[upper.tri(mat, diag = TRUE)] <- NA
  mat
}
#function to only keep the lower triangle so there is only one value per comparison

melt_lower <- function(mat) {
  var_order  <- rownames(mat)
  df         <- melt(mat, varnames = c("row_var", "col_var"), value.name = "value", na.rm = TRUE)
  df$row_var <- factor(df$row_var, levels = var_order)
  df$col_var <- factor(df$col_var, levels = var_order)
  df
}
#converts matrix to dataframe

apply_labels <- function(x) {
  x <- as.character(x)
  renamed <- ifelse(x %in% names(var_rename), var_rename[x], gsub("_", " ", x))
  renamed
}
#assigns the pretty labels to the variable names
#if there isnt one it falls back to swapping _ with space

IS_ZERO_PVAL <- function(x) !is.na(x) & x <= 2e-308
#captures p-values that were originally 0

sig_stars <- function(p) {
  ifelse(IS_ZERO_PVAL(p), "***",
         ifelse(p < ALPHA_THREE_STAR, "***",
                ifelse(p < ALPHA_TWO_STAR, "**",
                       ifelse(p < ALPHA_ONE_STAR, "*", "ns"))))
}
#assigns significance stars to p values

format_pval <- function(p, max_chars = 7) {
  vapply(p, function(x) {
    #sets 7 as maximum number of characters 
    #applys to every p-value in the vector
    if (is.na(x)) return(NA_character_)
    if (IS_ZERO_PVAL(x)) return("<2e-308")
    #adds on the < to specific p-values
    plain <- formatC(x, format = "fg", digits = 3, flag = "#")
    #uses three significant figures 
    plain <- sub("\\.$", "", sub("0+$", "", plain))
    #removes trailling zeros and decimal points
    if (nchar(plain) <= max_chars) {
      return(plain)
    }
    #checks if it fits within the previously defined 7 character max
    formatC(x, format = "e", digits = 1)
    #falls back on scientific notation if p-value has too many characters
  }, character(1))
}
#formats p-value to be a decimal with limit of characters, if cant falls back on scintific notation

heatmap_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      #custom override for the rest of the theme
      plot.tag          = element_text(face = "bold", size = base_size + 2,
                                       family = "sans"),
      #panels labels 
      axis.text.x       = element_text(angle = 40, hjust = 1, vjust = 1,
                                       size = base_size + 3, family = "sans"),
      axis.text.y       = element_text(size = base_size +3, family = "sans"),
      panel.grid        = element_blank(),
      #removes background grid
      legend.title      = element_text(size = base_size + 3, family = "sans"),
      legend.text       = element_text(size = base_size + 2, family = "sans"),
      legend.key.height = unit(1.4, "cm"),
      legend.key.width  = unit(0.45, "cm"),
      #makes sure the legend is correctly formatted
      plot.margin       = margin(10, 10, 10, 10)
      #add some space around the plot
    )
}
#formats the heatmap's style

plot_corr_values <- function(mat, title, global_lim = NULL) {
  tri        <- lower_tri_na(mat)
  df         <- melt_lower(tri)
  #calls functions
  df$row_var <- factor(df$row_var, levels = rev(levels(df$row_var)))
  #flips it so the triangle is with a top left orientation 
  
  lim <- if (!is.null(global_lim)) global_lim else
    ceiling(max(abs(df$value), na.rm = TRUE) * 10) / 10
  #sets the colour scale
  #uses the golbal limit so its the same across all plots 
  #if golbabl limit doesnt exit, compute it by taking the absolut maximum value
  #rounds up to nearest 0.1
  
  raw_label        <- sprintf("%.3f", df$value)
  #rounds the 3 significant figures in the tau values
  df$label         <- ifelse(raw_label == "-0.000", "0.000", raw_label)
  #if rounded to two decimal places goves -0.000 change to just 0.00
  df$text_colour   <- "black"
  #chooses text colou
  
  ggplot(df, aes(x = col_var, y = row_var, fill = value)) +
    #builds plot
    geom_tile(colour = "white", linewidth = 0.5) +
    #creates tiles
    geom_text(aes(label = label, colour = text_colour), size = 4,
              family = "sans") +
    #adds the correlation value inside the tile
    scale_colour_identity() +
    #use colours directly
    coord_fixed() +
    #aspect ratio
    scale_x_discrete(labels = apply_labels) +
    scale_y_discrete(labels = apply_labels) +
    #applys nice axis labels
    scale_fill_gradient2(
      #creates colour gradient
      low = "#d73027", mid = "#ffffff", high = "#4575b4",
      midpoint = 0, limits = c(-lim, lim), name = "\u03c4",
      #the name bit is the code for the tau symbol 
      guide = guide_colourbar(
        #for appearance of legends
        nbin        = 300,
        #lots pf bins mean smooth gradient
        draw.ulim   = TRUE,
        draw.llim   = TRUE,
        #draws arrows if values are outside the bar
        label.theme = element_text(size = 13, family = "sans"),
        barheight   = unit(8, "cm"),
        barwidth    = unit(0.8, "cm")
      )
    ) +
    labs(title = title, x = NULL, y = NULL) +
    #adds title and removes axis labels for the legend 
    heatmap_theme()
    #applies hte previously defined theme
}
#plotting for correlation heat map

plot_pval <- function(mat, title) {
  tri        <- lower_tri_na(mat)
  df         <- melt_lower(tri)
  #applies previously defined functions
  df$row_var <- factor(df$row_var, levels = rev(levels(df$row_var)))
  #reverts the y axis to keep the top left visual 
  
  df$stars <- sig_stars(df$value)
  #converts p values to stars for labelling 
  df$plab  <- format_pval(df$value)
  #uses earlier function to decided decimal or scientific notation
  df$label <- paste0(df$plab, "\n", df$stars)
  #combines the two into one label made of two lines
  
  df$text_colour <- "black"
  #chooses colour for text
  
  caption_text <- paste0(
    "*** p<", ALPHA_THREE_STAR,
    "   ** p<", ALPHA_TWO_STAR,
    "   * p<", ALPHA_ONE_STAR,
    "   ns p\u2265", ALPHA_ONE_STAR
  )
  #creates the legend for the star system
  
  pval_breaks <- c(0, 0.25, 0.5, 0.75, 1)
  pval_labels <- c("0", "0.25", "0.50", "0.75", "1.00")
  #defines the labels for the colour bar so that they dont over crowd
  
  #now to building of the plot
  ggplot(df, aes(x = col_var, y = row_var, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    #creates the tiles
    geom_text(aes(label = label, colour = text_colour),
              size = 3.8, lineheight = 1.2, family = "sans") +
    #sets setting for text in the tiles
    scale_colour_identity() +
    #uses text colour directly
    coord_fixed() +
    #aspect ratio
    scale_x_discrete(labels = apply_labels) +
    scale_y_discrete(labels = apply_labels) +
    #applays nice axis labels
    scale_fill_gradientn(
      colours = c("#2c7bb6",   #for most significant
                  "#abd9e9",   #for values around 0.05
                  "#d73027"),  #for not significant
      values  = scales::rescale(c(0, ALPHA_ONE_STAR, 1)),
      #maps the colurs to the p values
      limits  = c(0, 1),
      #ensures it stays within the limits
      name    = "p-value",
      #legend title
      breaks  = pval_breaks,
      labels  = pval_labels,
      #uses the breaks defined earlier
      guide   = guide_colourbar(
        nbin          = 300,
        #lots of bins for smooth gradient
        draw.ulim     = TRUE,
        draw.llim     = TRUE,
        label.theme   = element_text(size = 13, family = "sans"),
        barheight     = unit(8, "cm"),
        barwidth      = unit(0.8, "cm")
        #basically same as before from th correlation legend
      )
    ) +
    labs(
      title   = title,
      x       = NULL,
      y       = NULL,
      caption = caption_text
      #adds title and captions and removes axis labels from legend bar
    ) +
    heatmap_theme()
  #applies theme
}
#plots the p-value heat map

save_composite <- function(p, filename, w = 12, h = 18) {
  ggsave(
    #sets consistent width and height
    filename = file.path(output_dir, filename),
    #combines the output defined with the file name (defined later)
    plot     = p,
    width    = w,
    height   = h,
    dpi      = 300,
    bg       = "white"
  )
  message("  Saved: ", filename)
  #confirmation message 
}
#function for saving the files 

all_tau_vals <- c()
#collects all the tau values
for (model in model_names) {
  files <- file_list[[model]]
  #extracts files related to the current model
  for (key in c("tau", "partial_tau")) {
    #cycles through tau and partial tau values
    f <- files[[key]]
    #get file path for current file
    if (file.exists(f)) {
      #make sure it exists
      mat <- read_matrix(f)
      #reads the csv into a matrix
      tri <- mat
      #copies matrix into here
      tri[upper.tri(tri, diag = TRUE)] <- NA
      #removes upper triangle and diagonal
      all_tau_vals <- c(all_tau_vals, as.vector(tri[!is.na(tri)]))
      #converts all non-na values into vector and appends
    }
  }
}
global_tau_lim <- ceiling(max(abs(all_tau_vals), na.rm = TRUE) * 10) / 10
#computes the golbal colour scale for the legends, rounds up to nearest 0.1

for (model in model_names) {
  #loop through the models
  message("\n=== ", model, " ===")
  #add message in console on which one it is working
  files <- file_list[[model]]
  #gets file paths for current model
  
  if (file.exists(files$tau) && file.exists(files$pval)) {
    #if both files exist, proceed
    tau_mat  <- read_matrix(files$tau)
    pval_mat <- read_matrix(files$pval)
    #load in the data
    
    p_tau  <- plot_corr_values(tau_mat,  NULL,
                               global_lim = global_tau_lim)
    p_pval <- plot_pval(pval_mat,        NULL)
    #creates corerlation and p value plots
    
    composite <- ggarrange(p_tau, p_pval,
                           ncol   = 1, nrow = 2,
                           labels = c("a)", "b)"),
                           font.label = list(size = 20, face = "bold", family = "sans"))
    #arranges the two plots into one composite
    
    save_composite(composite, paste0(model, "_kendall_composite.png"))
  } else {
    if (!file.exists(files$tau))  warning("Missing: ", files$tau)
    if (!file.exists(files$pval)) warning("Missing: ", files$pval)
  }
  #saves both files 
  #if either is missing, prints out a warning
  
  #repeats for the partial tau
  if (file.exists(files$partial_tau) && file.exists(files$partial_pval)) {
    ptau_mat  <- read_matrix(files$partial_tau)
    ppval_mat <- read_matrix(files$partial_pval)
    
    p_ptau  <- plot_corr_values(ptau_mat,  NULL,
                                global_lim = global_tau_lim)
    p_ppval <- plot_pval(ppval_mat,         NULL)
    
    composite <- ggarrange(p_ptau, p_ppval,
                           ncol   = 1, nrow = 2,
                           labels = c("a)", "b)"),
                           font.label = list(size = 20, face = "bold", family = "sans"))
    
    save_composite(composite, paste0(model, "_partial_composite.png"))
  } else {
    if (!file.exists(files$partial_tau))  warning("Missing: ", files$partial_tau)
    if (!file.exists(files$partial_pval)) warning("Missing: ", files$partial_pval)
  }
}