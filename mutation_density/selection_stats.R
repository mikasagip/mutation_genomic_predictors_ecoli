library(ggplot2)
library(patchwork)
#import in libraries 

rates_foster_pos1 <- read.csv("rates_foster_pos1.csv")
rates_foster_pos2 <- read.csv("rates_foster_pos2.csv")
rates_foster_pos3 <- read.csv("rates_foster_pos3.csv")

rates_wei_pos1 <- read.csv("rates_wei_pos1.csv")
rates_wei_pos2 <- read.csv("rates_wei_pos2.csv")
rates_wei_pos3 <- read.csv("rates_wei_pos3.csv")
#reads in all the data from the csv files

out_df <<- data.frame(
  dataset = character(),
  diff_pos = character(),
  v_statistic = numeric(), 
  p_value = numeric(), 
  stringsAsFactors = FALSE
)

plot_list <<- list()
#list to store the plots

#define function to compare the rates
compare_rates <- function(df1, df2, label1, label2, dataset, diff_pos) {
  
  #match rows by trinucleotide so pairs align correctly
  merged <- merge(df1, df2, by = "trinucleotide",
                  suffixes = c("_1", "_2"))
  
  rate1 <- merged$mutation_rate_1
  rate2 <- merged$mutation_rate_2
  #takes the rates and puts them into lists 
  
  #calculates paired differences
  diffs <- rate1 - rate2
  diffs_df <- data.frame(diffs)
  #put into dataframe
  
  cat("Rows after merge:", nrow(merged), "\n")
  cat("Length of diffs:", length(diffs), "\n")
  #checks that all the rows are present
  
  #prints to console
  cat(label1, "vs", label2, "\n")
  
  file_stem <- paste0(label1, "_vs_", label2)
  #creates file name stem
  
  #normality test on differences
  sh <- shapiro.test(diffs)
  print(sh)
  
  if (sh$p.value > 0.05) {
    
    cat("Data consistent with normality\n")
    
  } else {
    
    cat("Data not consistent with normality\n")
  }
  #prints result and interpretation of p value
  
  result <- wilcox.test(diffs, mu = 0, exact = FALSE)
  #does the wilcoxon paired test as some not normal
  print(result)
  if (result$p.value > 0.05) {
    
    cat("Data not significantly different from 0\n")
    
  } else {
    
    cat("Data significantly different from 0\n")
  }
  
  out_df <<- rbind(out_df, data.frame(
    dataset = dataset,
    diff_pos = diff_pos,
    v_statistic = result$statistic, 
    p_value = result$p.value
  ))
  #write the results of the wilcoxon out to dataframe
  
  p_hist <- ggplot(diffs_df, aes(x=diffs)) +
    geom_histogram(binwidth = diff(range(diffs)) / 15, 
                   fill = colours[[colour_cnt]], 
                   color = "black", 
                   alpha = 0.9) +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(x = "Rate Difference",
         y = "Count"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.title   = element_text(face = "bold"),
      axis.line    = element_line(colour = "grey30")
    )
  
  p_qq <- ggplot(diffs_df, aes(sample = diffs)) +
    stat_qq(color = colours[[colour_cnt]], size = 2, alpha = 0.8) +
    stat_qq_line(color = "grey50", linewidth = 0.8) +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.line  = element_line(colour = "grey30")
    )
  
  key_base <- paste0(dataset, "_", diff_pos) 
  plot_list[[paste0(key_base, "_hist")]] <<- p_hist
  plot_list[[paste0(key_base, "_qq")]]   <<- p_qq
  #save plots for composite figures
  
  colour_cnt <<- colour_cnt + 1
  
}

colours = list("#FF5CCC", "#DC52BF", "#B947B1", "#973DA4", "#743296", "#512889")
#make list of colour codes 

colour_cnt = 1
#start colour counter

#foster comparisons
compare_rates(rates_foster_pos1, rates_foster_pos2,
              "foster_pos1", "foster_pos2", "foster", "1 and 2")

compare_rates(rates_foster_pos1, rates_foster_pos3,
              "foster_pos1", "foster_pos3", "foster", "1 and 3")

compare_rates(rates_foster_pos2, rates_foster_pos3,
              "foster_pos2", "foster_pos3", "foster", "2 and 3")


#wei comparisons
compare_rates(rates_wei_pos1, rates_wei_pos2,
              "wei_pos1", "wei_pos2", "wei", "1 and 2")

compare_rates(rates_wei_pos1, rates_wei_pos3,
              "wei_pos1", "wei_pos3", "wei", "1 and 3")

compare_rates(rates_wei_pos2, rates_wei_pos3,
              "wei_pos2", "wei_pos3", "wei", "2 and 3")

#write results to csv
write.csv(out_df, "selection_test_results_r.csv", row.names = FALSE)

#function to make histogram composite
make_composite_hist <- function(dataset_name1, dataset_name2, outfile) {
  hist1   <- plot_list[[paste0(dataset_name1, "_1 and 2_hist")]]
  hist2   <- plot_list[[paste0(dataset_name2, "_1 and 2_hist")]]
  
  hist3   <- plot_list[[paste0(dataset_name1, "_1 and 3_hist")]]
  hist4   <- plot_list[[paste0(dataset_name2, "_1 and 3_hist")]]
  
  hist5   <- plot_list[[paste0(dataset_name1, "_2 and 3_hist")]]
  hist6   <- plot_list[[paste0(dataset_name2, "_2 and 3_hist")]]
  #retrieves all the relevant plots
  
  composite <- (hist1 | hist2) /
    (hist3 | hist4) /
    (hist5 | hist6)
  #creates configuration
  
  composite <- composite + plot_annotation(
    tag_levels = list(c("a)", "b)", "c)", "d)", "e)", "f)")),
    theme = theme(plot.tag = element_text(face = "bold", size = 30))
  )
  #plots the composite
  
  ggsave(outfile, composite, width = 15, height = 15, dpi = 300, bg = "white")
  #saves the composite
}

#function to make qq composite
make_composite_qq <- function(dataset_name1, dataset_name2, outfile) {
  qq1   <- plot_list[[paste0(dataset_name1, "_1 and 2_qq")]]
  qq2   <- plot_list[[paste0(dataset_name2, "_1 and 2_qq")]]
  
  qq3   <- plot_list[[paste0(dataset_name1, "_1 and 3_qq")]]
  qq4   <- plot_list[[paste0(dataset_name2, "_1 and 3_qq")]]
  
  qq5   <- plot_list[[paste0(dataset_name1, "_2 and 3_qq")]]
  qq6   <- plot_list[[paste0(dataset_name2, "_2 and 3_qq")]]
  #retrieves all the relevant plots
  
  composite <- (qq1 | qq2) /
    (qq3 | qq4) /
    (qq5 | qq6)
  #creates configuration
  
  composite <- composite + plot_annotation(
    tag_levels = list(c("a)", "b)", "c)", "d)", "e)", "f)")),
    theme = theme(plot.tag = element_text(face = "bold", size = 30))
  )
  #plots the composite
  
  ggsave(outfile, composite, width = 15, height = 15, dpi = 300, bg = "white")
  #saves the composite
}

make_composite_hist("foster", "wei", "hist_composite.png")
make_composite_qq("foster", "wei", "qq_composite.png")
#run the composites 