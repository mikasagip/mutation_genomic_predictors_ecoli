library(ggplot2)
library(dplyr)
library(ggpubr)

df <- read.csv("all_properties.csv", stringsAsFactors = FALSE)
#imports in the data, maintains text as strings and not as factors

df_mut <- df %>%
  mutate(has_mutation = mutation_density > 0)
#creates a dataframe for just those genes that have had any mutatuion

tab <- table(df_mut$essentiality_pec, df_mut$has_mutation)
#makes them into a table

chisq_res <- chisq.test(tab, correct = FALSE)
chisq_res$p.value
chisq_res
#does a chisqr test on between the E and NE conditions 
#prints the result

pval_df <- data.frame(
  essentiality_pec = c("E", "NE"),         
  y.position = 1.05,                        
  label = ifelse(chisq_res$p.value < 0.001, "***",
                 ifelse(chisq_res$p.value < 0.01, "**",
                        ifelse(chisq_res$p.value < 0.05, "*", "ns")))
)
#creates the tags for significance

plot_a <- ggplot(df_mut, aes(x = essentiality_pec, fill = has_mutation)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("TRUE" = "dodgerblue3", "FALSE" = "red3")) +
  labs(
    x = "PEC Gene Essentiality",
    y = "Proportion",
    fill = "Mutation present"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold", size = 14, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(face = "bold", size = 10, family = "sans"),
    legend.text = element_text(size = 10, family = "sans"),
    legend.position = "top"
  ) +
  geom_text(
    data = pval_df,
    aes(x = essentiality_pec, y = y.position, label = label),
    inherit.aes = FALSE
  )
#plots out in a stacked bar based on presence or absence of any mutations 
#basically shows the propertion of data which is 0
#adds the significance tag from the chisqr test

  
df_nonzero <- df %>%
  filter(mutation_density > 0)
#creates a data frame for only those points with non-zero mutation density


plot_b <- ggplot(df_nonzero, aes(x = essentiality_pec, y = mutation_density, fill = essentiality_pec)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  scale_fill_manual(values = c("E" = "gold", "NE" = "darkorchid")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("E", "NE")),
    label.y = max(df_nonzero$mutation_density) * 1.05
  ) +
  labs(
    x = "PEC Gene Essentiality",
    y = "Mutation Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.position = "top",
    legend.text = element_text(color = "white"),
    legend.key = element_rect(fill = "white", color = "white")
  ) +
  guides(fill = guide_legend(title = "", override.aes = list(fill = "white", color = "white")))
#creates a violin plot based on the essentilaity data only non-zero values
#does a wilcoxon test and labels the graph with stars or ns depending on result


wilcox.test(mutation_density ~ essentiality_pec, data = df_nonzero)
#also runs a wilcoxon test seperatly so results can be reported

outlier_info <- df_nonzero %>%
  group_by(essentiality_pec) %>%
  summarise(
    Q1 = quantile(mutation_density, 0.25),
    Q3 = quantile(mutation_density, 0.75)
  ) %>%
  rowwise() %>%
  mutate(
    IQR = Q3 - Q1,
    lower = Q1 - 1.5*IQR,
    upper = Q3 + 1.5*IQR
  )
#defines the outliers and finds them

df_outliers <- df_nonzero %>%
  left_join(outlier_info, by = "essentiality_pec") %>%
  mutate(is_outlier = mutation_density < lower | mutation_density > upper)
#makes new data frame for outliers

tab_outliers <- table(df_outliers$essentiality_pec, df_outliers$is_outlier)
#makes into table for analysis

outlier_chi <- chisq.test(tab_outliers, correct = FALSE)
outlier_chi
#does chisqr test and prints the result

pval_df_out <- data.frame(
  essentiality_pec = c("E", "NE"),
  y.position = 1.05,
  label = ifelse(outlier_chi$p.value < 0.001, "***",
                 ifelse(outlier_chi$p.value < 0.01, "**",
                        ifelse(outlier_chi$p.value < 0.05, "*", "ns")))
)
#defines p-value ranges for labelling on graph

plot_c <- ggplot(df_outliers, aes(x = essentiality_pec, fill = is_outlier)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("TRUE" = "dodgerblue3", "FALSE" = "red3")) +
  labs(
    x = "PEC Gene Essentiality",
    y = "Proportion",
    fill = "Outlier"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold", size = 14, family = "sans"),
    axis.text = element_text(size = 10, family = "sans"),
    legend.title = element_text(face = "bold", size = 10, family = "sans"),
    legend.text = element_text(size = 10, family = "sans"),
    legend.position = "top"
  ) +
  geom_text(
    data = pval_df_out,
    aes(x = essentiality_pec, y = y.position, label = label),
    inherit.aes = FALSE
  )
#plots proportion of outliers for each groups

composite <- ggarrange(
  plot_a, plot_b, plot_c,
  labels = c("a)", "b)", "c)"),
  #labels them within the figure
  ncol = 3,  
  #two columns
  nrow = 1,  
  #one row
  heights = c(1, 1, 1),
  widths = c(1, 1, 1)
  #makes heights and widths in the figure equal
)
#combines plots into one figure

output_dir <- "figure_outputs"
if (!dir.exists(output_dir)) dir.create(output_dir)

ggsave(
  filename = paste0(output_dir, "/essentiality_pec.png"),
  plot     = composite,
  width    = 10,
  height   = 5,
  dpi      = 300,
  bg       = "white"
)