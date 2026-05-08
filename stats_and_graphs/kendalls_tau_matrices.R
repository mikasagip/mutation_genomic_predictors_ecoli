library(dplyr)
library(tidyr)
library(ppcor)
#load in libraries

df <- read.csv("all_properties.csv", stringsAsFactors = FALSE)
#imports in the data, maintains test as strings and not as factors

rownames(df) <- df$locus_tag
#set locus_tag as row identifier

cols_to_drop <- c(
  "mutation_count", "start", "middle", "end",
  "strand", "type", "asymmetry",
  "essentiality_gerdes", "essentiality_pec",
  "distance_from_oric",
  "gene_names"        
)
#defines columns to always drop (these are either non-numeric or not relevant)

df_working <- df %>% dplyr::select(-dplyr::any_of(cols_to_drop))
#save all other columns into new dataframe

core_cols <- c(
  "mutation_density",
  "avg_length",
  "binary_asymmetry",
  "relative_distance_from_oric",
  "transcript",
  "gc_content",
  "cpg_excess"
)
#defines core columns that will be in every model version 

model_cols <- list(
  model1 = c(core_cols,
             "relative_distance_oscillating_triangle",
             "binary_pec",
             "protein_abundance_pax"),
  
  model2 = c(core_cols,
             "relative_distance_oscillating_cosine",
             "binary_pec",
             "protein_abundance_pax"),
  
  model3 = c(core_cols,
             "relative_distance_oscillating_triangle",
             "binary_gerdes",
             "protein_abundance_pax"),
  
  model4 = c(core_cols,
             "relative_distance_oscillating_triangle",
             "binary_pec",
             "protein_abundance_lb")
)
#defines the columns to be used in every model variation (where model 1 is the default)

kendall_matrix <- function(data) {
  cols <- colnames(data)
  #saves columns
  n <- length(cols)
  #number of columns
  mat_tau <- matrix(NA, nrow = n, ncol = n, dimnames = list(cols, cols))
  mat_p   <- matrix(NA, nrow = n, ncol = n, dimnames = list(cols, cols))
  #creates empty nxn matrices for correlation and p-values
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        #loops through every column row combination
        mat_tau[i, j] <- 1
        mat_p[i, j]   <- 0
        #if same column title, corrleation =1 and p =0 becayse they are the same
      } else if (j > i) {
        #if not the same name, and the combination was not already tested
        test <- cor.test(data[[cols[i]]], data[[cols[j]]],
                         method = "kendall", exact = FALSE)
        #does kendalls tau correlation
        #exact = false means it can handle ties
        mat_tau[i, j] <- test$estimate
        mat_tau[j, i] <- test$estimate
        #return list of resulats for correlation
        mat_p[i, j]   <- test$p.value
        mat_p[j, i]   <- test$p.value
        #return list of results for p values
        
      }
    }
  }
  list(
    tau = as.data.frame(mat_tau),
    p   = as.data.frame(mat_p)
  )
  #converts the matrix frame into a dataframe
}

#kendalls tau correlation matrix function 


kendall_partial_matrix <- function(data) {
  result <- tryCatch(
    pcor(data, method = "kendall"),
    error = function(e) {
      cat(sprintf("  WARNING: Partial correlation failed — %s\n", e$message))
      return(NULL)
    }
  )
  
  if (is.null(result)) return(NULL)
  
  list(
    tau = as.data.frame(result$estimate),
    p   = as.data.frame(result$p.value)
  )
}
#creates kendalls partial correlation matrix function

results <- list()
#make list for results

for (model_name in names(model_cols)) {
  #iterates throught the 4 models
  
  cat(sprintf("\nRunning %s...\n", model_name))
  #prints in the console which model is running
  
  cols <- model_cols[[model_name]]
  #gets the columns
  model_data <- df_working %>%
    dplyr::select(all_of(cols)) %>%
    mutate(across(everything(), as.numeric))
  #selectes the columns and makes sure they are numeric
  
  #check column integrity
  if (any(duplicated(colnames(model_data)))) {
    stop("Duplicate column names detected!")
  }
  
  #ensure no columns became all NA
  bad_cols <- colnames(model_data)[colSums(is.na(model_data)) == nrow(model_data)]
  if (length(bad_cols) > 0) {
    stop(paste("Columns became all NA after numeric conversion:", paste(bad_cols, collapse = ", ")))
  }
  
  n_before <- nrow(model_data)
  model_data <- model_data %>% na.omit()
  n_after <- nrow(model_data)
  cat(sprintf("  Rows used: %d (dropped %d with NA)\n", n_after, n_before - n_after))
  #removes any rows with NA and reports this in the console
  
  matrices <- kendall_matrix(model_data)
  #compute kendalls tau matrices
  
  partial_matrices <- kendall_partial_matrix(model_data)
  #computes kendalls partial tau matrices
  
  results[[model_name]] <- list(full = matrices, partial = partial_matrices)
  #save the result in the list
  
  out_tau <- sprintf("%s_kendall_tau.csv", model_name)
  write.csv(matrices$tau, file = out_tau, row.names = TRUE)
  cat(sprintf("  Saved: %s\n", out_tau))
  
  out_p <- sprintf("%s_kendall_pvalues.csv", model_name)
  
  p_out <- matrices$p
  diag(p_out) <- NA
  p_out[p_out == 0] <- 2e-308
  diag(p_out) <- 0
  #if p=0 and is not on the diagonal, writes out as 2e-103 instead 
  
  write.csv(p_out, file = out_p, row.names = TRUE)
  cat(sprintf("  Saved: %s\n", out_p))
  
  if (!is.null(partial_matrices)) {
    #makes sure it wont error incase partial matrix fails 
    out_partial_tau <- sprintf("%s_kendall_partial_tau.csv", model_name)
    write.csv(partial_matrices$tau, file = out_partial_tau, row.names = TRUE)
    cat(sprintf("  Saved: %s\n", out_partial_tau))
    #writes out the partial matrix tau
    
    out_partial_p <- sprintf("%s_kendall_partial_pvalues.csv", model_name)
    partial_p_out <- partial_matrices$p
    diag(partial_p_out) <- NA
    partial_p_out[partial_p_out == 0] <- 2e-308
    diag(partial_p_out) <- 0
    #changes all 0 p-values not on diagonal with <2e-308
    write.csv(partial_p_out, file = out_partial_p, row.names = TRUE)
    cat(sprintf("  Saved: %s\n", out_partial_p))
    #writes out the partial matrix p-values
  }
}
#runs the correlation matrix function for each of the models

for (model_name in names(results)) {
  cat(sprintf("\n========== %s Kendall's Tau Matrix ==========\n", toupper(model_name)))
  print(round(results[[model_name]]$full$tau, 4))
  cat(sprintf("\n========== %s P-Value Matrix ==========\n", toupper(model_name)))
  print(round(results[[model_name]]$full$p, 4))
  if (!is.null(results[[model_name]]$partial)) {
    cat(sprintf("\n========== %s Kendall's Partial Tau Matrix ==========\n", toupper(model_name)))
    print(round(results[[model_name]]$partial$tau, 4))
    cat(sprintf("\n========== %s Partial P-Value Matrix ==========\n", toupper(model_name)))
    print(round(results[[model_name]]$partial$p, 4))
  }
}

#prints off the results (tau and p value) for each matrix 

cat("\nDone. All correlation matrices saved as CSVs.\n")