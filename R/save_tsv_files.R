save_tsv_files <- function(defs){
  save_dir <- gsub("//", "/",
                   paste0(defs$output.dir, "/tsv_files"),
                   fixed = TRUE)

  cat("\nSaving results to files: ")
  if (!dir.exists(save_dir)) dir.create(save_dir,
                                        recursive = TRUE,
                                        showWarnings = FALSE)


  print_results_correlations(correlations = defs$contrasts.corrected,
                             annotation   = defs$annotation.contrasts,
                             outputName   = paste0(save_dir,
                                                   "/contrasts_corrected.tsv"),
                             type         = "q_value")

  print_results_correlations(correlations = defs$contrasts,
                             annotation   = defs$annotation.contrasts,
                             outputName   = paste0(save_dir,
                                                   "/contrasts_raw.tsv"),
                             type         = "correlation")

  print_results_correlations(correlations = defs$sum,
                             annotation   = defs$annotation.sum,
                             outputName   = paste0(save_dir,
                                                   "/sum.tsv"),
                             type         = "sum")

  print_results_correlations(correlations = defs$sd,
                             annotation   = defs$annotation.sd,
                             outputName   = paste0(save_dir,
                                                   "/sd.tsv"),
                             type         = "sd")

  print_results_correlations(correlations = defs$cv,
                             annotation   = defs$annotation.cv,
                             outputName   = paste0(save_dir,
                                                   "/cv.tsv"),
                             type         = "cv")


  fieldnames <- c(paste0("correlations.",
                         c("pearson", "spearman", "kendall")),
                  paste0("results.correlations.pvalue.",
                         c("pearson", "spearman", "kendall")))

  filenames  <- paste0(save_dir, "/",
                       c(paste0(c("p", "s", "k"), "_corr_results.tsv"),
                         paste0(c("p", "s", "k"), "_corr_qvalues_results.tsv")))

  types     <- c(rep("correlation", 3), rep("q_value", 3))

  for (i in seq_along(fieldnames)){
    print_results_correlations(correlations = defs[[fieldnames[i]]],
                               annotation   = defs$annotation.cor,
                               outputName   = filenames[[i]],
                               type         = types[i])
  }

  cat(" done!\n")

  invisible(defs)
}
