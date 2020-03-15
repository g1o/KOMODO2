make_correlation_report <- function(defs){

  # ================== Sanity checks ==================
  assertthat::assert_that("linear.model.cutoff" %in% names(defs),
                          is.numeric(defs$linear.model.cutoff))


  # ================== Prepare report ==================
  cutoff <- defs$linear.model.cutoff

  # filter out those with no observations
  sumY    <- sapply(defs$y, sum)  # faster than apply() or colSums()!
  idx     <- (sumY != 0)
  Y       <- defs$y[, idx]
  sumY    <- sumY[idx]
  defs$sd <- defs$sd[idx]
  defs$cv <- defs$cv[idx]

  defs$greaterthanzero <- defs$greaterthanzero[idx]
  defs$heterogeneity   <- defs$heterogeneity[idx]

  # Prepare plotframe
  plotframe <- rbind(defs$contrasts.corrected[order(names(defs$contrasts.corrected))],
                     defs$results.correlations.pvalue.pearson[order(names(defs$results.correlations.pvalue.pearson))],
                     defs$results.correlations.pvalue.spearman[order(names(defs$results.correlations.pvalue.spearman))],
                     defs$results.correlations.pvalue.kendall[order(names(defs$results.correlations.pvalue.kendall))],
                     defs$correlations.pearson[order(names(defs$correlations.pearson))],
                     defs$correlations.spearman[order(names(defs$correlations.spearman))],
                     defs$correlations.kendall[order(names(defs$correlations.kendall))],
                     sumY[order(names(sumY))],
                     defs$sd[order(names(defs$sd))],
                     defs$cv[order(names(defs$cv))],
                     defs$greaterthanzero[order(names(defs$greaterthanzero))],
                     defs$heterogeneity[order(names(defs$heterogeneity))])
  plotframe <- rbind(plotframe,
                     Y[, order(colnames(Y))])

  rownames(plotframe)[1:12] <- c("corrected_contrasts",
                                 "Pearson_qvalue",
                                 "Spearman_qvalue",
                                 "Kendall_qvalue",
                                 "Pearson_cor",
                                 "Spearman_cor",
                                 "Kendall_cor",
                                 "size",
                                 "sd",
                                 "cv",
                                 "prevalence",
                                 "heterogeneity")

  plotframe <- as.data.frame(t(plotframe))

  description           <- unlist(defs$annotation.cor)
  plotframe$description <- description[order(names(description))]
  plotframe$name        <- rownames(plotframe)

  df_cutoff <- plotframe[plotframe$corrected_contrasts < cutoff, ]
  df_cutoff <- df_cutoff[df_cutoff$sd != 0, ] # remove trivial cases, constant values.

  cat("\nGenerating HTML5 report for results with\nphylogeny-aware q-values < ",
      cutoff, "\n(this may take a while).")
  # Prepare output folder
  od  <- normalizePath(defs$output.dir)
  cpd <- gsub("//", "/", paste0(od, "/correlation_Plots/"),
              fixed = TRUE)
  if(!dir.exists(cpd)) dir.create(cpd, recursive = TRUE)

  # Copy report template into output dir
  fp <- gsub("//", "/", paste0(normalizePath(defs$output.dir),
                               "/K2rep.Rmd"),
             fixed = TRUE)
  file.copy(system.file("extdata", "KOMODO2_correlation_report.Rmd",
                        package = "KOMODO2"), to = fp, overwrite = TRUE)

  suppressWarnings(rmarkdown::render(fp,
                                     output_file = "KOMODO2_report.html",
                                     quiet = TRUE))
  file.remove(fp)
  cat("done")

  invisible(defs)
}
