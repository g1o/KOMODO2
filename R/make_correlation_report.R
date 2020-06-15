make_correlation_report <- function(defs){

  # ================== Sanity checks ==================
  assertthat::assert_that("linear.model.cutoff" %in% names(defs),
                          is.numeric(defs$linear.model.cutoff))

  # "activate" packages that are used only in the .Rmd report generation.
  # (needed only to stop a NOTE in R CMD check)
  tmp <- dendextend::fac2num(factor(3:5))
  tmp <- heatmaply::BrBG(5)

  # ================== Prepare report ==================
  sumY    <- sapply(defs$y, sum)  # faster than apply() or colSums()!

  # filter out those with no observations (sum equals zero)
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
#                     defs$prevalence_per_group[order(names(defs$prevalence_per_group))],
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
#                                 "prevalence per group",
                                 "heterogeneity")

  plotframe <- as.data.frame(t(plotframe))

  description           <- unlist(defs$annotation.cor)
  plotframe$description <- description[order(names(description))]
  plotframe$name        <- rownames(plotframe)

  idx <- which(
      # q-value cutoffs
      (plotframe$corrected_contrasts < defs$linear_model.qvalue.cutoff) &
      (plotframe$Spearman_qvalue   < defs$spearman.qvalue.cutoff) &
      (plotframe$Pearson_qvalue    < defs$pearson.qvalue.cutoff) &
      (plotframe$Kendall_qvalue    < defs$kendall.qvalue.cutoff) &
      # correlation cutoffs
      (plotframe$Pearson_cor       > defs$pearson.cor.upper.cutoff) &
      (plotframe$Pearson_cor       < defs$pearson.cor.lower.cutoff) &
      (plotframe$Spearman_cor      > defs$spearman.cor.upper.cutoff) &
      (plotframe$Spearman_cor      < defs$spearman.cor.lower.cutoff) &
      (plotframe$Kendall_cor       > defs$kendall.cor.upper.cutoff) &
      (plotframe$Kendall_cor       < defs$kendall.cor.lower.cutoff) &
      # basic statistics cutoffs
      (plotframe$size              > defs$annotation_size.cutoff) &
      (plotframe$prevalence        > defs$prevalence.cutoff) &
      (plotframe$heterogeneity     > defs$heterogeneity.cutoff) &
      (plotframe$sd                > defs$sd.cutoff) &
      (plotframe$cv                > defs$cv.cutoff))

  df_cutoff <- plotframe[idx, ]

  if (isTRUE(defs$raw_data_sd_filter)) {
    # remove trivial cases, constant values.
    df_cutoff <- df_cutoff[df_cutoff$sd != 0, ]
  }

  defs$sig_IDs <- rownames(df_cutoff)

  cat("\nGenerating HTML5 report for results with\nphylogeny-aware q-values < ",
      defs$linear_model.qvalue.cutoff, "\n(this may take a while).")
  # Prepare output folder

#  od <- defs$output.dir

  #Uncomment the line below to generate full paths. Results may not work in servers.
  od  <- normalizePath(defs$output.dir)

  cpd <- gsub("//", "/", paste0(od, "/correlation_Plots/"),
              fixed = TRUE)

  if(!dir.exists(cpd)) dir.create(cpd, recursive = TRUE)

  # Copy report template into output dir
  fp1 <- gsub("//", "/", paste0((defs$output.dir),
                               "/_site.yml"),
             fixed = TRUE)

  file.copy(system.file("extdata", "_site.yml",
                        package = "KOMODO2"), to = fp1, overwrite = TRUE)

  fp2 <- gsub("//", "/", paste0((defs$output.dir),
                               "/index.Rmd"),
             fixed = TRUE)

  file.copy(system.file("extdata", "index.Rmd",
                        package = "KOMODO2"), to = fp2, overwrite = TRUE)

  fp3 <- gsub("//", "/", paste0((defs$output.dir),
                               "/heatmap_phylo_raw.Rmd"),
                               fixed = TRUE)

  file.copy(system.file("extdata", "heatmap_phylo_raw.Rmd",
                        package = "KOMODO2"), to = fp3, overwrite = TRUE)

  fp4 <- gsub("//", "/", paste0((defs$output.dir),
                                "/heatmap_phylo_norm.Rmd"),
              fixed = TRUE)

  file.copy(system.file("extdata", "heatmap_phylo_norm.Rmd",
                        package = "KOMODO2"), to = fp4, overwrite = TRUE)

  fp5 <- gsub("//", "/", paste0((defs$output.dir),
                                "/heatmap_phylo_perc.Rmd"),
              fixed = TRUE)

  file.copy(system.file("extdata", "heatmap_phylo_perc.Rmd",
                        package = "KOMODO2"), to = fp5, overwrite = TRUE)

  fp6 <- gsub("//", "/", paste0((defs$output.dir),
                                "/q_value_scatter.Rmd"),
              fixed = TRUE)

  file.copy(system.file("extdata", "q_value_scatter.Rmd",
                        package = "KOMODO2"), to = fp6, overwrite = TRUE)

  fp7 <- gsub("//", "/", paste0((defs$output.dir),
                                "/table.Rmd"),
              fixed = TRUE)

  file.copy(system.file("extdata", "table.Rmd",
                        package = "KOMODO2"), to = fp7, overwrite = TRUE)

  suppressMessages(
    suppressWarnings(rmarkdown::render_site(input = defs$output.dir,
                                            quiet = TRUE)))

#  suppressWarnings(rmarkdown::render(fp,
#                                     output_file = "KOMODO2_report.html",
#                                     quiet = TRUE))

#  suppressWarnings(rmarkdown::render(fp2,
#                                     output_file = "heatmap_phylo.html",
#                                     quiet = TRUE))


  file.remove(fp1)
  file.remove(fp2)
  file.remove(fp3)
  file.remove(fp4)
  file.remove(fp5)
  file.remove(fp6)
  file.remove(fp7)

  cat("done")

  invisible(defs)
}
