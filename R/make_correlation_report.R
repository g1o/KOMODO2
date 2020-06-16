make_correlation_report <- function(defs){

  # ================== Sanity checks ==================
  #assertthat::assert_that("linear.model.cutoff" %in% names(defs),
  #                        is.numeric(defs$linear.model.cutoff))

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

  cat("\nGenerating HTML5 report for results")
  cat("\n-----------------------------------")
  cat("\nUsing filters:\n")
  cat("\ncorrected_contrasts < ", defs$linear_model.qvalue.cutoff)
  cat("\nSpearman_qvalue     < ", defs$spearman.qvalue.cutoff)
  cat("\nPearson_qvalue      < ", defs$pearson.qvalue.cutoff)
  cat("\nKendall_qvalue      < ", defs$kendall.qvalue.cutoff)
  cat("\nPearson_cor         > ", defs$pearson.cor.upper.cutoff)
  cat("\nPearson_cor         < ", defs$pearson.cor.lower.cutoff)
  cat("\nSpearman_cor        > ", defs$spearman.cor.upper.cutoff)
  cat("\nSpearman_cor        < ", defs$spearman.cor.lower.cutoff)
  cat("\nKendall_cor         > ", defs$kendall.cor.upper.cutoff)
  cat("\nKendall_cor         < ", defs$kendall.cor.lower.cutoff)
  cat("\nsize                > ", defs$annotation_size.cutoff)
  cat("\nprevalence          > ", defs$prevalence.cutoff)
  cat("\nheterogeneity       > ", defs$heterogeneity.cutoff)
  cat("\nsd                  > ", defs$sd.cutoff)
  cat("\ncv                  > ", defs$cv.cutoff)
  cat("\n-----------------------------------")
  cat("\nThis may take a while...")


  # Prepare output folder

  #  od <- defs$output.dir

  #Uncomment the line below to generate full paths. Results may not work in servers.
  od  <- normalizePath(defs$output.dir)

  cpd <- gsub("//", "/", paste0(od, "/correlation_Plots/"),
              fixed = TRUE)

  if(!dir.exists(cpd)) dir.create(cpd, recursive = TRUE)

  # Copy report template into output dir
  files_to_copy <- dir(system.file("extdata", package = "KOMODO2"))
  fp <- files_to_copy

  for (i in seq_along(files_to_copy)){
    fp[i] <- gsub("//", "/", paste0(defs$output.dir, "/", files_to_copy[i]),
                  fixed = TRUE)
    file.copy(system.file("extdata", files_to_copy[i], package = "KOMODO2"),
              to = fp[i], overwrite = TRUE)
  }

  suppressWarnings(rmarkdown::render_site(input = defs$output.dir,
                                          quiet = TRUE))
  file.remove(fp)

  # Invoke browser andopen results
  myURL <- gsub("//", "/", paste0(defs$output.dir, "/index.html"), fixed = TRUE)
  myURL <- paste0("file:/",
                  normalizePath(gsub("./", "", myURL, fixed = TRUE)))
  utils::browseURL(myURL)
  cat("\nAnd we're done!")

  invisible(defs)
}
