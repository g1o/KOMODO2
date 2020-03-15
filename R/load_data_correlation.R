# Specific function to load data if type == "correlation" in the KOMODO2
# definition list. (Not exported to to the namespace)

load_data_correlation <- function(defs){

  # ================== Sanity checks ==================
  # TODO: write check function
  # defs <- check_inputs_correlation(defs)


  # ================== Process test.path  ==================

  defs$y.name <- utils::read.csv(defs$dataset.info,
                                 header       = FALSE,
                                 strip.white  = TRUE,
                                 comment.char = "",
                                 check.names  = FALSE,
                                 sep          = "\t",
                                 stringsAsFactors = FALSE)

  defs$x <- defs$y.name[, defs$x.column]

  if (defs$denominator.column == "") {
    defs$denominator <- ""
  } else {
    defs$denominator <- defs$y.name[, defs$denominator.column]
  }

  defs$y.name <- paste0(defs$annotation.files.dir, "/", defs$y.name[, 1])
  defs$y.name <- gsub(pattern = "//", replacement = "/", x = defs$y.name,
                      fixed = TRUE)

  cat("\nLoading data:\n")
  defs$y <- pbmcapply::pbmclapply(X              = defs$y.name,
                                  FUN            = utils::read.csv,
                                  sep            = "\t",
                                  header         = TRUE,
                                  colClasses     = "character",
                                  strip.white    = TRUE,
                                  comment.char   = "",
                                  check.names    = FALSE,
                                  mc.preschedule = FALSE,
                                  mc.cores       = defs$cores)

  return(defs)
}
