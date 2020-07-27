# Perform analysis for the KOMODO2 workflow
#
# This script implements the third step of the LCFD workflow of KOMODO2.
# It is responsible for performing the actual analysis and generating
# all the processed data and tables.
#
# The script expects enriched `KOMODO2`-type lists generated after running
# [load_data()] -> [clean_data()].
#
#
# @param defs an enriched KOMODO2-type list object (see Details).
#


do_analysis <- function(defs){

  # ================== Sanity checks ==================
  assert_that("tree.type" %in% names(defs),
              "tree.path" %in% names(defs))

  if("MHT.method" %in% names(defs)){
    assertthat::assert_that(defs$MHT.method %in% stats::p.adjust.methods)
  } else {
    defs$MHT.method <- "BH"
  }

  # Perform analysis
  #cat("\nMain Analysis:\n")
  defs <- switch(tolower(defs$type),
                 significance = do_analysis_significance(defs),
                 correlation  = do_analysis_correlation(defs))


}
