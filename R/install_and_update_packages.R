#' Install and update dependencies
#'
#' This function installs the latest versions of all packages
#' listed in the `Imports` and `Suggests` fields of KOMODO2.
#' It uses the standard `utils::install.packages()` for CRAN
#' packages, and `BiocManager::install()` for Bioconductor
#' dependencies. Further arguments to both functions are passed
#' as lists.
#'
#' @param which which dependencies to install/update? Accepts
#' "cran", "bioc" or "all".
#' @param cran.args list containing further arguments to be passed
#' down to `install.packages()`.
#' @param bioc.args list containing further arguments to be passed
#' down to `BiocManager::install()`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   install_and_update_packages(which = "all")
#' }

install_and_update_packages <- function(which = "all",
                                        cran.args = list(),
                                        bioc.args = list()){
  # ================== Sanity checks ==================
  assertthat::assert_that(is.character(which),
                          length(which) == 1,
                          which %in% c("cran", "bioc", "all"),
                          is.list(cran.args),
                          is.list(bioc.args))

  cran.args$pkgs <- c("assertthat",
                      "pbmcapply",
                      "ape",
                      "rmarkdown",
                      "nlme",
                      "BiocManager",
                      "taxize",
                      "ggplot2",
                      "plotly",
                      "DT",
                      "htmltools",
                      "htmlwidgets",
                      "pkgdown",
                      "knitr",
                      "dendextend",
                      "heatmaply")
  bioc.args$pkgs <- c("AnnotationDbi",
                      "KEGGREST",
                      "GO.db")



  if (which %in% c("cran", "all")){
    do.call(utils::install.packages, cran.args)
  }

  if (which %in% c("bioc", "all")){
    do.call(BiocManager::install, bioc.args)
  }
}
