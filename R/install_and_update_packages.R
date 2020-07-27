#' Install and update dependencies from Suggests field
#'
#' This function installs the latest versions of all packages
#' listed in the `Suggests` fields of KOMODO2.
#' It uses the standard `utils::install.packages()` for CRAN
#' packages, and `BiocManager::install()` for Bioconductor
#' dependencies. Further arguments to both functions are passed
#' as lists.
#'
#' @param cran.args list containing further arguments to be passed
#' down to `install.packages()`.
#' @param bioc.args list containing further arguments to be passed
#' down to `BiocManager::install()`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   install_and_update_packages()
#' }

install_and_update_packages <- function(cran.args = list(),
                                        bioc.args = list()){
  # ================== Sanity checks ==================
  assertthat::assert_that(is.list(cran.args),
                          is.list(bioc.args))

  cran.args$pkgs <- c("ggplot2",
                      "plotly",
                      "DT",
                      "htmltools",
                      "htmlwidgets",
                      "pkgdown",
                      "knitr")
  bioc.args$pkgs <- c("AnnotationDbi",
                      "KEGGREST",
                      "GO.db")

  do.call(utils::install.packages, cran.args)
  do.call(BiocManager::install, bioc.args)
}
