#' Run the KOMODO2 pipeline
#'
#' This script runs the complete LCFD workflow of KOMODO2 and generates the
#' HTML5 output pages and export files.
#'
#' The script expects a `KOMODO2`-type list, which is a list object containing
#' _at least_ the following fields:
#' \itemize{
#'    \item \code{test.path} (char string): path to the folder containing
#'    annotation files of the test group
#'    \item \code{back.path} (char string): path to the folder containing
#'    annotation files of the background group
#'    \item \code{x.path} (char string): path to the file containing
#'    the genomes' attributes (for correlation test)
#'    \item \code{y.path} (char string): path to the folder containing
#'    the the genomes and their annotations (for correlation test)
#'    \item \code{ontology} (char string): which ontology to use. Currently
#'    accepts "GO" or "Gene Ontology", "KEGG" and "other".
#'    \item \code{dict.path} (char string): file with the dictionary (terms and
#'    their meaning) of the ontology, if `ontology` is set as "other".
#'    \item \code{type} (char string): comparison module to use. Currently only
#'    the "correlation" type is implemented.
#' }
#'
#' The input definitions can also be passed as a file path. If that is the
#' case the file must be a text file with a `field = value` format.
#' Blank likes and lines starting with `#` are ignored. Required fields are the
#' same described for the `KOMODO2` list described above.
#'
#' @param defs either a KOMODO2-type list object or
#' a path to a text file containing the required definitions (see Details).
#' @param cores positive integer, how many CPU cores to use (multicore
#' acceleration does not work in Windows systems). Setting
#' this parameter overrides any `cores` field from `defs`. Multicore support is
#' currently implemented using the `parallel` package, which uses forking
#' (which means that multicore support is not available under Windows)
#' @param render.report logical: should a HTML5 report be generated?
#'
#' @importFrom assertthat assert_that is.count has_name
#' @importFrom dplyr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Download required files
#' retrieve_data_files(target.dir = "./data_folder")
#'
#' # Build an input list:
#' defs <- list(annotation.files.dir = "./data_folder/gene2GO",
#'              output.dir           = "./results/GO_Pan_proxy/",
#'              dataset.info         = "./data_folder/metadata/GO_metadata_Pan_proxy.txt",
#'              x.column             = 2,
#'              ontology             = "GO",
#'              dict.path            = "",
#'              column               = "GO",
#'              denominator.column   = "",
#'              tree.path            = "./data_folder/trees/tree_genome_IDs.nwk",
#'              tree.type            = "newick",
#'              linear.model.cutoff  = 0.5,
#'              type                 = "correlation")
#'
#' defs <- run_KOMODO2(defs, cores = parallel::detectCores() - 1)
#' }

# saveRDS(defs, "./results/test_defs.rds")

run_KOMODO2 <- function(defs, cores = NULL, render.report = TRUE){

  # ================== Sanity checks ==================
  assert_that(is.list(defs) || file.exists(defs),
              is.null(cores) || is.count(cores),
              is.logical(render.report),
              length(render.report) == 1,
              msg = "input error(s) in KOMODO2::run_KOMODO2()")


  # ================== Set up parallel processing ==================
  if(!is.null(cores)) {
    defs$cores <- cores
  } else if(!has_name(defs, "cores")) {
    defs$cores <- 1
  }

  if (.Platform$OS.type == "windows"){
    cat("\nMulticore support not currently available for Windows. Forcing cores = 1.")
    defs$cores <- 1
  } else {
    available.cores <- parallel::detectCores()
    if (defs$cores >= available.cores){
      cat("\nAttention: cores too large, we only have ", available.cores,
          " cores.\nUsing ", available.cores - 1,
          " cores for load_data().")
      defs$cores <- available.cores - 1
    }
  }

  defs <- defs %>%
    load_data() %>%       # Load required data
    clean_data() %>%      # Preliminary data cleaning
    do_analysis() %>%     # perform the analysis
    save_tsv_files() %>%  # Save results to .tsv files
    make_report(render.report = render.report) # make the report (if needed)

  invisible(defs)
}
