#' Run the KOMODO2 pipeline
#'
#' This function runs the complete workflow of KOMODO2 and generates the
#' HTML5 output pages and export files.
#'
#' The script expects a `KOMODO2`-type list, passed either as an actual list
#' object or as a file path. In the latter case, notice that
#' the file must be a text file with a `field = value` format.
#' Blank likes and lines starting with `#` are ignored. The function expects the
#'  input list to have the following fields:
#' \itemize{
#'    \item\code{annotation.files.dir} (required, string) - Folder where
#'    annotation files are located.
#'    \item\code{output.dir} (required, string) - output folder for results
#'    \item\code{dataset.info} (required, string) - genome metadata file, it
#'    should contain at least:
#'    \itemize{
#'      \item File names. Please notice this information should be the first
#'      column in metadata file;
#'      \item Phenotype data (numeric, this is the value KOMODO2 uses to rank
#'      species when searching for associations)
#'      \item Normalization data (numeric, this is the value KOMODO2 uses as a
#'      denominator to compute annotation term frequencies to remove potential
#'      biases caused by, for instance, overannotation of model organisms or
#'      large differences in the counts of genomic elements). Please notice that
#'      KOMODO2 does not require normalization data for GO, as it computes the
#'      total number of GO terms per species and uses it as a normalizing factor.
#'    }
#'    \item\code{x.column} (required, numeric) - which column in "dataset.info"
#'    contains the phenotype data?
#'    \item\code{ontology} (required, string)  - which dictionary data type to
#'    use? Possible values are "GO" and "other". For GO, KOMODO2 can compute
#'    normalization data.
#'    \item\code{dict.path} (required, string)  - file for dictionary file
#'    (two-column file containing annotation IDs and their descriptions. Not
#'    needed for GO.
#'    \item\code{column} (required, string)  - which column in annotation files
#'    should be used (column name)
#'    \item\code{denominator.column} (optional, numeric) - which column contains
#'    normalization data (column number)
#'    \item\code{tree.path} (required, string)  - path for tree file in either
#'    newick or nexus format
#'    \item\code{tree.type} (required, string) - tree file type (either "nexus"
#'    or "newick")
#'    \item\code{cores} (optional, numeric) - how many cores to use? If not
#'    provided the function defaults to 1.
#'    \item\code{linear.model.cutoff} (required, numeric) - parameter that
#'    regulates how much graphical output is produced. We configure it to
#'    generate plots only for annotation terms with corrected q-values for
#'    phylogenetically independent contrasts (standard: smaller than 0.5).
#'    \item\code{MHT.method} (optional, string) - type of multiple hypothesis
#'    correction to be used. Accepts all methods listed by
#'    `stats::p.adjust.methods()`. If not provided the function defaults to
#'    "BH".
#' }
#'
#' @param defs either a KOMODO2-type list object or
#' a path to a text file containing the required definitions (see Details).
#' @param type type of analysis to perform. Currently only "correlation" is
#' supported.
#' @param cores positive integer, how many CPU cores to use (multicore
#' acceleration does not work in Windows systems). Setting
#' this parameter overrides any `cores` field from `defs`. Multicore support is
#' currently implemented using the `parallel` package, which uses forking
#' (which means that multicore support is not available under Windows)
#' @param render.report logical: should a HTML5 report be generated?
#'
#' @importFrom assertthat assert_that is.count has_name
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Install packages for report generation (only needs to be done once)
#' install_and_update_packages()
#'
#' # Download data files
#' retrieve_data_files(target.dir = "./data_folder")
#' defs <- "./data_folder/parameters_validation/parameters_gene2GO_Pan_proxy.txt"
#'
#' # Run KOMODO2
#' res <- run_KOMODO2(defs, cores = parallel::detectCores() - 1)
#' }

# saveRDS(defs, "./results/test_defs.rds")

run_KOMODO2 <- function(defs, type = "correlation",
                        cores = NULL, render.report = TRUE){

  # ================== Sanity checks ==================
  assert_that(is.list(defs) || file.exists(defs),
              is.null(cores) || is.count(cores),
              is.null(type) || is.character(type),
              is.null(type) || length(type) == 1,
              is.logical(render.report),
              length(render.report) == 1,
              msg = "input error(s) in KOMODO2::run_KOMODO2()")

  # If defs is a file path, read it into list
  if(!is.list(defs)) {
    defs <- read_komodo2_file(defs)
  }

  if(is.null(defs$type)) defs$type <- type

  # ================== Set up parallel processing ==================
  if(!is.null(cores)) {
    defs$cores <- cores
  } else if(!has_name(defs, "cores")) {
    defs$cores <- 1
  }

  available.cores <- parallel::detectCores()
  if (defs$cores >= available.cores){
    cat("\nAttention: cores too large, we only have ", available.cores,
        " cores.\nUsing ", available.cores - 1,
        " cores for load_data().")
    defs$cores <- available.cores - 1
  }
  if (.Platform$OS.type == "windows"){
    defs$cl <- parallel::makeCluster(defs$cores, setup_strategy = "sequential")
  }

  defs <- load_data(defs)      # Load required data
  defs <- clean_data(defs)     # Preliminary data cleaning
  defs <- do_analysis(defs)    # perform the analysis
  defs <- save_tsv_files(defs) # Save results to .tsv files
  defs <- make_report(defs, render.report = render.report) # generate HTML page

  if (.Platform$OS.type == "windows"){
    ## Stop the cluster
    parallel::stopCluster(defs$cl)
  }
  invisible(defs)
}
