#' Retrieve data files from the Github repository
#'
#' This script retrieves relevant data files from the KOMODO2 project
#' repository. It will download the data into a folder containing
#' directories related to dictionary files, Gene Ontology annotation
#' files, tree files, etc.
#'
#' If the `target.dir` provided does not exist it is created
#' (recursively) by the function.
#'
#' @param target.dir path to the folder where the files will be saved (
#' accepts relative and absolute paths)
#' @param method Method to be used for downloading files. Current download
#' methods are "internal", "wininet" (Windows only) "libcurl", "wget" and
#' "curl", and there is a value "auto": see ‘Details’ and ‘Note’ in the
#' documentation of \code{utils::download.file()}.
#' @param unzip The unzip method to be used. See the documentation of
#' \code{utils::unzip()} for details.
#' @param url repository URL. Do not change unless you really know what
#' you're doing.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   retrieve_data_files(target.dir = "./data_folder")
#' }

retrieve_data_files <- function(target.dir,
                                method = "auto",
                                unzip  = getOption("unzip"),
                                url = "https://github.com/fcampelo/KOMODO2-CRAN/"){

  assertthat::assert_that(is.character(target.dir),
                          length(target.dir) == 1,
                          is.character(url),
                          length(url) == 1)

  # Remove trailing slash
  target.dir <- gsub("/$", "", target.dir)
  url <- gsub("/$", "", url)

  if(!dir.exists(target.dir)){
    dir.create(target.dir, recursive = TRUE)
  } else {
    filelist <- dir(target.dir, full.names = TRUE)
    unlink(filelist, recursive = TRUE, force = TRUE)
  }

  url_full <- gsub("//", "/", paste0(url, "/archive/master.zip"),
                   fixed = TRUE)
  cat("\nRetrieving online resources... ")
  res1 <- utils::download.file(url_full,
                               quiet    = TRUE,
                               destfile = paste0(target.dir, "/tmpdata.zip"),
                               cacheOK  = FALSE,
                               method   = method)
  if(res1 != 0) stop("Error downloading file \n", url_full)

  utils::unzip(paste0(target.dir, "/tmpdata.zip"),
               unzip = unzip,
               exdir = target.dir)

  file.remove(paste0(target.dir, "/tmpdata.zip"))

  srcdir <- dir(target.dir, full.names = TRUE)
  tmpdir <- paste0(srcdir, "/data_files/")

  res2 <- all(file.copy(from      = paste0(tmpdir, dir(tmpdir)),
                        to        = target.dir,
                        overwrite = TRUE,
                        recursive = TRUE))

  if (!res2) stop("Error processing downloaded file.")

  unlink(srcdir, recursive = TRUE)

  cat("done!\n")

  invisible(TRUE)
}
