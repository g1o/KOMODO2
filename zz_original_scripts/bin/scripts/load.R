# load.R - first step of the LCFD workflow of KOMODO2.
# Responsible for the loading of all the data required.
# Separates the data loading, which can be the longest step of a workflow,
#  from the analysis itself, which can be redone multiple times.
# 
# input (KOMODO2 list):
#   type: (char) which comparison module to use:
#                - "significance": compares two groups of genomes within an 
#                                  ontology
#                - "correlation": establishes how much a variable explains the 
#                                 variations seen in the genomes
#
#   test.path: (char) path of dir with annotation files of the test group.
#   back.path: (char) path of dir with annotation files of the background one.
#
#   x.path: (char) file with genomes' attributes for correlation test.
#   y.path: (char) directory with the genomes and their annotations 
#                  (correlation test).
#
#   ontology: (char) which ontology to use: "GO" (or "Gene Ontology"),
#                    "KEGG" and "other"
# 
#   dict.path: (char) file with the dictionary (terms and their meaning) 
#                     of an ontology, used if ontology is set as "other"
#
# output:
#   test.anno: (list) list of data.frames with annotation of genomes
#                     (test group)
#   back.anno: (list) same, but for the genomes of the background group.
#
#   attributes: (data.frame) table with the numeric attributes of the genomes,
#                           used in correlation tests.
#   annotations: (list) list of data.frames with annotation for each genome,
#                       used in correlation tests.
#
#   dictionary: (data.frame) table with the terms of an ontology and their
#                            meaning


# Load data files

KOMODO2$type <- "correlation"

if (KOMODO2$type == "significance") {

  # If the paths point to a file with the paths of the genome files
  # May have the total number of elements in that genome in the second column
  # in case the files themselves have missing (with no annotated) elements
  if (isTRUE(file_test("-f", KOMODO2$test.path))) {
    KOMODO2$test.name <- read.table(KOMODO2$test.path, sep = "\t",
                             strip.white = TRUE, comment.char = "",
                             check.names = FALSE, header = FALSE)
    if (ncol(KOMODO2$test.name) == 2) {
      KOMODO2$testElementCount <- KOMODO2$test.name[, 2]
    }
    KOMODO2$test.name <- as.character(KOMODO2$test.name[, 1])
  }
  if (isTRUE(file_test("-f", KOMODO2$back.path))) {
    KOMODO2$back.name <- read.table(KOMODO2$back.path, sep = "\t",
                             strip.white = TRUE, comment.char = "",
                             check.names = FALSE, header = FALSE)
    if (ncol(KOMODO2$back.name) == 2) {
      KOMODO2$backElementCount <- KOMODO2$back.name[, 2]
    }
    KOMODO2$back.name <- as.character(KOMODO2$back.name[, 1])
  }
  
  
  # If the paths point to directories, get the relative path of its files
  if (isTRUE(file_test("-d", KOMODO2$test.path))) {
    KOMODO2$test.name <- as.character(list.files(path = KOMODO2$test.path,
                                      all.files = FALSE, full.names = TRUE,
                                      recursive = FALSE))
  }
  if (isTRUE(file_test("-d", KOMODO2$back.path))) {
    KOMODO2$back.name <- as.character(list.files(path = KOMODO2$back.path,
                                      all.files = FALSE, full.names = TRUE,
                                      recursive = FALSE))
  }
  
  # Content of files
  KOMODO2$test <- lapply(KOMODO2$test.name, read.table, sep = "\t",
                         header = TRUE, colClasses = "character",
                         strip.white = TRUE, comment.char = "",
                         row.names = 1, check.names = FALSE)
  KOMODO2$back <- lapply(KOMODO2$back.name, read.table, sep = "\t",
                         header = TRUE, colClasses = "character",
                         strip.white = TRUE, comment.char = "",
                         row.names = 1, check.names = FALSE)
                         


} else if (KOMODO2$type == "correlation") {


#   KOMODO2$x <- read.table(KOMODO2$x.path, sep = "\t", 
#                           header = TRUE, colClasses = "character",
#                           strip.white = TRUE, comment.char = "")

#   KOMODO2$y <- read.table(KOMODO2$y.path, sep = "\t",
#                           header = TRUE, colClasses = "character",
#                           strip.white = TRUE, comment.char = "")
  
  
  KOMODO2$y.name <- read.csv(KOMODO2$dataset.info, header = FALSE,
                             strip.white = TRUE, comment.char = "",
                             check.names = FALSE, sep = "\t")
  KOMODO2$x <- KOMODO2$y.name[, KOMODO2$x.column]
  if (KOMODO2$denominator.column == "") {
    KOMODO2$denominator <- ""
  } else {
    KOMODO2$denominator <- KOMODO2$y.name[, KOMODO2$denominator.column]
  }
  
  KOMODO2$y.name <- as.character(KOMODO2$y.name[, 1])
  KOMODO2$y.name <- paste0(KOMODO2$annotation_files_dir, KOMODO2$y.name)
  KOMODO2$y <- lapply(KOMODO2$y.name, read.csv, sep = "\t", header = TRUE, 
                      colClasses = "character", strip.white = TRUE, 
                      comment.char = "", check.names = FALSE)
  
}

# Load ontology's dictionary, if it isn't currently supported by KOMODO2
if (KOMODO2$ontology == "other" & KOMODO2$dict.path != "") {

  KOMODO2$dictionary <- read.csv(KOMODO2$dict.path, sep = "\t", quote = "", 
                                 colClasses = "character", 
                                 strip.white = TRUE, comment.char = "#")

}

