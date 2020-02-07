# clean.R - second step of the LCFD workflow of KOMODO2.
# Responsible for dealing with troubling aspects of the data,
#  such as missing values, merging, outliers and undesired characteres.
# Can also adapt data to allow for more flexible inputs form the user,
#  such as automatically converting common annotation output to a single
#  standard format. 
#
# input:
#   test.anno: (list) list of data frames, each one with the annotation
#                     of a genome of the test group.
#   back.anno: (list) list of data frames, each one with the annotation
#                     of a genome of the background group.
# output:
#   test.anno: (list) list of treated data frames, each one with the
#                     annotation of a genome of the test group.
#   back.anno: (list) list of treated data frames, each one with the
#                     annotation of a genome of the background group.


# Function to parse the tab format from Uniprot
KOMODO2$GenomeMap <- function(genome, column) {
  genomeMap <- strsplit(genome[, column], " *; *")
  names(genomeMap) <- rownames(genome)
  return(genomeMap)
}

if (KOMODO2$type == "significance") {

  KOMODO2$test.name <- basename(KOMODO2$test.name)
  KOMODO2$back.name <- basename(KOMODO2$back.name)
  
  names(KOMODO2$test) <- KOMODO2$test.name
  names(KOMODO2$back) <- KOMODO2$back.name
  
  # Safety check, removes the counts if any value is missing
  if (!is.null(KOMODO2$testElementCount)) {
    names(KOMODO2$testElementCount) <- KOMODO2$test.name
    if (isTRUE(any(is.na(KOMODO2$testElementCount)))) {
      KOMODO2$testElementCount <- NULL
    }
  }
  if (!is.null(KOMODO2$backElementCount)) {
    names(KOMODO2$backElementCount) <- KOMODO2$back.name
    if (isTRUE(any(is.na(KOMODO2$backElementCount)))) {
      KOMODO2$backElementCount <- NULL
    }
  }
  
  KOMODO2$test.anno <- lapply(KOMODO2$test, KOMODO2$GenomeMap, KOMODO2$column)
  KOMODO2$back.anno <- lapply(KOMODO2$back, KOMODO2$GenomeMap, KOMODO2$column)
  
  KOMODO2$GenomeMap <- NULL
  KOMODO2$test <- NULL
  KOMODO2$back <- NULL
  
} else if (KOMODO2$type == "correlation") {
  
  KOMODO2$y.name <- basename(KOMODO2$y.name)
  names(KOMODO2$y) <- KOMODO2$y.name
  names(KOMODO2$x) <- KOMODO2$y.name

  if (length(KOMODO2$denominator) <= 1) {
    KOMODO2$denominator <- NULL
  } else {
    names(KOMODO2$denominator) <- KOMODO2$y.name
  }

  # Safety check, returns error if any value is missing
  if (isTRUE(any(is.na(KOMODO2$x)))) {
    stop("Values missing for the (x variable), check your dataset info file.")
  }
  
  KOMODO2$x <- as.data.frame(KOMODO2$x)
  
  KOMODO2$y.anno <- lapply(KOMODO2$y, KOMODO2$GenomeMap, KOMODO2$column)
  KOMODO2$GenomeMap <- NULL
  KOMODO2$y <- NULL
  

}

if (KOMODO2$ontology == "other" & KOMODO2$dict.path != "") {
  KOMODO2$dictionary <- unique(KOMODO2$dictionary)

  # Conversion to a named list
  KOMODO2$temp <- as.list(KOMODO2$dictionary[, 2])
  names(KOMODO2$temp) <- KOMODO2$dictionary[, 1]
  KOMODO2$dictionary <- KOMODO2$temp
  KOMODO2$temp <- NULL
}


