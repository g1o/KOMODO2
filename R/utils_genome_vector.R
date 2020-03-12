# Small genome vector functions
# (not exported to the package namespace)


# Count number of genes in each genome of a group.
#
# Input:
#   someAnno: list of genomes, each with a data frame that maps
#               each gene to its annotations (GO, KO).
#   genome.names: char vector of names of the genomes to count elements.
#                   May be used to restrict which genomes to count in this
#                   function.
#   mode: char, defines whether KOMODO2 must consider all elements in all
#           genomes (default), or treat each element independently of
#           others (experiment). The latter is an experiment for
#           alignments.
# Returns an integer vector with the number of elements in each genome
#   of someAnno.
GroupElementCount <- function(someAnno, genome.names = NULL, mode = "default") {

  if (is.null(genome.names) | length(genome.names) == 0) {
    # TODO: should this return zero or throw an error?
    return(0)
  }

  if (mode == "default") {
    elementCount <- sapply(X   = someAnno[genome.names],
                           FUN = length)

  } else if (mode == "experiment") {
    elementCount        <- rep(1, length(someAnno[genome.names]))
    names(elementCount) <- names(someAnno[genome.names])

  } else {
    stop("'", mode, "' is not a recognized mode in KOMODO2::GroupElementCount()")
  }

  return(elementCount)
}


# Obtains the annotation from a file about a specific genome.
#
# Input:
#   file: char, path to the file with the annotation.
#   column: char label or integer index of column with the desired annotation.
# Returns list with the annotation of that genome.
FileAddition <- function(file, column) {

  # ===== Sanity check =====
  assertthat::assert_that(is.character(file),
                          length(file) == 1,
                          file.exists(file),
                          is.character(column) || is.numeric(column),
                          length(column) == 1)


  anno <- utils::read.table(file,
                            sep = "\t", header = TRUE,
                            colClasses = "character",
                            strip.white = TRUE, comment.char = "",
                            row.names = 1, check.names = FALSE)

  # ===== Sanity check =====
  if (is.numeric(column)){
    assertthat::assert_that(column %% floor(column) == 0,
                            column >= 1, column <= ncol(anno))
  } else {
    assertthat::assert_that(column %in% names(anno))
  }

  # Parsing the tab format from Uniprot
  genomeMap <- strsplit(anno[, column], split = " *; *")
  names(genomeAnno) <- rownames(anno)

  return(genomeAnno)
}


# Generates a phyletic vector for a specific genome.
# Input:
#   genomeAnno: list, annotation of the genome.
#   ontologyInfo: list, wrapper for the ontology. Since it may be GO, KO, or an
#                   arbitrary one, the information may come in different
#                   variables, hence the wrapper.
# Returns list with the frequency count of each term in that genome.
GenerateGenomeVector <- function(genomeAnno, ontologyInfo, column) {

  # ===== Sanity check =====
  assertthat::assert_that(is.list(genomeAnno) || is.character(genomeAnno),
                          is.list(ontologyInfo),
                          "name" %in% names(ontologyInfo))

  # Check input format (file or annotation list), adapt if possible
  if (is.character(genomeAnno)) {
    if (utils::file_test("-f", genomeAnno)) {
      genomeAnno <- FileAddition(genomeAnno, column)
    } else {
      # TODO: Should it throw an error in this case?
      return(NULL)
    }
  }
  if (!is.list(genomeAnno)) {
    # TODO: Should it throw an error in this case?
    return(NULL)
  }

  genomeVector <- rep.int(0, times = length(ontologyInfo$name))
  names(genomeVector) <- ontologyInfo$name

  if (is.element(ontologyInfo$ontology, c("go", "gene ontology"))) {
    genomeAnno <- lapply(genomeAnno, RemoveObsoleteAndAlternative,
                         ontologyInfo$allObsolete, ontologyInfo$allSynonym)
    genomeAnno <- lapply(genomeAnno, ObtainGeneGOancestors,
                         ontologyInfo$allAncestor)
  }

  # If the genome has any annotation, count terms; just return otherwise
  # 'countIDs' has the names of ontologic terms found in a genome and the
  # number of genes/proteins/elements in which they appeared
  # countIDs[, 1] = ontologic terms
  # countIDs[, 2] = occurrences of the term
  genomeAnno <- unlist(genomeAnno)
  if (length(genomeAnno) > 0) {
    countIDs <- as.data.frame(table(genomeAnno), stringsAsFactors = FALSE)
    countIDs <- countIDs[is.element(countIDs[, 1], ontologyInfo$name), ]

    genomeVector[countIDs[, 1]] <- genomeVector[countIDs[, 1]] +
      countIDs[, 2]
  }

  return(genomeVector)
}



