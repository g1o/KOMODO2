# func.R - part of the last step of the LCFD workflow of KOMODO2.
# Contains any function needed to perform the actual analysis.
# Sourced by the file "do.R", reloading a modified version of this file doesn't
#  require calling steps 1 and 2 again.
#
# input:
#   none.
# output:
#   none.

GenerateTree <- function(taxonIds) {
  taxize_class <- classification(taxonIds, db = "ncbi")
  taxize_tree <- class2tree(taxize_class, check = TRUE)
  plot(taxize_tree)
#  return(tree)
}
  


InstallPackages <- function() {
  # Install the latest version of Bioconductor and all necessary libraries.
  # 
  # Args:
  #   none.
  # Returns:
  #   nothing, just sets everything for use.
  #
  # Additional notes:
  #   optional function, use it if LoadBioconductor() fails to run.

  source("http://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("GO.db")
  biocLite("KEGGREST")
  install.packages("boot")
  install.packages("Cairo")
  install.packages("vcd")
  install.packages("gplots")
  install.packages("ggplot2")
  install.packages("foreach")
  install.packages("doParallel")
  install.packages("iterators")
  install.packages("itertools")
  install.packages("CHNOSZ")
  install.packages("taxize")
  install.packages("plotly") #for visualization of results
  install.packages("ggplot2")#for visualization of results
  install.packages("scales") #for visualization of results
  install.packages("DT")     #for visualization of results
	
}


LoadBioconductor <- function() {
  # Loads Bioconductor and all KOMODO2's dependencies.
  # 
  # Args:
  #   none.
  # Returns:
  #   nothing, just sets Bioconductor for use.
  #
  # Additional notes:
  #   used at the start of the program, required for ObtainGeneGOancestors()
  #   and similar functions.


  suppressMessages(source("http://bioconductor.org/biocLite.R"))
  suppressMessages(library("GO.db"))
  suppressMessages(library("taxize"))
  suppressMessages(library("ape"))
  suppressMessages(library("nlme"))
  suppressMessages(library("KEGGREST"))
  suppressMessages(library("boot"))
  suppressMessages(library("Cairo"))
  suppressMessages(library("vcd"))
  suppressMessages(library("gplots"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("foreach"))
  suppressMessages(library("doParallel"))
  suppressMessages(library("iterators"))
  suppressMessages(library("itertools"))
  suppressMessages(library("CHNOSZ"))
  suppressMessages(library("plotly"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("scales"))
  suppressMessages(library("DT"))

}


ReloadFunc <- function() {
  # Function for debugging, updates the funcions in func.R mid-analysis
  #  without erasing any data carried so far. 
  #
  # Args:
  #   none.
  # Returns:
  #   none.
  
  KOMODO2.env <- new.env()
  suppressMessages(sys.source("func.R", envir = KOMODO2.env))
  suppressMessages(attach(KOMODO2.env))
  rm(KOMODO2.env)
}



# -=-=-=- Ontology manipulation -=-=-=-

ListAncestors <- function() {
  # Prepares a listing of ancestors for each GO ID. Used to vectorize later 
  # operations.
  # 
  # Args:
  #   none.
  # Returns:
  #   allAncestor: (list) All GOXXANCESTOR combined.
  #                       Used to vectorize ObtainGeneGOancestors().

  allAncestor <- c(as.list(GOBPANCESTOR), as.list(GOCCANCESTOR),
                   as.list(GOMFANCESTOR))
  
  return(allAncestor)
}


ListSynonyms <- function() {
  # Prepares a listing of synonymous for GO IDs with alternative ids. Used to
  # vectorize later operations.
  # 
  # Args:
  #   none.
  # Returns:
  #   allSynonym: (char) GOSYNONYM's mapping in a single vector.
  #                      Used to vectorize RemoveObsoleteAndAlternative().

  allSynonym <- as.character(GOSYNONYM)
  return(allSynonym)
}


ListObsoletes <- function() {
  # Prepares a listing of obsolete GO IDs. Used to vectorize later operations.
  # 
  # Args:
  #   none.
  # Returns:
  #   allObsolete: (vector) vector with GOOBSOLETE's mapping.
  #                         Used to vectorize RemoveObsoleteAndAlternative().

  allObsolete <- as.character(GOOBSOLETE)
  names(allObsolete) <- NULL
  return(allObsolete)
}


ListKOs <- function() {
  # Prepares a listing of KO and their annotation.
  # 
  # Args:
  #   none.
  # Returns:
  #   allKOs: (vector) char vector of KOs and their annotation.

  allKOs <- keggList("ko")
  names(allKOs) <- substr(names(allKOs), 4, 9)  # remove "ko:" part from names
  return(allKOs)
}


CreateDictionary <- function(test.anno, back.anno = NULL) {
  # Creates an dictionary of valid terms from the data; does not contain an
  #  annotation description.
  # 
  # Args:
  #   test.anno: (list) genomes in the test group, each with a data frame that
  #                     maps each genomic element to its annotations. 
  #   back.anno: (list) genomes in the background group, each with a data frame
  #                     that maps each genomic element to its annotations.
  # Returns:
  #   dict: (char) dictionary of terms.

  # Parsing the tab format from Uniprot
  dict <- unique(unlist(sapply(test.anno, unlist)))
  if (!is.null(back.anno)){
    dict <- unique(c(dict, unlist(sapply(back.anno, unlist)))) 
  }
  names(dict) <- dict

  return (dict)
}


RemoveObsoleteAndAlternative <- function(geneIDs, allObsolete, allSynonym) {
  # Removes GOs that are obsolete and replace alternative ids for main ids.
  # 
  # Args:
  #   geneIDs: (vector) vector of char factors mapping gene IDs to GO ID's.
  #   allObsolete: (vector) vector with GOOBSOLETE's mapping.
  #                         Used to vectorize RemoveObsoleteAndAlternative().
  #   allSynonym: (char) GOSYNONYM's mapping in a single vector.
  #                      Used to vectorize RemoveObsoleteAndAlternative().
  # Returns:
  #   geneIDs: (vector) vector of char factors mapping gene IDs to GO ID's
  #                     without obsolete and alternative GO IDs.

  geneIDs <- setdiff(geneIDs, allObsolete)

  alternativeIDs <- intersect(geneIDs, names(allSynonym))
  geneIDs <- setdiff(geneIDs, alternativeIDs)
  newIDs <- allSynonym[alternativeIDs]
  geneIDs <- unique(c(geneIDs, newIDs))

  return(geneIDs)
}



ObtainGeneGOancestors <- function(geneIDs, allAncestor) {
  # Finds all GO ID ancestors for a given gene
  # 
  # Args:
  #   geneIDs: (vector) vector of char factors mapping gene IDs to GO ID's.
  #   allAncestor: (list) All GOXXANCESTOR combined.
  #                       Used to vectorize ObtainGeneGOancestors().
  # Returns:
  #   geneIDs: (char) vector with all GO IDs found for the gene.
  
  geneAncestors <- unlist(allAncestor[geneIDs], use.names = FALSE)
  geneAncestors <- geneAncestors[!is.null(geneAncestors)]
  geneAncestors <- geneAncestors[!geneAncestors == "all"]
  geneIDs <- unique(c(geneIDs, geneAncestors))

  geneIDs <- as.character(geneIDs)
  return(geneIDs)
}


# -=-=-=- Genome Vector internal functions -=-=-=-


GroupElementCount <- function(someAnno, genome.names=NULL, mode = "default") {
  # Returns a vector with the number of genes in each genome of a group.
  # 
  # Args:
  #   someAnno: (list) list of genomes, each with a data frame that maps
  #                    each gene to its annotations (GO, KO). 
  #   genome.names: (vector) names of the genomes to count elements. May
  #                          restrict which genomes to count in this function.
  #   mode: (char) defines whether KOMODO2 must consider all elements in all
  #                genomes (default), or treat each element independently of
  #                others (experiment). The latter is an experiment for 
  #                alignments.
  # Returns:
  #   elementCount: (integer) vector with the number of elements in each genome
  #                           of someAnno.
  
  if (is.null(genome.names) | length(genome.names) == 0) {
    return(0)
  }
  
  if (mode == "default") {
    elementCount <- sapply(someAnno[genome.names], length)
  } else if (mode == "experiment") {
    elementCount <- rep(1, length(someAnno[genome.names]))
    names(elementCount) <- names(someAnno[genome.names])
  }
  
  return(elementCount)
}


FileAddition <- function(file, column) {
  # Obtains the annotation from a file about a specific genome.
  # 
  # Args:
  #   file: (char) path to the file with the annotation.
  #   column: (char) label (or index) of column with the desired annotation.
  # Returns:
  #   genomeAnno: (list) list with the annotation of that genome.

  anno <- read.table(file, sep = "\t", header = TRUE, colClasses = "character",
                   strip.white = TRUE, comment.char = "", row.names = 1,
                   check.names = FALSE)
  
  # Parsing the tab format from Uniprot
  genomeMap <- strsplit(anno[, column], " *; *")
  names(genomeAnno) <- rownames(anno)

  return(genomeAnno)
}



GenerateGenomeVector <- function(genomeAnno, ontologyInfo) {
  # Generates a phyletic vector for a specific genome.
  # 
  # Args:
  #   genomeAnno: (list) annotation of the genome.
  #   ontologyInfo: (list) wrapper the ontology. Since it may be GO, KO, or an
  #                        arbitrary one, the information may come in different
  #                        variables, hence the wrapper.
  # Returns:
  #   genomeVector: (list) list with the frequency count of each term
  #                        in that genome.
  
  # Check input format (file or annotation list), adapt if possible
  if (is.character(genomeAnno)) {
    if (file_test("-f", genomeAnno)) {
      genomeAnno <- FileAddition(file, KOMODO2$column)
    }
  }
  if (!is.list(genomeAnno)) {
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
    countIDs <- as.data.frame(table(genomeAnno))
    countIDs <- countIDs[is.element(countIDs[, 1], ontologyInfo$name), ]

    # countIDs store terms as factor, convert to character with rownames()
    rownames(countIDs) <- countIDs[, 1]
    genomeVector[rownames(countIDs)] <- genomeVector[rownames(countIDs)] +
                                          countIDs[, 2]
  }

  return(genomeVector)
}


# -=-=-=- Genome Vector operations -=-=-=-


AddGenomeVectors <- function(someAnno, genome.names, someGV = NULL) {
  # Generates phyletic vectors of each genome in a group, which stores the
  # count of appearences of a GO term in that genome.
  # 
  # Args:
  #   someAnno: (list) annotation from the original input data.
  #   genome.names: (char) a vector with the names of the genomes to select.
  #   someGV: (data.frame) a previously processed set of genome vectors.
  # Returns:
  #   genomeVectors: (list) list with the phyletic vector of each genome of
  #                         the input group.

  # To generate genome vectors, we need some information about the ontology.
  # Since it may be GO, KO or an arbitrary one, the information may come 
  # from different variables, hence the following wrapper.
  ontologyInfo <- list(ontology = tolower(KOMODO2$ontology))
  if (is.element(ontologyInfo$ontology, c("go", "gene ontology"))) { 
    ontologyInfo$allAncestor <- KOMODO2$allAncestor
    ontologyInfo$allObsolete <- KOMODO2$allObsolete
    ontologyInfo$allSynonym <- KOMODO2$allSynonym
    ontologyInfo$name <- names(KOMODO2$allAncestor)
  } else if (ontologyInfo$ontology == "kegg") {
    ontologyInfo$name <- names(KOMODO2$allKOs)
  } else if (ontologyInfo$ontology == "other") {
    ontologyInfo$name <- names(KOMODO2$dictionary)
  }  
  
  # Select genomes from someAnno that are in genome.names and not in someGV,
  # so we don't add genome vectors that already exist.
  genome.names <- setdiff(genome.names, rownames(someGV))
  someAnno <- someAnno[genome.names]
  
  # Add the genome vectors and fit them into a data.frame format.
  genomeVectors <- foreach(genomeAnno = iter(someAnno), .combine = "rbind",
                           .export = "GenerateGenomeVector") %dopar% {
    GenerateGenomeVector(genomeAnno, ontologyInfo)
  }
  rownames(genomeVectors) <- names(someAnno)
  genomeVectors <- as.data.frame(genomeVectors, optional = TRUE)
  
  # Merge with the existing genome vectors, if given.
  if (!is.null(someGV)) {
    genomeVectors <- rbind(someGV, genomeVectors)
  }
  
  return(genomeVectors)
}


SelectGenomeVectors <- function(someGV, genome.names) {
  # Select genomes from someAnno that are in genome.names. We don't do it 
  # directly because some genomes of genome.names may be missing in someGV.
  # 
  # Args:
  #   someGV: (data.frame) a previously processed set of genome vectors.
  #   genome.names: (char) a vector with the names of the genomes to select.
  # Returns:
  #   someGV: (data.frame) a previously processed set of genome vectors, with
  #                        only the existing genomes in genome.names.
  
  return(someGV[intersect(genome.names, rownames(someGV)), ])
}



UpdateGenomeVector <- function(oldGenomeAnno, newGenomeAnno, genomeVector,
                               ontologyInfo) {
  # Generates a difference phyletic vector for a genome annotation update.
  # 
  # Args:
  #   oldGenomeAnno: (list) genome annotation from a previous run.
  #   newGenomeAnno: (list) genome annotation to be analyzed.
  #   genomeVector: (data.frame) genome vector of an input genome, with the
  #                              frequency of each term.
  #   ontologyInfo: (list) wrapper the ontology. Since it may be GO, KO, or an
  #                        arbitrary one, the information may come in different
  #                        variables, hence the wrapper.
  # Returns:
  #   genomeVector: (data.frame) updated genome vector.

  if (is.character(newGenomeAnno)) {
    if (file_test("-f", newGenomeAnno)) {
      newGenomeAnno <- FileAddition(file, KOMODO2$column)
    }
  }
  if (!is.list(newGenomeAnno) || !is.list(oldGenomeAnno)) {
    return(NULL)
  }
  
  # For common elements
  inCommon <- intersect(names(oldGenomeAnno), names(newGenomeAnno))
  isEqual <- Map("setequal", oldGenomeAnno[inCommon],
                             newGenomeAnno[inCommon])
  differentAnno <- names(isEqual[isEqual == FALSE])
  for (i in differentAnno) {
    remove <- oldGenomeAnno[[i]]
    add <- newGenomeAnno[[i]]
    if (is.element(ontologyInfo$ontology, c("go", "gene ontology"))) {
      remove <- RemoveObsoleteAndAlternative(remove, 
                                             ontologyInfo$allObsolete,
                                             ontologyInfo$allSynonym)
      remove <- ObtainGeneGOancestors(remove, ontologyInfo$allAncestor)
      add <- RemoveObsoleteAndAlternative(add, 
                                          ontologyInfo$allObsolete,
                                          ontologyInfo$allSynonym)
      add <- ObtainGeneGOancestors(add, ontologyInfo$allAncestor)
    }
    genomeVector[remove] <- genomeVector[remove] - 1
    genomeVector[add] <- genomeVector[add] + 1
  }

  # For removed elements
  removedAnno <- setdiff(names(oldGenomeAnno), names(newGenomeAnno))
  for (i in removedAnno) {
    remove <- oldGenomeAnno[[i]]
    if (is.element(ontologyInfo$ontology, c("go", "gene ontology"))) {
      remove <- RemoveObsoleteAndAlternative(remove, 
                                             ontologyInfo$allObsolete,
                                             ontologyInfo$allSynonym)
      remove <- ObtainGeneGOancestors(remove, ontologyInfo$allAncestor)
    }
    genomeVector[remove] <- genomeVector[remove] - 1
  }

  # For new elements
  introducedAnno <- setdiff(names(newGenomeAnno), names(oldGenomeAnno))
  for (i in introducedAnno) {
    add <- newGenomeAnno[[i]]
    if (is.element(ontologyInfo$ontology, c("go", "gene ontology"))) {
      add <- RemoveObsoleteAndAlternative(add, 
                                          ontologyInfo$allObsolete,
                                          ontologyInfo$allSynonym)
      add <- ObtainGeneGOancestors(add, ontologyInfo$allAncestor)
    }
    genomeVector[add] <- genomeVector[add] + 1
  }
  
  return(genomeVector)
}


# -=-=-=- Parameter Vector for parametric tests (Fisher currently) -=-=-=-



ParameterVectors <- function(testGV, backGV, 
                             testElementCount, backElementCount) {
  # Creates a summary of both test and background groups about the number of
  # annotated genes for each GO term.
  #
  # Args: 
  #   testGV (list) list of the number of occurrences of each term among the
  #                 genomes in the test group.
  #   backGV (list) list of the number of occurrences of each term among the
  #                 genomes in the background group.
  #   testElementCount: (integer) vector with the number of elements in each
  #                               genome of testGV.
  #   backElementCount: (integer) vector with the number of elements in each
  #                               genome of backGV.
  #
  # Returns:
  #   parameterVectors: (list) list with the phyletic vector of each genome of
  #                          the input group.


  t1 <- colSums(testGV)
  b1 <- colSums(backGV)

  nonZero <- t1 > 0 | b1 > 0
  t1 <- t1[nonZero, drop=FALSE]
  b1 <- b1[nonZero, drop=FALSE]

  testGroupNumber <- sum(testElementCount)
  backGroupNumber <- sum(backElementCount)

  t0 <- testGroupNumber - t1
  b0 <- backGroupNumber - b1

  parameterVectors <- Map("c", t1, b1, t0, b0)
  return(parameterVectors)
}


# -=-=-=- Bootstrap Setup -=-=-=-


BootGenomeVectorsIndices <- function(testGV, backGV, R = 100) {
  # Bootstraps the genome listings and return a list with the indices of the
  # bootstrap.
  #
  # Args: 
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  #   R: (numeric) number of samples to run in the bootstrap.
  #   
  # Returns:
  #   bootIndices: (list) list with the indices on the genome listings used by
  #                       the bootstrap.

  require(boot)

  bootIndicesFunc <- function(data, indices) {
    return(indices)
  }

  columns <- colSums(testGV) > 0 | colSums(backGV) > 0

  bootTestGroup <- boot(testGV[, columns], statistic = bootIndicesFunc, R = R)
  bootBackGroup <- boot(backGV[, columns], statistic = bootIndicesFunc, R = R)

  bootIndices <- list(test = bootTestGroup$t, back = bootBackGroup$t,
                      columns = names(columns[columns == TRUE]))
  return(bootIndices)
}


BootParameterVectors <- function(testGV, backGV, testElementCount,
                                 backElementCount, bootIndices) {
  # Creates a summary of both test and background groups about the number of
  # genes that are annotated for each GO term.
  #
  # Args: 
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  #   testElementCount: (integer) vector with the number of elements in each
  #                               genome of testGV.
  #   backElementCount: (integer) vector with the number of elements in each
  #                               genome of backGV.
  #   bootIndices: (list) list with the indices on the genome listings used by
  #                       the bootstrap.
  #
  # Returns:
  #   bParameterVectors: (list) list with the phyletic vector of each genome of
  #                             the input group.

  # This function is quicker with matrices than with data frames
  # Not the case of future functions, so we're making it local
  testGV <- as.matrix(testGV)
  backGV <- as.matrix(backGV)
  
  # Accessing data from the bootstrap with the indices
  t1 <- foreach(sample = iter(bootIndices$test, by = "row"), 
                .combine = "rbind") %dopar% {
    colSums(testGV[sample, bootIndices$columns, drop=FALSE])
  }
  b1 <- foreach(sample = iter(bootIndices$back, by = "row"),
                .combine = "rbind") %dopar% {
    colSums(backGV[sample, bootIndices$columns, drop=FALSE])
  } 
   
  t1 <- t1[, (colSums(t1) > 0 | colSums(b1) > 0), drop=FALSE]
  b1 <- b1[, colnames(t1), drop=FALSE]

  bParameterVectors <- array(0, c(nrow(t1), ncol(t1), 4))
  dimnames(bParameterVectors) <- list(row = 1:nrow(t1),
                                    col = colnames(t1),
                                    parameter = c("t1", "b1", "t0", "b0"))
  bParameterVectors[, , "t1"] <- t1
  bParameterVectors[, , "b1"] <- b1
  
  bParameterVectors[, , "t0"] <- sum(testElementCount) - t1
  bParameterVectors[, , "b0"] <- sum(backElementCount) - b1
  
  return(bParameterVectors)
}



# -=-=-=- Simple group comparison -=-=-=-

StatisticalTest <- function(parameterVectors, test = "fisher.over") {
  # Computes the chosen statistical test between two phyletic vectors.
  #
  # Args:
  #   parameterVectors: (list) the phyletic vector with the number of 
  #                            occurrences of each term among the genomes
  #                            in the groups. Format: (t0, t1, b0, b1).
  #   test: (char) the statistical test to perform.
  # Returns:
  #   results: (numeric) p-value of all terms.


  results <- vector(mode = "numeric", length = length(parameterVectors))
  names(results) <- names(parameterVectors)
  
  chunkSize <- ceiling(length(parameterVectors)/getDoParWorkers())
  
  if (test == "chi2") {
    # Separating data in chunks to improve memory usage when allocating
    # and managing multiple cores
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        parameters <- matrix(parameterVectors[[term]], nrow = 2)
        results[[term]] <- chisq.test(parameters)$p.value
      }
      results[chunk]
    }
  } else if (test == "fisher.over") {
    # Separating data in chunks to improve memory usage when allocating
    # and managing multiple cores
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        parameters <- matrix(parameterVectors[[term]], nrow = 2)
        results[[term]] <- fisher.test(parameters,
                                       alternative = "greater")$p.value
      }
      results[chunk]
    }
  } else if (test == "fisher.under") {
    # Separating data in chunks to improve memory usage when allocating
    # and managing multiple cores
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        parameters <- matrix(parameterVectors[[term]], nrow = 2)
        results[[term]] <- fisher.test(parameters, 
                                       alternative = "less")$p.value
      }
      results[chunk]
    }
  }
  
  results <- sort(results)
  return(results)
}



BootStatisticalTest <- function(bParameterVectors, test = "fisher.over") {
  # Computes the chosen statistical test between two phyletic vectors.
  #
  # Args:
  #   bParameterVectors: (list) the phyletic vector with the number of 
  #                             occurrences of each term among the genomes
  #                             in the groups. Format: (t1, b1, t0, b0).
  #   test: (char) the statistical test to perform.
  # Returns:
  #   bresults: (list) p-value of all terms.

  R <- nrow(bParameterVectors)
  chunkSize <- ceiling(R/getDoParWorkers())

  # In each sample, ensure that cases where both t1 = 0 and b1 = 0 won't be
  # evaluated.
  bresults <- vector(mode = "list", length = R)
  names(bresults) <- c(1:R)
  nonZero <- as.matrix((bParameterVectors[, , "t1"] + 
                        bParameterVectors[, , "b1"]) > 0)
  for (sample in 1:R) {
    bresults[[sample]] <- bParameterVectors[sample, nonZero[sample, ], "t1",
                                            drop=FALSE]
  }

  if (test == "chi2") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"),
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        for (term in colnames(bresults[[sample]])) {
          parameters <- matrix(bParameterVectors[sample, term, ], nrow = 2)
          bresults[[sample]][[term]] <- chisq.test(parameters)$p.value
        }
      }
      bresults[chunk]
    }
  } else if (test == "fisher.over") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"),
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        for (term in colnames(bresults[[sample]])) {
          parameters <- matrix(bParameterVectors[sample, term, ], nrow = 2)
          bresults[[sample]][[term]] <- fisher.test(parameters,
                                             alternative = "greater")$p.value
        }
      }
      bresults[chunk]
    }
  } else if (test == "fisher.under") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"),
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        for (term in colnames(bresults[[sample]])) {
          parameters <- matrix(bParameterVectors[sample, term, ], nrow = 2)
          bresults[[sample]][[term]] <- fisher.test(parameters,
                                             alternative = "less")$p.value
        }
      }
      bresults[chunk]
    }
  }

  for (sample in 1:R) {
    bresults[[sample]] <- sort(bresults[[sample]])
  }
  return(bresults)
}



# -=-=-=- Distribution comparison -=-=-=-

CompareDistributions <- function(parameterVectors, testGV, backGV,
                                 test = "ks.over") {
  # Computes the chosen statistical test between two phyletic vectors.
  #
  # Args:
  #   parameterVectors: (list) list with the genome vector of each genome of
  #                            the input group.
  #   testGV: (list) number of occurrences of each term among the genomes in
  #                  the test group.
  #   backGV: (list) number of occurrences of each term among the genomes in
  #                  the background group.
  #   test: (char) the statistical test to perform.
  # Returns:
  #   results: (list) p-value of all terms.


  results <- vector(mode = "numeric", length = length(parameterVectors))
  names(results) <- names(parameterVectors)

  chunkSize <- ceiling(length(parameterVectors)/getDoParWorkers())
  
  if (test == "ks.over") {
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        test <- testGV[[term]]
        back <- backGV[[term]]
        results[[term]] <- ks.test(test, back, alternative = "less")$p.value
      }
      results[chunk]
    }
  } else if (test == "ks.under") {
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        test <- testGV[[term]]
        back <- backGV[[term]]
        results[[term]] <- ks.test(test, back, alternative = "greater")$p.value
      }
      results[chunk]
    }
  } else if (test == "wilcox.over") {
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        test <- testGV[[term]]
        back <- backGV[[term]]
        results[[term]] <- wilcox.test(test, back,
                                       alternative = "greater")$p.value
      }
      results[chunk]
    }
  } else if (test == "wilcox.under") {
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        test <- testGV[[term]]
        back <- backGV[[term]]
        results[[term]] <- wilcox.test(test, back, 
                                       alternative = "less")$p.value
      }
      results[chunk]
    }
  } else if (test == "kruskal") {
    results <- foreach(chunk = ichunk(names(parameterVectors), chunkSize,
                       mode = "character"), .combine = "c") %dopar% {
      for (term in chunk) {
        test <- testGV[[term]]
        back <- backGV[[term]]
        results[[term]] <- kruskal.test(list(test, back))$p.value
      }
      results[chunk]
    }
  }

  results <- sort(results)
  return(results)
}


BootDistributions <- function(bParameterVectors, bootIndices, testGV, backGV,
                              test = "ks.over") {
  # Computes the chosen statistical test between two phyletic vectors.
  #
  # Args:
  #   bParameterVectors: (list) list with the genome vector of each genome of
  #                             the input group.
  #   bootIndices: (list) list with the indices on the genome listings used by
  #                       the bootstrap.
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  #   test: (char) the statistical test to perform.
  # Returns:
  #   results: (list) p-value of all terms.

  R <- nrow(bootIndices$test)
  chunkSize <- ceiling(R/getDoParWorkers())

  # In each sample, ensure that cases where both t1 = 0 and b1 = 0 won't be
  # evaluated.
  bresults <- vector(mode = "list", length = R)
  names(bresults) <- c(1:R)
  nonZero <- as.matrix((bParameterVectors[, , "t1"] + 
                        bParameterVectors[, , "b1"]) > 0)
  for (sample in 1:R) {
    bresults[[sample]] <- bParameterVectors[sample, nonZero[sample, ], "t1",
                                            drop=FALSE]
  }
  
  if (test == "ks.over") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"), 
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        bterms <- colnames(bresults[[sample]])
        btestGV <- testGV[bootIndices$test[as.numeric(sample), ], bterms,
                          drop=FALSE]
        bbackGV <- backGV[bootIndices$back[as.numeric(sample), ], bterms,
                          drop=FALSE]
        for (term in bterms) {
          test <- btestGV[[term]]
          back <- bbackGV[[term]]
          bresults[[sample]][[term]] <- ks.test(test, back,
                                                alternative = "less")$p.value
        }
      }
      bresults[chunk]
    }
  } else if (test == "ks.under") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"), 
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        bterms <- colnames(bresults[[sample]])
        btestGV <- testGV[bootIndices$test[as.numeric(sample), ], bterms,
                          drop=FALSE]
        bbackGV <- backGV[bootIndices$back[as.numeric(sample), ], bterms,
                          drop=FALSE]
        for (term in bterms) {
          test <- btestGV[[term]]
          back <- bbackGV[[term]]
          bresults[[sample]][[term]] <- ks.test(test, back,
                                           alternative = "greater")$p.value
        }
      }
      bresults[chunk]
    }
  } else if (test == "wilcox.over") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"), 
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        bterms <- colnames(bresults[[sample]])
        btestGV <- testGV[bootIndices$test[as.numeric(sample), ], bterms,
                          drop=FALSE]
        bbackGV <- backGV[bootIndices$back[as.numeric(sample), ], bterms,
                          drop=FALSE]
        for (term in bterms) {
          test <- btestGV[[term]]
          back <- bbackGV[[term]]
          bresults[[sample]][[term]] <- wilcox.test(test, back,
                                               alternative = "greater")$p.value
        }
      }
      bresults[chunk]
    }
  } else if (test == "wilcox.under") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"), 
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        bterms <- colnames(bresults[[sample]])
        btestGV <- testGV[bootIndices$test[as.numeric(sample), ], bterms,
                          drop=FALSE]
        bbackGV <- backGV[bootIndices$back[as.numeric(sample), ], bterms,
                          drop=FALSE]
        for (term in bterms) {
          test <- btestGV[[term]]
          back <- bbackGV[[term]]
          bresults[[sample]][[term]] <- wilcox.test(test, back,
                                               alternative = "less")$p.value
        }
      }
      bresults[chunk]
    }
  } else if (test == "kruskal") {
    bresults <- foreach(chunk = ichunk(1:R, chunkSize, mode = "character"), 
                        .combine = "c") %dopar% {
      for (sample in chunk) {
        bterms <- colnames(bresults[[sample]])
        btestGV <- testGV[bootIndices$test[as.numeric(sample), ], bterms,
                          drop=FALSE]
        bbackGV <- backGV[bootIndices$back[as.numeric(sample), ], bterms,
                          drop=FALSE]
        for (term in bterms) {
          test <- btestGV[[term]]
          back <- bbackGV[[term]]
          bresults[[sample]][[term]] <- kruskal.test(list(test, back))$p.value
        }
      }
      bresults[chunk]
    }
  }

  for (sample in 1:R) {
    bresults[[sample]] <- sort(bresults[[sample]])
  }
  return(bresults)
}



# -=-=-= Post initial statistical analysis =-=-=-


MultipleHypothesisCorrection <- function(pvalue, method = "BH") {
  # Applies the multiple hypothesis correction to the results list.
  #
  # Args:
  #   pvalue: (list) p-value of all terms.
  #   method: (char) method of correction to use. Accepts:
  #                   "bonferroni" - Bonferroni correction
  #                   "holm" - Holm (1979)
  #                   "hochberg" - Hochberg (1988)
  #                   "hommel" - Hommel (1988)
  #                   "BH" or "fdr" - Benjamini & Hochberg (1995)
  #                   "BY" - Benjamini & Yekutieli (2001)
  # Returns: 
  #   results: (list) p-value of all terms, corrected for multple
  #                   hypothesis.
  # 

  results <- p.adjust(pvalue, method = method)
  return(results)
}


Boot_MHCorrection <- function(bpvalue, method = "BH") {
  # Applies the multiple hypothesis correction to the results list.
  #
  # Args:
  #   pvalue: (list) p-value of all terms.
  #   method: (char) method of correction to use. Accepts:
  #                   "bonferroni" - Bonferroni correction
  #                   "holm" - Holm (1979)
  #                   "hochberg" - Hochberg (1988)
  #                   "hommel" - Hommel (1988)
  #                   "BH" or "fdr" - Benjamini & Hochberg (1995)
  #                   "BY" - Benjamini & Yekutieli (2001)
  # Returns: 
  #   results: (list) p-value of all terms, corrected for multple
  #                   hypothesis.
  # 

  R <- length(bpvalue)

  bresults <- vector(mode = "list", length = R)
  names(bresults) <- c(1:R)

  for (sample in 1:R) {
    bresults[[sample]] <- p.adjust(bpvalue[[sample]], method = method)
  }

  return(bresults)
}


SignificantSamples <- function(results, parameterVectors,
                               criticalValue = 0.05) {
  # Counts the number of significant occurrences of an ontology term among the
  # bootstrap samples.
  #
  # Args:
  #   results: (list) p-value of all terms, corrected for multple hypothesis.
  #   parameterVectors: (list) phyletic vector with the number of occurrences
  #                            of each term among the genomes in the groups.
  #                            Format: (t0, t1, b0, b1).
  #   criticalValue: (numeric) cutoff for the significance criteria. The terms
  #                            must show a value below it to be significant.
  # Returns: 
  #   count: (vector) number of significant occurrences among the bootstrap 
  #                   samples. Evaluates terms available in parameterVectors.
  # 
  
  R <- length(results)

  if (is.null(names(parameterVectors))) {  # Case when using the TS algorithm.
    terms <- character(0)
    for (i in 1:R) {
      terms <- union(terms, names(parameterVectors[[i]])) 
    }
    count <- rep.int(0, times = length(terms))
    names(count) <- terms
  } else {
    count <- rep.int(0, times = length(parameterVectors))
    names(count) <- names(parameterVectors)
  }
  
  for (sample in 1:R) {
    terms <- names(results[[sample]][results[[sample]] < criticalValue])
    count[terms] <- count[terms] + 1
  }

  return(count)
}


GetPhiCoefficient <- function(parameterVectors) {
  # Adds information of size effect to the results. Used in parametric tests.
  #
  # Args:
  #   parameterVectors: (list) phyletic vector with the number of occurrences
  #                            of each term among the genomes in the groups.
  #
  # Returns:
  #   results: (list) effect size (Phi Coefficient) of all terms.

  results <- vector(mode = "numeric", length = length(parameterVectors))
  names(results) <- names(parameterVectors)
  
  for (i in names(results)) {
    results[[i]] <- assocstats(matrix(parameterVectors[[i]], nrow = 2))$phi
  }

  return(results)

}


Get.r.Coefficient <- function (testGV, backGV, pvalue) {
  # Informs the size effect to the results. Used in non-parametric tests. The
  # effect size (r) is r = z/sqrt(n), where n is the total number of genomes.
  #
  # Args:
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  #   pvalue: (list) p-value of all terms.
  #
  # Returns:
  #   results: (numeric) effect size (r coefficient) of each term.


  testLen <- nrow(testGV)
  backLen <- nrow(backGV)

  wMean <- min(testLen, backLen) * (testLen + backLen + 1) / 2
  wSError <- sqrt(testLen * backLen * (testLen + backLen + 1) / 12)
  sqrtn <- sqrt(testLen + backLen)

  terms <- names(pvalue)
  ranks <- Map("c", testGV[, terms, drop=FALSE], backGV[, terms, drop=FALSE])
  ranks <- sapply(ranks, rank)
  
  # wStatistic is what we need to calculate z and, therefore, r
  if (testLen < backLen) {
    ranks <- ranks[1:testLen, , drop=FALSE]
    wStatistic <- colSums(ranks)
  } else if (backLen < testLen) {
    ranks <- ranks[(testLen + 1):(testLen + backLen), , drop=FALSE]
    wStatistic <- colSums(ranks)
  } else {
    ranks1 <- ranks[1:testLen, , drop=FALSE]
    ranks2 <- ranks[(testLen + 1):(testLen + backLen), , drop=FALSE]
    wStatistic <- unlist(Map("min", colSums(ranks1), colSums(ranks2)))
  }
  
  z <- (wStatistic - wMean) / wSError
  results <- abs(z/sqrtn)
  
  return(results)
  
#   label <- factor(c(rep("test", length(testGV[[names(pvalue)[1]]][, 1])),
#                     rep("back", length(backGV[[names(pvalue)[1]]][, 1]))))
#   sqrtn <- sqrt(length(label))
# 
#   for (i in names(pvalue)) {
#     values <- c(testGV[[i]][, 1], backGV[[i]][, 1])
# 
#     z <- statistic(wilcox_test(values ~ label), type = "standardized")
# #     z <- slot(w@statistic, "teststatistic")
#     pvalue[[i]] <- abs(z/sqrtn)
#   }
}


ConfidenceInterval <- function(effectSize, sd, type, test.name = NULL,
                               back.name = NULL, parameterVectors = NULL) {
  # Calculates the confidence interval for a specific effect size.
  #
  # Args:
  #   effectSize: (list) effect size (Phi Coefficient) of all terms.
  #   sd: (numeric) number of standard deviations for confidence.
  #   type: (char) whether the result was from a parametric (P) or from a 
  #                non-parametric (N) test.
  #   test.name: (char) name of the genomes used in the test group.
  #   back.name: (char) name of the genomes used in the background group.
  #   parameterVectors: (list) phyletic vector with the number of occurrences
  #                            of each term among the genomes in the groups.
  # Returns:
  #   confInterval: (matrix) confidence interval for each term.

  confInterval <- vector(mode = "list", length = length(effectSize))
  names(confInterval) <- names(effectSize)

#   colnames(parameterVectors) <- c("t1", "b1", "t0", "b0")
  
  if (type == "P") {
    ntest <- parameterVectors[[1]][1] + parameterVectors[[1]][3]
    nback <- parameterVectors[[1]][2] + parameterVectors[[1]][4]
  } else if (type == "N") {
    ntest <- length(test.name)
    nback <- length(back.name)
  }

  sigma <- sqrt(((ntest + nback) / (ntest * nback)) +
                 (effectSize^2 / (2 * (ntest + nback))))
  margin <- sd * sigma
  confInterval <- cbind(effectSize - margin, effectSize + margin)
                       
  return(confInterval)

}


CalculateFoldChange <- function(parameterVectors, type = "P",
                                test.name = NULL, back.name = NULL) {
  # Calculates the fold change between the two genome groups.
  #
  # Args:
  #   parameterVectors: (list) phyletic vector with the number of occurrences 
  #                          of each term among the genomes in the groups.
  #   type: (char) whether the result was from a parametric (P) or from a 
  #                non-parametric (N) test.
  #   test.name: (char) name of the genomes used in the test group, used if 
  #                     type is "N".
  #   back.name: (char) name of the genomes used in the background group,
  #                     used if type is "N".
  # Returns:
  #   foldChange: (vector) vector with the fold change for each GO.

  # Chech bug: altering names if they contain the character ":"
  # Problem with the make.names() function

  name <- names(parameterVectors)
  parameterVectors <- t(as.data.frame(parameterVectors, optional = TRUE))
  colnames(parameterVectors) <- c("t1", "b1", "t0", "b0")

  # For parametric (P) tests, we compare the change in proportion p = x/(n-x)
  # For non-parametric (N) ones, we compare the change in the population mean
  if (type == "P") {
    foldChange <- (parameterVectors[, "t1"] / parameterVectors[, "t0"]) / 
                  (parameterVectors[, "b1"] / parameterVectors[, "b0"])
  } else if (type == "N") {
    foldChange <- (parameterVectors[, "t1"] / length(test.name)) / 
                  (parameterVectors[, "b1"] / length(back.name))
  }

  names(foldChange) <- name
  return(foldChange)
}



# -=-=-=- Phylogenetically Independent Contrast analysis -=-=-=-


FindContrasts <- function(x, y, tree, method = "gls", denominator = 1) {



  # Produces a vector with the correlation of each ontology term with the 
  # attribute in question after correcting for phylogenetic bias (see
  # Felsenstein 1985 and APE package for details).
  #
  # Args:
  #   x: (vector) variable with the counting of an attribute of 
  #               interest, like G+C, gene count or longevity.
  #   y: (data.frame) table counting the presence of annotations
  #                   of an ontology in each genome.
  #   method: (char) method to use, allows "pearson", "spearman" and "kendall".
  #   denominator: (numeric) parameter for normalization of the y variable.
  # Returns:
  #   correlations: (vector) correlation of all listed ontology terms for the
  #                          attribute in question.
  
  if (method == "pic") {
    tmp_x <- as.vector(as.numeric(KOMODO2$x[,1]))
    names(tmp_x) <- rownames(KOMODO2$x)
    contrast_x <- pic(tmp_x, tree)
#  correlations <- vector(mode = "numeric", length = ncol(y))
    models <- vector(mode="numeric", length=ncol(y))
#  models2 <- vector(length=ncol(y))
#  models3 <- vector(mode="numeric", length=ncol(y))
    names(models) <- colnames(y)
#  names(models2) <- colnames(y)
#  names(models3) <- colnames(y)
#  names(correlations) <- colnames(y)
  
  # Normalizing 
    if (!is.null(denominator)) {
#    y <- as.data.frame(t(t(y) / denominator))
      y <- y / denominator
    }
  
    for (i in 1:ncol(y)) {
      tmp_y <- as.vector(as.numeric(y[, i]))
      names(tmp_y) <- rownames(x)
      contrast_y <- pic(tmp_y, tree)
      model <- lm(contrast_y ~ contrast_x + 0)
      models[[i]] <- summary(model)$coefficients[1,4]
    }
    models <- sort(models, decreasing = FALSE)
    return(models)
  }
  if (method == "gls") {
    tmp_x <- as.vector(as.numeric(KOMODO2$x[,1]))
    names(tmp_x) <- rownames(KOMODO2$x)
    models <- vector(mode="numeric", length=ncol(y))
    names(models) <- colnames(y)

  # Normalizing 
    if (!is.null(denominator)) {
      y_norm <- y / denominator
      y_norm = y_norm * 10^6 #getting counts per million to avoid false convergence (8) error from gls function for small values, see http://r.789695.n4.nabble.com/quot-False-convergence-quot-in-LME-td860675.html
    }

    for (i in 1:ncol(y_norm)) {
      tmp_y <- as.vector(as.numeric(y_norm[, i]))
      names(tmp_y) <- rownames(x)
      df1 <- as.data.frame(cbind(tmp_x, tmp_y))
      if (sum(tmp_y == 0) > 0) {
        model <-gls(tmp_y~tmp_x, data=df1, correlation=corPagel(1,tree), control = list(singular.ok = TRUE))
#    models2[[i]] <- model
        print(summary(model)$coefficients[2])
        models[[i]] <- as.numeric(summary(model)$coefficients[2])
      } else {
        models[[i]] <- 1
      }
    }
    models <- sort(models, decreasing = FALSE)
    return(models)
  #  return(correlations)
  }
}

# -=-=-=- Correlation analysis -=-=-=-


# Consider adding differentiation for x and y denominator column.
FindCorrelations <- function(x, y, method = "pearson", denominator = 1) {
  # Produces a vector with the correlation of each ontology term with the 
  # attribute in question.
  #
  # Args:
  #   x: (vector) variable with the counting of an attribute of 
  #               interest, like G+C, gene count or longevity.
  #   y: (data.frame) table counting the presence of annotations
  #                   of an ontology in each genome.
  #   method: (char) method to use, allows "pearson", "spearman" and "kendall".
  #   denominator: (numeric) parameter for normalization of the y variable.
  # Returns:
  #   correlations: (vector) correlation of all listed ontology terms for the
  #                          attribute in question.
  correlations <- vector(mode = "numeric", length = ncol(y))
  correlations.pvalue <- vector(mode = "numeric", length = ncol(y))
  names(correlations) <- colnames(y)
  names(correlations.pvalue) <- colnames(y)

  # Normalizing 
  if (!is.null(denominator)) {
    y <- as.data.frame(t(t(y) / denominator))
#    y <- y / denominator
  }
  for (i in 1:ncol(y)) {
    correlations[[i]] <- cor(x[rownames(y), 1], y[, i], method = method)
    correlations.pvalue[[i]] <- (cor.test(x[rownames(y), 1], y[, i], method = method))$p.value
  }

  correlations <- sort(correlations, decreasing = TRUE)
  correlations.pvalue <- sort(correlations.pvalue, decreasing = TRUE)
  
  results <- list("cor" = correlations, "cor.pvalues" = correlations.pvalue)
  
  return(results)
#  return(correlations)
}


AddVariableX <- function(x.path) {
  # Produces a vector with the variable x of a correlation analysis from a
  # tab-separated file with two columns: genome names and variable, no header. 
  #
  # Args:
  #   x.path: (char) path to the file with the variable data. 
  # 
  # Returns:
  #   x: (vector) vector with the data of variable x for correlation analysis.

  if (is.null(x.path) | file.exists(x.path) == FALSE) {
    return(NULL)
  }
  
  x.table <- read.table(x.path, sep = "\t", header = FALSE, colClasses = "character",
                        strip.white = TRUE, comment.char = "")
  x <- as.numeric(x.table[, 2])
  names(x) <- x.table[, 1]
  x <- as.data.frame(x)
  
  return(x)
}


VariableCorrelation <- function(x, results, method = "pearson", 
                                testGV = NULL, backGV = NULL) {
  # Produces a vector with the correlation of each ontology term with the 
  # variable x in question.
  #
  # Args:
  #   x: (vector) variable with the counting of an attribute of 
  #               interest, like G+C, gene count or longevity.
  #   results: (list) p-value of all terms.
  #   method: (char) method to use, allows "pearson", "spearman" and "kendall".
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  # Returns:
  #   correlations: (vector) correlation of all listed ontology terms for the
  #                          attribute in question.

  if (is.null(testGV) && is.null(backGV)) {
    return (NULL)
  }
  
  # check variable x: if it is a path to a tabular file, read it.
  if (is.character(x) == TRUE) {
    x <- AddVariableX(x)
  }
  if (is.null(x)) {
    return(NULL)
  }
  
  # obtain the frequency of every term among the genomes
  if (!is.null(testGV) && !is.null(backGV)) {
    y <- rbind(testGV[names(results)], backGV[names(results)])
  } else if (!is.null(testGV)) {
    y <- testGV[names(results)]
  } else if (!is.null(backGV)) {
    y <- backGV[names(results)]
  }
  
  variableCor <- FindCorrelations(x, y, method)
  return (variableCor)

}


PrintVariableCor <- function (x, ontology, results, annotation,
                              method = "pearson", testGV = NULL, 
                              backGV = NULL, outputName) {
  # Wrapper that prints the result of VariableCorrelation. Analogous to 
  # PrintResults().
  #
  # Args:
  #   x: (vector) variable with the counting of an attribute of 
  #               interest, like G+C, gene count or longevity.
  #   ontology: (char) which ontology is used: GO, KO, other.
  #   results: (list) p-value of all terms.
  #   method: (char) method to use, allows "pearson", "spearman" and "kendall".
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  #   outputName: (char) name of the output file.
  # Returns:
  #   correlations: (vector) correlation of all listed ontology terms for the
  #                          attribute in question.

  corr.dir <- file.path(KOMODO2$output.dir, "correlations", method)
  dir.create(corr.dir, recursive = TRUE, showWarnings = FALSE)
  
  variableCor <- VariableCorrelation(x, results, method, testGV, backGV)
  
  if (length(variableCor) != 0) {
    outputName <- file.path("correlations", method, outputName)
    PrintCResults(variableCor, annotation, outputName, ontology)    
  }
}




TermCorrelation <- function(term, results, method = "pearson",
                            testGV = NULL, backGV = NULL) {
  # Produces a vector with the correlation of each ontology term with the 
  # term in question.
  #
  # Args:
  #   term: (char) which term to correlate with the rest. 
  #   results: (list) p-value of all terms.
  #   method: (char) method to use, allows "pearson", "spearman" and "kendall".
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  # 
  # Returns:
  #   termCor: (vector) correlation of all listed ontology terms for the
  #                     term in question.

  if (is.null(testGV) && is.null(backGV)) {
    return (NULL)
  }

  # obtain the main term's frequency in each genome
  # obtain the frequency of every term among the genomes
  if (!is.null(testGV) && !is.null(backGV)) {
    termFreq <- rbind(testGV[term], backGV[term])
    allFreq <- rbind(testGV[names(results)], backGV[names(results)])
  } else if (!is.null(testGV)) {
    termFreq <- testGV[term]
    allFreq <- testGV[names(results)]
  } else if (!is.null(backGV)) {
    termFreq <- backGV[term]
    allFreq <- backGV[names(results)]
  }
  
  termCor <- FindCorrelations(termFreq, allFreq, method)
  return (termCor)

}


PrintSignTermCor <- function (ontology, results, annotation,
                              method = "pearson", testGV = NULL, 
                              backGV = NULL, criticalValue = 0.05) {
  # TermCorrelation() applied to all significant terms of an analysis. Prints
  # a file for each significant term.
  #
  # Args:
  #   ontology: (char) which ontology is used: GO, KO, other.
  #   results: (list) p-value of all terms.
  #   annotation: (list) translation of the ontology's terms.
  #   method: (char) method to use, allows "pearson", "spearman" and "kendall".
  #   testGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the test group.
  #   backGV: (list) list of the number of occurrences of each term among the
  #                  genomes in the background group.
  #   criticalValue: (numeric) cutoff for the significance criteria. The terms
  #                            must show a value below it to be significant.
  # Returns:
  #   none.


  if (!is.null(criticalValue)) {
    termsToPrint <- names(results[results < criticalValue])
  } else {
    termsToPrint <- names(results)
  }
  
  if (is.null(termsToPrint)) {
    return(NULL)
  }

  corr.dir <- file.path(KOMODO2$output.dir, "correlations", method)
  dir.create(corr.dir, recursive = TRUE, showWarnings = FALSE)
  
  for (term in termsToPrint) {
    termCor <- TermCorrelation(term, results, method, testGV, backGV)
    if (length(termCor) > 1) {
      outputName <- file.path("correlations", method, 
                              paste0(tolower(gsub(":", "", term)), "_cor.tsv"))
      PrintCResults(termCor, annotation, outputName, ontology)
    }
  }
}


# -=-=-=- Display output -=-=-=-


AnnotateResults <- function(results, ontology) {
  # Adds the field "Term" of each GO/KO ID to the results.
  #
  # Args:
  #   results: (list) p-value of all GO/KO IDs, corrected for multple
  #                   hypothesis.
  #   ontology: (char) which ontology is used: GO, KO, other.
  # Returns:
  #   annotation: (list) translation of the ontology's terms.

  ontology <- tolower(ontology)
  
  if (ontology == "go" | ontology == "gene ontology") {
    annotation <- as.list(Term(GOTERM[names(results)]))
  } else if (ontology == "kegg") {
    annotation <- as.list(KOMODO2$allKOs[names(results)])
  } else if (ontology == "other") {
    annotation <- KOMODO2$dictionary[names(results)]
  }
  
  return(annotation)

}


PrintResults <- function(KOMODO2, test, outputName) {
  # Generates a file with the results of the analysis.
  #
  # Args:
  #   KOMODO2: (list) list with all data that was generated.
  #   test: (char) which test's results to print.
  #   outputName: (char) name of the output file.
  # Returns:
  #   none.
  
  pvalue <- switch(test, chi2 = "pvalue.chi2",
                         f.over = "pvalue.f.over",
                         f.under = "pvalue.f.under")
  results <- switch(test, chi2 = "results.chi2",
                          f.over = "results.f.over",
                          f.under = "results.f.under")
  count <- switch(test, chi2 = "count.chi2",
                        f.over = "count.f.over",
                        f.under = "count.f.under")


  output <- vector(mode = "list", length = length(KOMODO2[[results]]))
  names(output) <- names(KOMODO2[[results]])

  if (KOMODO2$bootstrap > 1) {
    for (i in names(output)) {
      output[[i]] <- c(i, KOMODO2[[results]][[i]], KOMODO2[[pvalue]][[i]],
                       paste(KOMODO2[[count]][[i]], '/', KOMODO2$bootstrap,
                             sep = ""), KOMODO2$phi[[i]],
                       KOMODO2$foldChange.p[[i]],
                       KOMODO2$parameterVectors[[i]], KOMODO2$annotation[[i]])
    }
  } else if (KOMODO2$taxonSampling > 1) {
    for (i in names(output)) {
      output[[i]] <- c(i, KOMODO2[[results]][[i]], KOMODO2[[pvalue]][[i]],
                       paste(KOMODO2[[count]][[i]], '/', KOMODO2$taxonSampling,
                             sep = ""), KOMODO2$phi[[i]],
                       KOMODO2$foldChange.p[[i]],
                       KOMODO2$parameterVectors[[i]], KOMODO2$annotation[[i]])
    }
  } else {
    for (i in names(output)) {
      output[[i]] <- c(i, KOMODO2[[results]][[i]], KOMODO2[[pvalue]][[i]],
                       KOMODO2$phi[[i]], KOMODO2$foldChange.p[[i]],
                       KOMODO2$parameterVectors[[i]], KOMODO2$annotation[[i]])
    }
  }

  df.output <- data.frame(matrix(unlist(output), nrow = length(output),
                          byrow = T, dimnames = list(names(output))))

  ontology <- tolower(KOMODO2$ontology)
  if (ontology == "go" | ontology == "gene ontology") {
    id <- "GO"
  } else if (ontology == "kegg") {
    id <- "KO"
  } else if (ontology == "other") {
    id <- "ID"
  }

  if (KOMODO2$bootstrap > 1) {
    output.columns <- c(id, "qvalue", "pvalue",
                        paste("bootstrap(", KOMODO2$criticalValue, ')',
                        sep = ""), "ES", "fold change", "t1", "b1", "t0", "b0",
                        "annotation")
  } else if (KOMODO2$taxonSampling > 1) {
    output.columns <- c(id, "qvalue", "pvalue",
                        paste("TS(", KOMODO2$criticalValue, ')',
                        sep = ""), "ES", "fold change", "t1", "b1", "t0", "b0",
                        "annotation")
  } else {
    output.columns <- c(id, "qvalue", "pvalue", "ES", "fold change", "t1",
                        "b1", "t0", "b0", "annotation")
  }

  if (!is.null(KOMODO2$output.dir)) { 
    outputName <- file.path(KOMODO2$output.dir, outputName)
    dir.create(KOMODO2$output.dir, recursive = TRUE, showWarnings = FALSE)
  }
  write.table(df.output, file = outputName, quote = FALSE,
              row.names = FALSE, col.names = output.columns, sep = "\t")

}


PrintNPResults <- function(KOMODO2, test, outputName) {
  # Generates a file with the results of the analysis.
  #
  # Args:
  #   KOMODO2: (list) list with all data that was generated.
  #   test: (char) which test's results to print.
  #   outputName: (char) name of the output file.
  # Returns:
  #   none.
  
  pvalue <- switch(test, ks.over = "pvalue.ks.over",
                         ks.under = "pvalue.ks.under",
                         w.over = "pvalue.w.over",
                         w.under = "pvalue.w.under",
                         kruskal = "pvalue.kruskal")
  results <- switch(test, ks.over = "results.ks.over",
                          ks.under = "results.ks.under",
                          w.over = "results.w.over",
                          w.under = "results.w.under",
                          kruskal = "results.kruskal")
  count <- switch(test, ks.over = "count.ks.over",
                        ks.under = "count.ks.under",
                        w.over = "count.w.over",
                        w.under = "count.w.under",
                        kruskal = "count.kruskal")


  output <- vector(mode = "list", length = length(KOMODO2[[results]]))
  names(output) <- names(KOMODO2[[results]])

  if (KOMODO2$bootstrap > 1) {
    for (i in names(output)) {
      output[[i]] <- c(i, KOMODO2[[results]][[i]], KOMODO2[[pvalue]][[i]],
                       paste(KOMODO2[[count]][[i]], '/', KOMODO2$bootstrap,
                             sep = ""), KOMODO2$r[[i]],
                       KOMODO2$foldChange.n[[i]], KOMODO2$annotation[[i]])
    }
  } else if (KOMODO2$taxonSampling > 1) {
    for (i in names(output)) {
      output[[i]] <- c(i, KOMODO2[[results]][[i]], KOMODO2[[pvalue]][[i]],
                       paste(KOMODO2[[count]][[i]], '/', KOMODO2$taxonSampling,
                             sep = ""), KOMODO2$r[[i]],
                       KOMODO2$foldChange.n[[i]], KOMODO2$annotation[[i]])
    }
  } else {
    for (i in names(output)) {
      output[[i]] <- c(i, KOMODO2[[results]][[i]], KOMODO2[[pvalue]][[i]],
                       KOMODO2$r[[i]], KOMODO2$foldChange.n[[i]],
                       KOMODO2$annotation[[i]])
    }
  }

  df.output <- data.frame(matrix(unlist(output), nrow = length(output),
                          byrow = T, dimnames = list(names(output))))

  ontology <- tolower(KOMODO2$ontology)
  if (ontology == "go" | ontology == "gene ontology") {
    id <- "GO"
  } else if (ontology == "kegg") {
    id <- "KO"
  } else if (ontology == "other") {
    id <- "ID"
  }

  if (KOMODO2$bootstrap > 1) {
    output.columns <- c(id, "qvalue", "pvalue",
                        paste("bootstrap(", KOMODO2$criticalValue, ')',
                        sep = ""), "ES", "fold change", "annotation")
  } else if (KOMODO2$taxonSampling > 1) {
    output.columns <- c(id, "qvalue", "pvalue",
                        paste("TS(", KOMODO2$criticalValue, ')',
                        sep = ""), "ES", "fold change", "annotation")
  } else {
    output.columns <- c(id, "qvalue", "pvalue", "ES", "fold change",
                        "annotation")
  }

  if (!is.null(KOMODO2$output.dir)) { 
    outputName <- file.path(KOMODO2$output.dir, outputName)
    dir.create(KOMODO2$output.dir, recursive = TRUE, showWarnings = FALSE)
  }
  write.table(df.output, file = outputName, quote = FALSE,
              row.names = FALSE, col.names = output.columns, sep = "\t")

}


PrintCResults <- function(correlations, annotation, outputName, ontology) {
  # Generates a file with the results of the analysis.
  #
  # Args:
  #   correlations: (vector) value of the correlations between ontology's term
  #                          and variable.
  #   annotation: (list) translation of the ontology's terms.
  #   outputName: (char) name of the output file.
  #   ontology: (char) which ontology was used for correlation.
  # Returns:
  #   none.

#   output <- Map("c", names(correlations), correlations,
#                      annotation)

  output <- annotation[names(correlations)]
#  colnames(output)
  for (term in names(correlations)) {
#    term <- names(correlations)[i]
    output[[term]] <- c(term, correlations[[term]], annotation[[term]])
  }


  df.output <- data.frame(matrix(unlist(output), nrow = length(output),
                          byrow = T, dimnames = list(names(output))))

  ontology <- tolower(ontology)
  if (ontology == "go" | ontology == "gene ontology") {
    id <- "GO"
  } else if (ontology == "kegg") {
    id <- "KO"
  } else if (ontology == "other") {
    id <- "ID"
  }
  output.columns <- c(id, "CC", "annotation")

  if (!is.null(KOMODO2$output.dir)) { 
    outputName <- file.path(KOMODO2$output.dir, outputName)
    dir.create(KOMODO2$output.dir, recursive = TRUE, showWarnings = FALSE)
  }
  write.table(df.output, file = outputName, quote = FALSE,
              row.names = FALSE, col.names = output.columns, sep = "\t")

}



HierarchicalClustering <- function(KOMODO2, test = "f.over",
                                   annot = TRUE, criticalValue = 0.01,
                                   outputName = NULL, dendrogram = "both",
                                   rowv = TRUE, colv = TRUE, width = 3840,
                                   height = 3840, res = 200, units = "px",
                                   margins = c(5, 5), nResults = NULL) {
  # Generates a hierarchical clustering and presents it in a heatmap.
  # 
  # Args:
  #   KOMODO2: contains supertestGV and backGV.
  #   test: (char) name of the statistical test whose results this function
  #                will print.
  #   annot: (char) replaces GO terms by their annotation if TRUE. Does
  #                 nothing if FALSE. 
  #   criticalValue: (float) p-value cutoff of each term; those above it won't
  #                          be presented in the heatmap.
  #   outputName: (char) name of the PNG file to be produced. If it receives
  #                       no value (NULL), then KOMODO2 will only show the
  #                       heatmap in screen.
  #   dendrogram: (char) parameter for which dendrograms will appear, allows 
  #                      "both", "row", "column" and "none".
  #   rowv: (logical) whether to show a dendrogram for rows.
  #   colv: (logical) whether to show a dendrogram for columns.
  #   width: (integer) width of the PNG file. Default to 3840 units.
  #   height: (integer) height of the PNG file. Default to 3840 units.
  #   res: (integer) resolution of the PNG file. Default to 200.
  #   units: (char) units for width adn height. Default to "px" (pixels), can
  #                 also be "in" (inches), "cm" or "mm".
  #   margins: (vector) space for column and row names (terms and genomes).
  #   nResults: (numeric) number of significant results to show in the heatmap.
  # 
  # Returns:
  #   nothing, just prints the heatmap in a file.

  
  results <- switch(test, chi2 = "results.chi2",
                          f.over = "results.f.over",
                          f.under = "results.f.under",
                          ks.over = "results.ks.over",
                          ks.under = "results.ks.under",
                          w.over = "results.w.over",
                          w.under = "results.w.under",
                          kruskal = "results.kruskal")
  results <- KOMODO2[[results]]
  
  # Filtering results for clear visualization and sanity check.
  # The heatmap algorithm requires at least two terms to work.
  results <- results[results < criticalValue]
  if (length(results) < 2) {
    return(NULL)
  }
  
  # Ensures no error from trying to print more values than results available,
  # or from bypassing the above sanity check. Also, no floats, if they appear
  # for some reason.
  if (!is.null(nResults)) {
    nResults <- as.integer(max(2, nResults))
    results <- results[1:min(nResults, length(results))]
  }

  terms <- names(results)
  testData <- KOMODO2$testGV[, terms]
  backData <- KOMODO2$backGV[, terms]
  
  # For cases when some genomes have no element
  testData <- testData[rowSums(testData) > 0, ]
  backData <- backData[rowSums(backData) > 0, ]  
  
  type <- switch(test, chi2 = "P", f.over = "P", f.under = "P", ks.over = "N",
                       ks.under = "N", w.over = "N", w.under = "N",
                       kruskal = "N")
   
  if (type == "P") {
    testData <- testData / (KOMODO2$testElementCount[rownames(testData)]
                            - testData)
    backData <- backData / (KOMODO2$backElementCount[rownames(backData)]
                            - backData)
  }
  
#   rownames(testData) <- basename(KOMODO2$test.name)
#   rownames(backData) <- basename(KOMODO2$back.name)

  if (isTRUE(annot)){ 
    names(testData) <- unlist(KOMODO2$annotation[names(testData)])
    names(backData) <- unlist(KOMODO2$annotation[names(backData)])
  }

  fullData <- rbind(testData, backData)
  
  # Colors to help identifying the group of each genome in the heatmap
  groupColors <- c(rep("blue", times = nrow(testData)),
                   rep("red", times = nrow(backData)))
  names(groupColors) <- rownames(fullData)

  # preparing to print heatmap in a PNG file
  if (!is.null(outputName)) {
    if (!is.null(KOMODO2$output.dir)) { 
      output.dir <- file.path(KOMODO2$output.dir, "heatmaps")
      dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
      outputName <- file.path(output.dir, outputName)
    }
    png(outputName, width = width, height = height, res = res, units = units)
  }
  
  heatmap.2(as.matrix(fullData), scale = "column", trace = "none",
            lhei = c(2,8), col = greenred(75), symkey = FALSE,
            RowSideColors = groupColors, cexRow = 1, cexCol = 2,
            margins = margins, dendrogram = dendrogram,
            srtRow = NULL,
            srtCol = 45,
            Rowv = rowv, Colv = colv)

  if (!is.null(outputName)) {
    dev.off()
  }
}


HierarchicalClusteringCor <- function(y, KOMODO2, annot = TRUE,
                                      criticalValue, outputName = NULL,
                                      dendrogram = "both", rowv = TRUE, 
                                      colv = TRUE, width = 3840, height = 3840,
                                      res = 200, units = "px",
                                      margins = c(5, 5), nResults = NULL) {
  # Generates a hierarchical clustering and presents it in a heatmap.
  # 
  # Args:
  #   y: (data.frame) table counting the presence of annotations
  #                   of an ontology in each genome.
  #   KOMODO2: contains testGV and backGV.
  #   annot: (char) replaces GO terms by their annotation if TRUE. Does
  #                 nothing if FALSE. 
  #   criticalValue: (numeric) cutoff for the annotation; columns without 
  #                            elements above it won't appear in the heatmap.
  #   outputName: (char) name of the PNG file to be produced. If it receives
  #                      no value (NULL), then KOMODO2 will only show the
  #                      heatmap in screen.
  #   dendrogram: (char) parameter for which dendrograms will appear, allows 
  #                      "both", "row", "column" and "none".
  #   rowv: (logical) whether to show a dendrogram for rows.
  #   colv: (logical) whether to show a dendrogram for columns.
  #   width: (integer) width of the PNG file. Default to 3840 units.
  #   height: (integer) height of the PNG file. Default to 3840 units.
  #   res: (integer) resolution of the PNG file. Default to 200.
  #   units: (char) units for width adn height. Default to "px" (pixels), can
  #                 also be "in" (inches), "cm" or "mm".
  #   margins: (vector) space for column and row names (terms and genomes).
  #   nResults: (numeric) number of significant results to show in the heatmap.
  # 
  # Returns:
  #   nothing, just prints the heatmap in a file.

  
  # Filtering columns according to the criteria given by criticalValue
  termOccurrence <- colSums(y >= criticalValue)
  removal <- as.numeric(names(termOccurrence[termOccurrence == 0]))
  y <- y[, -removal]
  if (length(y) < 2) {
    return(NULL)
  }
  
  # Normalization with Z-score
  colSd <- apply(y, 2, sd)
  normalizedAnno <- (y - colMeans(y)) / colSd

  if (isTRUE(annot)){ 
    names(normalizedAnno) <- unlist(KOMODO2$annotation[names(y)])
  }

  # preparing to print heatmap in a PNG file
  if (!is.null(outputName)) {
    if (!is.null(KOMODO2$output.dir)) { 
      output.dir <- file.path(KOMODO2$output.dir, "heatmaps")
      dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
      outputName <- file.path(output.dir, outputName)
    }
    png(outputName, width = width, height = height, res = res, units = units)
  }

#   colorMap <- function(genomeName) { 
#     if (is.element(genomeName, basename(KOMODO2$test.name))) "blue" else "red"
#   }
#   groupColors <- sapply(rownames(fullData), colorMap))

  heatmap.2(as.matrix(normalizedAnno), scale = "column", trace = "none",
            lhei = c(2,8), col = greenred(75), symkey = FALSE,
            cexRow = 1, cexCol = 1, margin = c(25, 10), 
            dendrogram = dendrogram, Rowv = rowv, Colv = colv)

  if (!is.null(outputName)) {
    dev.off()
  }
}


PrintVolcanoPlot <- function(foldChange, results, outputName = NULL) {
  # Prints a volcano plot from the results of one of the statistic tests.
  #
  # Args:
  #   foldChange: (vector) vector with the fold change for each GO.
  #   results: (list) p-value of all GO IDs, corrected for multple hypothesis.
  #   outputName: (char) name of the output image file. If NULL, will only
  #                      print the volcano plot on the screen.
  # Returns:
  #   none.

  if (is.null(foldChange) | length(results) == 0) {
    return(NULL)
  }
  
  logFoldChange <- log2(foldChange)
  threshold <- as.factor(abs(logFoldChange) > 2 & unlist(results) < 0.01)
  info <- data.frame(pvalue = unlist(results), logFC = logFoldChange,
                     threshold = threshold, annotation = names(logFoldChange))

#   nameInThreshold <- function(name) { 
#     if (threshold[name] == TRUE) name else NA
#   }
#   info$annotation <- unlist(lapply(names(logFoldChange), nameInThreshold))

#   Cairo(file = outputName)
#   png(outputName)
  g <- ggplot(data = info, aes(x = logFC, y = -log10(pvalue),
                               colour = threshold)) +
       geom_point(alpha = 0.4, size = 1.75) +
       theme(legend.position = "none") +
       xlab("log2 fold change") + ylab("-log10 p-value")
#  +
#        geom_text(data = info, aes(x = logFC, y = jitter(-log10(pvalue), 200),
#                  label = annotation, size = 0.8),
#                  colour = "black", hjust = 1, vjust = 0)

  if (!is.null(outputName)) {
    if (!is.null(KOMODO2$output.dir)) { 
      output.dir <- file.path(KOMODO2$output.dir, "volcanoplots")
      dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
      outputName <- file.path(output.dir, outputName)
    }
    ggsave(g, filename = outputName, type = "cairo")
  }
  dev.off()

}



PrintVenn <- function() {
  # Prints a Venn diagram from the results of one of the statistic tests.
  # Internal use.
  #
  # Args:
  #
  # Returns:
  #   none.

  listSig <- list(Fisher = names(KOMODO2$results.f.over
                                 [KOMODO2$results.f.over < 0.01]),
                   KS = names(KOMODO2$results.ks.over
                             [KOMODO2$results.ks.over < 0.01]),
                   Wilcoxon = names(KOMODO2$results.w.over
                                 [KOMODO2$results.w.over < 0.01]))

#   listSig <- list(fisher.under = names(KOMODO2$results.f.under
#                                  [KOMODO2$results.f.under < 0.05]),
#                    ks.under = names(KOMODO2$results.ks.under
#                              [KOMODO2$results.ks.under < 0.05]),
#                    wilcox.under = names(KOMODO2$results.w.under
#                                  [KOMODO2$results.w.under < 0.05]))
# 
#   listSig <- list(chi2 = names(KOMODO2$results.chi2
#                           [KOMODO2$results.chi2 < 0.05]),
#                    f.over = names(KOMODO2$results.f.over
#                             [KOMODO2$results.f.over < 0.05]),
#                    f.under = names(KOMODO2$results.f.under
#                              [KOMODO2$results.f.under < 0.05]))
# 
#   listCor <- list(pearson = names(KOMODO2$correlations.pearson
#                              [KOMODO2$correlations.pearson >= 0.8]),
#                    spearman = names(KOMODO2$correlations.spearman
#                               [KOMODO2$correlations.spearman >= 0.8]),
#                    kendall = names(KOMODO2$correlations.kendall
#                              [KOMODO2$correlations.kendall >= 0.8]))
# 
#   listCor <- list(pearson = names(KOMODO2$correlations.pearson
#                              [1:200]),
#                    spearman = names(KOMODO2$correlations.spearman
#                               [1:200]),
#                    kendall = names(KOMODO2$correlations.kendall
#                              [1:200]))
}


PrintHist <- function(term, testGV, backGV, outputName, annotation = NULL,
                      title_size = 20, labels = c("test", "back")) {
  # Prints a density histogram showing the distribution of the two groups for
  # a specific term and prints it in a file.
  #
  # Args:
  #   term: (char) ontology term to print.
  #   testGV: KOMODO2$testGV, has the distribution for all terms in the test
  #           group.
  #   backGV: KOMODO2$backGV, has the distribution for all terms in the back
  #           group.
  #   outputName: (char) name of the file to be produced, recommended to end 
  #                      with the '.eps', '.pdf' or '.jpeg' file format.
  #   annotation: (list) list that maps ontology terms to an annotation; 
  #                      typically, the parameter is KOMODO2$annotation.
  #   title_size: (numeric) size of the title in the figure, purely cosmetic.
  #   labels:  (char) how to name the test and back groups, respectively.
  # Returns:
  #   none.
  
  if (!is.numeric(KOMODO2$testGV[[term]]) | 
      !is.numeric(KOMODO2$backGV[[term]])) {
    return (NULL)
  }

  final <- rbind(data.frame(term = testGV[[term]],
                            from = rep("test", length(testGV[[term]]))),
                 data.frame(term = backGV[[term]],
                            from = rep("back", length(backGV[[term]]))))

  if (is.null(annotation)) {
    ggplot(final, aes(term, fill = from)) + geom_density(alpha = 0.3) + 
    ggtitle(term) + theme(plot.title = element_text(size = title_size)) +
    xlab("Occurrences") + ylab("Density") +
    scale_fill_manual(values = c("blue", "red"),
                      name = "Genomic\nGroup",
                      breaks = c("test", "back"),
                      labels = labels)
  } else {
    ggplot(final, aes(term, fill = from)) + geom_density(alpha = 0.3) +
    ggtitle(paste0(term, " - ",  annotation[[term]])) + 
    theme(plot.title = element_text(size = title_size)) +
    xlab("Occurrences") + ylab("Density") +
    scale_fill_manual(values = c("blue", "red"),
                      name = "Genomic\nGroup",
                      breaks = c("test", "back"),
                      labels = labels)
  }

  if (!is.null(KOMODO2$output.dir)) { 
    output.dir <- file.path(KOMODO2$output.dir, "histograms")
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    outputName <- file.path(output.dir, outputName)
  }
  ggsave(file = outputName, type = "cairo")
  dev.off()

}

PrintSignificantHist <- function(results, criticalValue = 0.01, testGV, backGV,
                                 annotation = NULL, title_size = 20, 
                                 labels = c("test", "back")) {
  # Prints a density histogram showing the distribution of the two groups for
  # every significant term and prints each in a separated file.
  #
  # Args:
  #   results: (list) p-value of all GO/KO IDs.
  #   criticalValue: (numeric) critical significance value for a term to be
  #                            printed.
  #   testGV: KOMODO2$testGV, has the distribution for all terms in the test
  #           group.
  #   backGV: KOMODO2$backGV, has the distribution for all terms in the back
  #           group.
  #   annotation: (list) list that maps ontology terms to an annotation; 
  #                      typically, the parameter is KOMODO2$annotation.
  #   title_size: (numeric) size of the title in the figure, purely cosmetic.
  #   labels:  (char) how to name the test and back groups, respectively.
  # Returns:
  #   none.

  if (length(results[results < criticalValue]) == 0) {
    return(NULL)
  }
  
  for (term in names(results[results < criticalValue])) {
    outputName <- paste0(tolower(gsub(":", "", term)), ".jpeg")
    PrintHist(term, testGV, backGV, outputName,
              annotation, title_size, labels)
  }
}



# GenomeDendrogram <- function() {
# 
#   test <- rbind(KOMODO2$testGV, KOMODO2$backGV)
#   for (i in colnames(KOMODO2$testGV)) {
#     test[, i] <- as.numeric(as.logical(test[, i]))
#   }
#   
#   library("ape")
# #   hc <- hclust(dist(rbind(KOMODO2$testGV, KOMODO2$backGV)))
#   hc <- hclust(dist(test))
#   color_label <- c(rep("blue", 33), rep("red", 17))
# 
#   png("phylogenetic_KO_binary.png")
#   plot(as.phylo(hc), type = "cladogram", label.offset = 1,
#        cex = 0.6, tip.color = color_label)
#   dev.off()
# 
# }



# -=-=-=- Taxonomy Sampling (TS) -=-=-=-



TS_TaxonomyData <- function(idsFile, nodes) {
  # A function that provides the taxonomic information about the input IDs,
  # using the pre-processed information (nodes) from a NCBI taxonomy file
  # (downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/) referred by
  # the parameter taxondir. Generated by "getnodes(KOMODO2$taxondir)". 
  #
  # Args:
  #   idsFile: (char) either a path to the file with the Taxonomy IDs to
  #                   sample, or a character vector with these IDs.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure.
  # Returns:
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon.
  
  # Get ids to search
  if (!is.null(idsFile) & isTRUE(file.exists(idsFile))) {
    ids <- read.table(idsFile, sep = "\t", header = FALSE, comment.char = "",
                      colClasses = "integer", strip.white = TRUE)
    ids <- ids[, 1]
  } else {
    ids <- as.integer(idsFile)  # assume it's a test.name or back.name, for now
  }
  
  
  # Sanity check: unique ID inputs, remove duplicates.
  if (any(duplicated(ids))) {
    cat("Warning: some ids are repeated:", ids[duplicated(ids)], "\n")
    ids <- unique(ids)
  }
  
  # Sanity check: filter IDs that aren't part of NCBI notation.
  if (!all(is.element(ids, nodes$id))) {
    cat("Warning: the following inputs are not part of NCBI taxonomy IDs:",
        ids[!is.element(ids, nodes$id)], "\n")
    ids <- ids[is.element(ids, nodes$id)]
  }

  # Counting how often each taxonomy ID occurs in the dataset
  countIDs <- rep.int(0, nrow(nodes))
  names(countIDs) <- nodes$id
  nodes_parent <- nodes$parent
  names(nodes_parent) <- nodes$id
  searchIDs <- ids
  countIDs[as.character(searchIDs)] <- 1
  while (length(searchIDs) > 0) {
    searchIDs <- nodes_parent[as.character(searchIDs)]
    parentage <- table(searchIDs)
    countIDs[names(parentage)] <- countIDs[names(parentage)] + parentage
    searchIDs <- searchIDs[searchIDs != 1]
  }
  
  countIDs <- countIDs[countIDs > 0]
  return(countIDs)
  
}


TS_Algorithm_Strict <- function(taxon, m, nodes, countIDs, replacement = "no",
                                method = "balanced") {
  # An algorithm that receives a group of Taxonomy IDs and the size m of the
  # sample to obtain from them. Returns a vector with a maximized taxonomy 
  # diversity. Assumes that every input id is unique.
  # If replacement = "yes", then it may return repeated IDs if it increases
  # the taxonomy diversity.
  #
  # Args:
  #   taxon: (char/integer) Taxon from which to start sampling children taxa.
  #   m: (integer) size of the sample to generate.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure.
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  #   replacement: (char) whether the algorithm allows to repeat IDs in order
  #                       to maximize taxonomy diversity and to reach m IDs
  #                       in the output.
  #   method: (char) "balanced" - ensures that the taxa allocation differ by
  #                               at most 1 (when with replacement).
  #                  "randomized" - samples children uniformly, without trying
  #                                 to level child allocation.
  # Returns:
  #   outputIDs: (vector) vector of IDs with maximized taxonomy diversity.
  
  # Sanity check
  if (m <= 0) {
    Print("Error: m less or equal than zero during recursion.\n")
    exit(0)
  }
  
  # First step: Find the sub-taxa (children nodes) of the current taxon
  # Basically, nodes[nodes[, 2] == taxon, 1]
  taxon <- as.integer(taxon)
  children <- nodes[nodes[, 2] == taxon & nodes[, 1] != taxon, 1]
  children <- intersect(children, names(countIDs))
  
  # Condition to end recursion
  if (length(children) == 0) {
    if (replacement == "no") {
      return(as.character(taxon))
    } else {
      return(rep(as.character(taxon), m))
    }
  }
  
  m_i <- rep.int(0, length(children))
  names(m_i) <- children
  
  #  balanced - ensure two taxa allocation (m_i) differ at most by 1 
  #  randomized - fully random child sampling, uniform distribution among taxa
  if (method == "balanced") {
    m_i <- m_i + floor(m / length(children))
    sampledChildren <- sample(children, m - sum(m_i))
    m_i[sampledChildren] <- m_i[sampledChildren] + 1
  } else if (method == "randomized") {
    m_i <- table(sample(children, m, replace = TRUE))
  }
  
  outputIDs <- character(0)
  for (id in names(m_i)) {
    if (m_i[id] == 0) { next; }
    outputIDs <- c(outputIDs,
                   TS_Algorithm_Strict(id, m_i[id], nodes, countIDs,
                                       replacement, method))
  }
  
  return(outputIDs)
}



TS_Algorithm <- function(taxon, m, nodes, countIDs, replacement = "no", 
                         method = "balanced") {
  # An algorithm that receives a group of Taxonomy IDs and the size m of the
  # sample to obtain from them. Returns a vector with a maximized taxonomy 
  # diversity. Assumes that every input id is unique.
  # If replacement = "yes", then it may return repeated IDs if it increases
  # the taxonomy diversity.
  #
  # Args:
  #   taxon: (char/integer) Taxon from which to start sampling children taxa.
  #   m: (integer) size of the sample to generate.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure.
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  #   replacement: (char) whether the algorithm allows to repeat IDs in order
  #                       to maximize taxonomy diversity and to reach m IDs
  #                       in the output.
  #   method: (char) "balanced" - ensures that the taxa allocation differ by
  #                               at most 1 (when with replacement).
  #                  "randomized" - samples children uniformly, without trying
  #                                 to level child allocation.
  # Returns:
  #   outputIDs: (vector) vector of IDs with maximized taxonomy diversity.
  
  # Sanity check
  if (m <= 0) {
    Print("Error: m less or equal than zero during recursion.\n")
    exit(0)
  }
  
  # First step: Find the sub-taxa (children nodes) of the current taxon
  # Basically, nodes[nodes[, 2] == taxon, 1]
  taxon <- as.integer(taxon)
  children <- nodes[nodes[, 2] == taxon & nodes[, 1] != taxon, 1]
  children <- intersect(children, names(countIDs))
  
  # Condition to end recursion
  if (length(children) == 0) {
    if (replacement == "no") {
      return(as.character(taxon))
    } else {
      return(rep(as.character(taxon), m))
    }
  }
  
  childrenCount <- countIDs[as.character(children)]
  m_i <- rep.int(0, length(childrenCount))
  names(m_i) <- names(childrenCount)
  
  
  #  Balanced - only sample if 
  #   length(childrenCount[childrenCount > m_i]) < remaining_m
  #
  #  Semi-balanced - as it is
  #  Child sampling - fully random child sampling
  
  if (method == "balanced") {
    while (m > 0 & length(childrenCount[childrenCount > m_i]) <= m) {
      child <- names(childrenCount[childrenCount > m_i])
      m_i[child] <- m_i[child] + 1
      m <- m - length(child)
    }
    child <- sample(names(childrenCount[childrenCount > m_i]), m)
    m_i[child] <- m_i[child] + 1
    # m <- 0
    
  } else if (method == "randomized") {
    # If we don't have the m fully distributed over m_i (children),
    # choose a random child from childrenCount that still has taxa
    # available to choose.
    while (m > 0 & length(childrenCount[childrenCount > m_i]) > 0) {
      child <- sample(names(childrenCount[childrenCount > m_i]), 1)
      m_i[child] <- m_i[child] + 1
      m <- m - 1
    }
  }
  
  outputIDs <- character(0)
  for (id in names(m_i)) {
    if (m_i[id] == 0) { next; }
    outputIDs <- c(outputIDs,
                   TS_Algorithm(id, m_i[id], nodes, countIDs, replacement,
                                method))
  }
  
  return(outputIDs)
}




TS_ParAlgorithm <- function(m, nodes, countIDs, replacement = "yes",
                            taxonSampling = 100) {
  # Wrapper for parallel processing of multiple instances of TS_Algorithm.
  #
  # Args:
  #   m: (integer) size of the sample to generate.
  #   nodes: (data.frame) pre-processed information about the NCBI taxonomy
  #                       structure.
  #   countIDs: (vector) count of how many taxnomoy IDs belong to each taxon,
  #                      created by TS_TaxonomyData().
  #   replacement: (char) whether the algorithm allows to repeat IDs in order
  #                       to maximize taxonomy diversity and to reach m IDs
  #                       in the output.
  #   taxonSampling: (integer) number of occurrences
  # Returns:
  #   some.names: (list) sample vectors of IDs maximizing taxonomy diversity.

  some.names <- foreach (sample = iter(1:taxonSampling), 
                         .export = "TS_Algorithm") %dopar% {
    TS_Algorithm(1, m, nodes, countIDs, replacement)
  }
  
  return(some.names)
}


TSAddGenomeVectors <- function(someAnno, some.names) {
  # Wrapper for parallel processing of multiple instances of AddGenomeVectors
  # with some.names list (created by TS_ParAlgorithm) as input.
  #
  # Args:
  #   someAnno: (list) annotation from the original input data.
  #   some.names: (list) sample vectors of IDs with maximized taxonomy diversity.
  # Returns:
  #   someGV: (list) list with the phyletic vector of each genome of the input
  #                  group.
  
  ids <- unique(unlist(some.names))
  someGV <- AddGenomeVectors(someAnno, ids)
  
  return(someGV)

}


TSGroupElementCount <- function(someAnno, genome.names = NULL,
                                method = "default") {
  # Wrapper for GroupElementCount when using the TS algorithm. Returns a list
  # of vectors with the number of genes in each genome of a group. 
  # 
  # Args:
  #   someAnno: list of genomes, each with a data frame that maps
  #             each gene to its annotations (GO, KO). 
  #   genome.names: (vector) names of the genomes to count elements. In this
  #                          case, this is produced by TS_ParAlgorithm().
  #   method: method defining whether KOMODO2 must consider all elements in all
  #           genomes (default), or treat each element independently of others
  #           (experiment). The latter is an experiment for alignments.
  # Returns:
  #   tsElementCount: vector with the number of elements in each genome of 
  #                   someAnno, listed by the samples of genome.names.
  
  if (is.null(genome.names) | length(genome.names) == 0) {
    return(0)
  }
  
  tsElementCount <- foreach (some.name = iter(genome.names),
                             .export = "GroupElementCount") %dopar% {
    GroupElementCount(someAnno, some.name, method)
  }
  
  return(tsElementCount)
}


TSParameterVectors <- function(testGV, backGV,
                               tsTestElementCount, tsBackElementCount) {
  # Wrapper for ParameterVectors when using the TS algorithm. Creates a summary
  # of both test and background groups about the number of annotated genes for
  # each GO term.
  #
  # Args: 
  #   testGV: (data.frame) number of occurrences of each term among the genomes
  #                        of the test group.
  #   backGV: (data.frame) number of occurrences of each term among the genomes
  #                        of the background group.
  #   tsTestElementCount: vector with the number of elements in the genomes of 
  #                       the samples from the test group.
  #   tsBackElementCount: vector with the number of elements in the genomes of 
  #                       the samples from the background group.
  #
  # Returns:
  #   tsParameterVectors: (list) list with the parameter vectors for each
  #                              sample in a group.
  
  tsParameterVectors <- foreach (i = iter(1:length(tsTestElementCount)),
                                 .export = "ParameterVectors") %dopar% {
    ParameterVectors(testGV[names(tsTestElementCount[[i]]), ],
                     backGV[names(tsBackElementCount[[i]]), ],
                     tsTestElementCount[[i]], tsBackElementCount[[i]])
  }

  return(tsParameterVectors)
}


TSStatisticalTest <- function(tsParameterVectors, test = "fisher.over") {
  # Computes the chosen statistical test between two phyletic vectors.
  #
  # Args:
  #   tsParameterVectors: (list) list with the parameter vectors for each
  #                              sample in a group.
  #   test: (char) the statistical test to perform.
  # Returns:
  #   tsPvalues: (numeric) p-value of all terms.
  
  tsPvalues <- foreach (parameterVectors = iter(tsParameterVectors),
                        .export = "StatisticalTest", 
                        .packages = c("itertools", "foreach")) %dopar% {
    StatisticalTest(parameterVectors, test)
  }
  
  return(tsPvalues)
}


TSDistributions <- function(tsParameterVectors, testGV, backGV,
                            test.name, back.name, test = "ks.over") {
  # Computes the chosen statistical test between two phyletic vectors.
  #
  # Args:
  #   tsParameterVectors: (list) list with the parameter vectors for each
  #                              sample in a group.
  #   testGV: (data.frame) number of occurrences of each term among the genomes
  #                        of the test group.
  #   backGV: (data.frame) number of occurrences of each term among the genomes
  #                        of the background group.
  #   test.name: (char) name of the genomes used in the test group.
  #   back.name: (char) name of the genomes used in the background group.
  #   test: (char) the statistical test to perform.
  # Returns:
  #   results: (list) p-value of all terms.
  
  results <- foreach (i = iter(1:length(tsParameterVectors)),
                      .export = "CompareDistributions", 
                      .packages = c("itertools", "foreach")) %dopar% {
    CompareDistributions(tsParameterVectors[[i]], 
                         testGV[test.name[[i]], ],
                         backGV[back.name[[i]], ],
                         test)
  }
  
  return(results)
}


TS_MHCorrection <- function(tsPvalue, method = "BH") {
  # Applies the multiple hypothesis correction to the results list.
  #
  # Args:
  #   tsPvalue: (list) p-value of all terms.
  #   method: (char) method of correction to use. Accepts:
  #                   "bonferroni" - Bonferroni correction
  #                   "holm" - Holm (1979)
  #                   "hochberg" - Hochberg (1988)
  #                   "hommel" - Hommel (1988)
  #                   "BH" or "fdr" - Benjamini & Hochberg (1995)
  #                   "BY" - Benjamini & Yekutieli (2001)
  # Returns: 
  #   results: (list) p-value of all terms, corrected for multple
  #                   hypothesis.
  # 
  
  results <- foreach (pvalue = iter(tsPvalue)) %dopar% {
    p.adjust(pvalue, method = method)
  }
  
  return(results)
}



# -=-=-=- Annotation Enrichment Analysis (AEA) -=-=-=-


#ListOffspring <- function() {
#  # Prepares a listing of ancestors for each GO ID. Used to vectorize later 
#  # operations.
#  # 
#  # Args:
#  #   none.
#  # Returns:
#  #   allOffspring: (list) All GOXXOFFSPRING combined.
#
#  allOffspring <- c(as.list(GOBPOFFSPRING), as.list(GOCCOFFSPRING),
#                    as.list(GOMFOFFSPRING))
#  
#  return(allOffspring)
#}


# sketch: 
# - Create a list of offspring GO terms (branch, all descendants)
# - for every GO term, search offspring and determine parameters:
#   - number of genes in signature
#   - number of annotations to signature
#   - number of terms in branch
#   - number of annotations to branch
#   - number of annotations between signature and branch (Mgt)
#   - number of annotations between top signature and branch (~Mgt)
# - estimate pvalue as p(Mgt) = P(~Mgt >= Mgt)

# - search for significant GO terms only?
# - search using pre or pos FindAncestor()?
# - Remember to compare for KEGG Orthology later


#AnnotationEnrichmentAnalysis <- function(testGV, backGV, allOffspring, 
#                                         reps = 100) {
#  # Implementation of the Annotation Enrichment Analysis algorithm 
#  # (see doi:10.1038/srep04191). 
#  # 
#  # Args:
#  #   testGV: (data.frame) number of occurrences of each term among the genomes
#  #                        of the test group.
#  #   backGV: (data.frame) number of occurrences of each term among the genomes
#  #                        of the background group.
#  #   allOffspring: (list) All GOXXOFFSPRING combined.
#  #   reps: (integer) how many samples to generate in order to compute the 
#  #                   pvalue.
#  # Returns:
#  #   
#
#
##  # Filtering terms in testGV and backGV without occurrences
##  t1 <- colSums(testGV)
##  b1 <- colSums(backGV)
##
##  nonZero <- t1 > 0 | b1 > 0
##  t1 <- t1[nonZero, drop=FALSE]
##  b1 <- b1[nonZero, drop=FALSE]
#  
#  # The algorithm itself
#  genesInSignature <- nrow(testGV)
#  annotationsToSignature <- sum(rowSums(testGV))
#  
#  fullGV <- rbind(testGV, backGV)
#  nAnnotations <- colSums(fullGV)
#  names(nAnnotations) <- rownames(fullGV)
#  
#  result <- rep(0, ncol(testGV))
#  names(result) <- colnames(testGV)
#
#  for (term in colnames(testGV)) {
#    termsInBranch <- unique(c(term, unlist(allOffspring[term])))
#    termsInBranch <- termsInBranch[!is.na(termsInBranch)]
#    annotationsToBranch <- sum(rowSums(fullGV[, termsInBranch]))
#    
#    genomeShuffle <- sample(c(rownames(testGV), rownames(backGV)))
#    termShuffle <- sample(colnames(testGV))
#    
#    observed <- sum(rowSums(testGV[, termsInBranch]))
#    extremeValues <- 0
#    
#    for (k in 1:reps) {
#      # Consider logarithmic approach
#      i <- 1
#      j <- 1
#      while (sum(nAnnotations[genomeShuffle[1:i]]) < annotationsToSignature) {
#        i <- i + 1
#      }
#      while (sum(termShuffle[1:j]) < annotationsToBranch) {
#        j <- j + 1
#      }
#      annoToShuffle <- sum(rowSums(fullGV[genomeShuffle[1:i], termShuffle[1:j]]))
#      if (annoToShuffle >= observed) {
#        extremeValues <- extremeValues + 1
#      }
#    }
#    result[term] <- extremeValues / reps
#  }
#  
#  
#}


