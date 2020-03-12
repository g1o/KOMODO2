# Small auxiliary functions (not exported to the package namespace)

# Read a komodo2 list from a key-value file.
read_komodo2_file <- function(file.path){
  df <- utils::read.table(file = file.path, header = FALSE,
                          strip.white = TRUE,
                          comment.char = "#", sep = '=' ,
                          col.names = c('Key', 'Value'),
                          stringsAsFactors = FALSE)

  outlist <- as.list(df$Value)
  names(outlist) <- as.character(df$Key)
  return(outlist)
}


# Function to parse the tab format from Uniprot
parse_GenomeMap <- function(genome, column, split = " *; *") {
  genomeMap <- strsplit(x = genome[, column], split = split)
  names(genomeMap) <- rownames(genome)
  return(genomeMap)
}


# Function to generate and plot a taxonomic tree
GenerateTree <- function(taxonIds, db = "ncbi") {
  taxize_class <- taxize::classification(taxonIds, db = db)
  taxize_tree  <- taxize::class2tree(taxize_class, check = TRUE)
  # taxize::plot.classtree(taxize_tree)
  invisible(taxize_tree)
}


# Adds the field "Term" of each GO/KO ID to the results.
#
# Args:
#   defs: (list) a KOMODO2-type list object (generated internally by the
#   KOMODO2 functions)
#   results: (list)  MHT-corrected p-values of all GO/KO IDs
#   ontology: (char) which ontology is used: GO, KO, other.
# Returns:
#   annotation: (list) translation of the ontology's terms.
AnnotateResults <- function(defs, results, ontology) {

  x <- switch(tolower(ontology),
              "go"            = as.list(AnnotationDbi::Term(GO.db::GOTERM[names(results)])),
              "gene ontology" = as.list(AnnotationDbi::Term(GO.db::GOTERM[names(results)])),
              "kegg"          = as.list(defs$allKOs[names(results)]),
              "other"         = defs$dictionary[names(results)])

  return(x)
}
