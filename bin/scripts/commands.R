#get children of GO term, useful to grep for gene symbols in raw interpro files

GO_ID <- c("GO:0050793")

outfile <- paste0(GO_ID, "_")

if (is.character(GOMFOFFSPRING[[GO_ID]])) {
  outfile <- paste0(outfile, "GOMFOFFSPRING.txt")
  GO_ids <- GOMFOFFSPRING[[GO_ID]]

} else if (is.character(GOBPOFFSPRING[[GO_ID]])) {
  outfile <- paste0(outfile, "GOBPOFFSPRING.txt")
  GO_ids <- GOBPOFFSPRING[[GO_ID]]

} else if (is.character(GOCCOFFSPRING[[GO_ID]])) {
  outfile <- paste0(outfile, "GOCCOFFSPRING.txt")
  GO_ids <- GOCCOFFSPRING[[GO_ID]]
} else {
  outfile <- paste0(outfile, "DEURUIM.txt")
}

print(outfile)
print(GO_ids)
