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




#to get the intersection of some usefult flat files produced by KOMODO2
awk -F"\t" '{if ($2 < 0.2){print $1}}' contrasts_corrected.tsv | sort > cont_2
awk -F"\t" '{if ($2 < 0.2){print $1}}' p_corr_qvalues_results.tsv | sort > pear_2
awk -F"\t" '{if ($2 > 20){print $1}}' sum.tsv | sort > sum
awk -F"\t" '{if ($2 > 1){print $1}}' sd.tsv | sort > sd

comm -1 -2 cont_2 pear_2 > comm_cont_pear_2
comm -1 -2 comm_cont_pear_2 sum > comm_cont_pear_sum_2
comm -1 -2 comm_cont_pear_sum_2 sd > comm_cont_pear_sum_sd_2

awk -F"\t" '{if ($2 < 0.1){print $1}}' contrasts_corrected.tsv | sort > cont_1
awk -F"\t" '{if ($2 < 0.1){print $1}}' p_corr_qvalues_results.tsv | sort > pear_1
awk -F"\t" '{if ($2 > 20){print $1}}' sum.tsv | sort > sum
awk -F"\t" '{if ($2 > 1){print $1}}' sd.tsv | sort > sd

comm -1 -2 cont_1 pear_1 > comm_cont_pear_1
comm -1 -2 comm_cont_pear_1 sum > comm_cont_pear_sum_1
comm -1 -2 comm_cont_pear_sum_1 sd > comm_cont_pear_sum_sd_1

