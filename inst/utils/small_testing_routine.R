# library(devtools)
# devtools::install_github("fcampelo/KOMODO2")
#
# library(KOMODO2)
# # Install/update Bioconductor package dependencies for KOMODO2
# KOMODO2::install_and_update_packages()
#
# # Download required files
# KOMODO2::retrieve_data_files(target.dir = "./data_folder")
#
# defs <- "./data_folder/parameters_validation/parameters_gene2GO_Pan_proxy.txt"
# res  <- KOMODO2::run_KOMODO2(defs, cores = parallel::detectCores() - 1,
#                              render.report = FALSE)
#
# # Un-comment to remove the non-CRAN packages that were installed (if you want)
# # remove.packages(c("KOMODO2", "AnnotationDbi", "KEGGREST", "GO.db"))
