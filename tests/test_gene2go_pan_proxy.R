# Reinstall first:
devtools::install_github("fcampelo/KOMODO2")

library(KOMODO2)

install_and_update_packages()

# Retrieve fresh data files:
retrieve_data_files("./data_folder")

# Run the example:
defs <- "./data_folder/parameters_validation/parameters_gene2GO_Pan_proxy.txt"
defs <- run_KOMODO2(defs, cores = parallel::detectCores() - 1)
