# Run the package with all test scripts available in the
# folder data_files/parameters_validation

defs.list <- dir("./data_files/parameters_validation",
                 pattern = "parameters_",
                 full.names = TRUE)

for (i in 8:length(defs.list)){ #seq_along(defs.list)){
  cat("\n\n-----\nRunning file", i, "\n")
  defs <- run_KOMODO2(defs.list[i], cores = 3)
  fn   <- paste0(defs$output.dir, "final_defs.rds")
  saveRDS(defs, fn)
}
