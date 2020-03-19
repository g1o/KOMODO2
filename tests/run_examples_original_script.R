# Assumes that we're working under the KOMODO2 development environment
# (i.e., that the current working directory is the main KOMODO2 project folder)

# get list of parameter files
defs.list <- normalizePath(dir("./zz_original_scripts/bin/scripts/parameters_validation",
                               pattern = "parameters_",
                               full.names = TRUE))

cwd <- getwd()
setwd("./zz_original_scripts/bin/scripts/")
for (i in 2:length(defs.list)){#seq_along(defs.list)){
  cat("\n\n-----\nRunning file", i, "\n")

  # Load specs file
  source(defs.list[i])
  source("load.R")
  source("clean.R")
  source("do.R")

  saveRDS(KOMODO2, paste0(KOMODO2$output.dir, "final_KOMODO2_list.rds"))
}

setwd(cwd)
