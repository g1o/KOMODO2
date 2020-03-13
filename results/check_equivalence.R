# Validate equivalence of implementations

data1 <- readRDS("./results/test_defs.rds")    # created using package
data2 <- readRDS("./results/test_KOMODO2.rds") # created using original scripts

allnames <- unique(c(names(data1), names(data2)))

# Check which names are present/missing in each list
df <- data.frame(Name = allnames, data1 = "", data2 = "",
                 stringsAsFactors = FALSE)
for (i in seq_along(allnames)){
  df$data1[i] <- df$Name[i] %in% names(data1)
  df$data2[i] <- df$Name[i] %in% names(data2)
}

df

# Fix fields that changed name between inplementations
names(data2)[which(names(data2) == "annotation_files_dir")] <- "annotation.files.dir"
names(data2)[which(names(data2) == "tree_path")] <- "tree.path"
names(data2)[which(names(data2) == "tree_type")] <- "tree.type"
names(data2)[which(names(data2) == "linear_model_cutoff")] <- "linear.model.cutoff"
names(data2)[which(names(data2) == "contrasts.cor")] <- "annotation.contrasts"
names(data2)[which(names(data2) == "sum.cor")] <- "annotation.sum"
names(data2)[which(names(data2) == "cv.cor")] <- "annotation.cv"
names(data2)[which(names(data2) == "sd.cor")] <- "annotation.sd"

# Fix fields that were added/removed and do not change the final result:
data1$MHT.method <- NULL
data2$cl <- NULL
data2$yElementCount <- NULL

# Check that all naming inconsistencies were treated:
all(sort(names(data1)) == sort(names(data2)))

# put both list in the same order of field names
data1 <- data1[names(data2)]
all(names(data1) == names(data2))



# Now check if data types are different (prior to testing for equivalence):
df <- data.frame(Name = names(data1),
                 ClassInData1 = "",
                 ClassInData2 = "",
                 stringsAsFactors = FALSE)

df$ClassInData1 <- sapply(data1, class)
df$ClassInData2 <- sapply(data2, class)

df$SameClass <- (df$ClassInData1 == df$ClassInData2)
df


# Fix the two known class differences (no consequence in final analysis):
tmp <- names(data2$heterogeneity)
data2$heterogeneity <- as.numeric(data2$heterogeneity)
names(data2$heterogeneity) <- tmp

tmp <- names(data2$mode)
data2$mode <- as.numeric(data2$mode)
names(data2$mode) <- tmp

all(sapply(data2, class) == sapply(data1, class))


# Now finally check equivalence:
which(!mapply(identical, data1, data2))

# Ignore the following fields (differences related to file/folder paths
# and other non-consequential differences)
# - annotation.files.dir
# - output.dir
# - dataset.info
# - tree.path
# - cores

# Verify the remaining ones manually
# - ontology: OK
# (difference is: "GO" vs "go")
data1$ontology
data2$ontology

# - x: OK
# (difference is: column name, not used in the analysis)
data1$x
data2$x
colnames(data1$x) <- colnames(data2$x)
identical(data1$x, data2$x)

# - greaterthanzero: OK
# Differences are:
#    1) the numbers are now normalized by the length of their respective "y"
#       entry)
#    2) the way greaterthanzero was calculated in the original scripts was
#       wrong (see code below):
#
# tmp <- as.numeric(lapply(KOMODO2$y, greaterthanzero))
# tmp <- (mapply("-", length(KOMODO2$y.name), tmp)) # <------ DOESN'T MAKE SENSE
# names(tmp) <- names(KOMODO2$y)
# KOMODO2$greaterthanzero <- tmp

data1$greaterthanzero[1:5]
data2$greaterthanzero[1:5]

ly <- unique(sapply(defs$y, length)) # check that all y entries have the same length
tmp <- (ly - data2$greaterthanzero) / ly

identical(tmp, data1$greaterthanzero)


# - heterogeneity: OK
# Difference is: the numbers are now normalized by the length of their
# respective "y" entry)

identical(data1$heterogeneity, data2$heterogeneity / ly)



# - mode: OK
# (difference is: some entries have more than 1 modal value, routines are
# picking up different signals. Should not have any effect on the final
# analysis)
idx <- which(data1$mode != data2$mode)

# Get all modes of the observations detected as different above
tmp <- lapply(data2$y[idx],
              function(v){
                tmp <- sort(table(v), decreasing = TRUE)
                tmp <- as.numeric(names(tmp[tmp == tmp[1]]))
              })

# Check that in all cases we have more than one mode
min(sapply(tmp, length))

# Check that the mode detected in data1 is amongst the modes from data2
all(mapply(function(a,b){a %in% b},
           a = data1$mode[idx],
           b = tmp))



# Done!
