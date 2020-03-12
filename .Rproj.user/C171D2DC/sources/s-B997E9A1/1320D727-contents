#This is a sample parameter file for KOMODO2. Most variables are
# self-explanatory, and a formal definition of each can be found in KOMODO2
#documentation.
annotation.files.dir = "" # (required, string) - Folder where annotation files are located.
output.dir = ""           # (required, string) - output folder for results
dataset.info = ""         # (required, string) - genome metadata file, it should contain at least:
                          #    1) path for annotation data (if "annotation.files.dir" not provided OR file names (if "annotation.files.dir) is provided. Please notice this information should be the first column in metadata file;
                          #    2) phenotype data (numeric, this is the value KOMODO2 uses to rank species when searching for associations)
                          #    3) normalization data (numeric, this is the value KOMODO2 uses as a denominator to compute annotation term frequencies to remove potential biases caused by, for instance, overannotation of model organisms or large differences in the counts of genomic elements). Please notice KOMODO2 does not require normalization data for GO, as it computes the total number of GO terms per species and uses it as a normalizing factor.
x.column =                # (required, numeric) - which column in "dataset.info" contains the phenotype data?
ontology = ""             # (required, string)  - which dictionary data type to use? Possible values are "GO" and "other". For GO, KOMODO2 can compute normalization data.
dict.path = ""            # (required, string)  - file for dictionary file (two-column file containing annotation IDs and their descriptions, not needed for GO
column = ""               # (required, string)  - which column in annotation files should be used (column name)
denominator.column =      # (optional, numeric) - which column contains normalization data
tree.path = ""            # (required, string)  - path for tree file in either newick or nexus format
tree.type = ""            # (required, string)  - tree file type (either "nexus" or "newick", case-sensitive)
cores =                   # (optional, numeric) - how many cores to use?
linear.model.cutoff = 0.5 # (required, numeric) - Parameter that regulates how much graphical output is produced. We configure it to generate plots only for annotation terms with corrected q-values for phylogenetically independent contrasts smaller than 0.5.
