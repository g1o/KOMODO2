# This is a sample parameter file for KOMODO2.
# NOTICE: It assumes that the relevant data files (including this one) were
# retrieved using KOMODO2::retrieve_data_files("./data_folder").
# Change the target folder/file paths below if needed

annotation.files.dir = "./data_folder/domain2GO/" #Directory where annotation files are located. If not provided file described in variable "dataset.info" should contain absolute paths to annotation files

output.dir = "./results/domain2GO_Pan_proxy/" #output directory for results

dataset.info = "./data_folder/metadata/Pfam_metadata_Pan_proxy.txt" #genome metadata file it should contain at least for each genome: 1) path for annotation data; 2) phenotype data (numeric); 3) normalization data (numeric)

x.column = 2 #which column contains phenotype data?

ontology = "GO" #which dictionary data type?

dict.path = "" #file for dictionary file (two-column file containing annotation IDs and their descriptions not needed for GO

column = "Pfam" #which column in annotation file should be used (column name)

denominator.column = 4 #which column contains normalization data (numeric)

tree.path = "./data_folder/trees/tree_genome_IDs.nwk" #path for tree file

tree.type = "newick" #tree file type (either "nexus" or "newick" case-sensitive)

# cores = 4 #how many cores to use?

linear.model.cutoff = 0.5 #basically to tell KOMODO2 how much graphical output it should produce. We configure it to generate plots only for annotation terms with corrected q-values for phylogenetically independent contrasts smaller than 0.5.
