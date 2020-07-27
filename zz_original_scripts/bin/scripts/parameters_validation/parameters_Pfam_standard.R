#This is a sample parameter file for KOMODO2. Most variables are self-explanatory, and a formal definition of each can be found in KOMODO2 documentation.

KOMODO2 <- list(annotation_files_dir = "../../data/Pfam/", #Directory where annotation files are located. If not provided, file described in variable "dataset.info" should contain absolute paths to annotation files
                output.dir = "../../results/Pfam_standard/", #output directory for results

                dataset.info = "../../data/metadata/Pfam_metadata_standard.txt", #genome metadata file, it should contain at least, for each genome: 1) path for annotation data; 2) phenotype data (numeric); 3) normalization data (numeric)

                x.column = 2, #which column contains phenotype data?

                ontology = "other", #which dictionary data type?

                dict.path = "../validation/data/dics/Pfam.dic", #file for dictionary file (two-column file containing annotation IDs and their descriptions, not needed for GO

                column = "Pfam", #which column in annotation file should be used (column name)
                denominator.column = 4, #which column contains normalization data (numeric)
                tree_path = "../validation/data/trees/tree_genome_IDs.nwk", #path for tree file

                tree_type = "newick", #tree file type (either "nexus" or "newick", case-sensitive)
                cores = 4, #how many cores to use?


                #ADVANCED PARAMETERS (not really...)

                #These parameters are basically to tell KOMODO2 how much graphical output it should produce. We configure it to generate plots only for annotation terms with corrected q-values for phylogenetically independent contrasts smaller than 0.5. 
                
                linear_model_cutoff = 0.5

                )
