#This is a sample parameter file for KOMODO2. Most variables are self-explanatory, and a formal definition of each can be found in KOMODO2 documentation.

KOMODO2 <- list(output.dir = "/Users/chico/projects/KOMODO2/development/KOMODO2/results/EUK/GO/", #output directory for results

                dataset.info = "/Users/chico/projects/KOMODO2/development/KOMODO2/data/validation/EUK/metadata/GO_metadata_KOMODO2.txt", #genome metadata file, it should contain at least, for each genome: 1) path for annotation data; 2) phenotype data (numeric); 3) normalization data (numeric)

                x.column = 2, #which column contains phenotype data?

                ontology = "GO", #which dictionary data type?

                dict.path = "", #file for dictionary file (two-column file containing annotation IDs and their descriptions, not needed for GO

                column = "GO", #which column in annotation file should be used (column name)
                denominator.column = 3, #which column contains normalization data (numeric)
                tree_path = "/Users/chico/projects/KOMODO2/development/KOMODO2/data/validation/EUK/tree/eukarya_tree.nex", #path for tree file

                tree_type = "nexus", #tree file type (either "nexus" or "newick", case-sensitive)
                cores = 4 #how many cores to use?


                #ADVANCED PARAMETERS (not really...)

                #These parameters are basically to tell KOMODO2 how much graphical output it should produce. We configure it to generate plots only for annotation terms with corrected q-values for phylogenetically independent contrasts greater than 0.5. 

                )
