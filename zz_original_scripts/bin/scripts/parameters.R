#This is a sample parameter file for KOMODO2. Most variables are self-explanatory, and a formal definition of each can be found in KOMODO2 documentation.

KOMODO2 <- list(annotation_files_dir = "", #(mandatory, string, use quotes) - Directory where annotation files are located. If not provided, file described in variable "dataset.info" should contain absolute paths to annotation files
                output.dir = "", #(mandatory, string, use quotes) - output directory for results

                dataset.info = "", #(mandatory, string, use quotes) genome metadata file, it should contain at least:
                                   # 1) path for annotation data (if "annotation_files_dir" not provided OR file names (if "annotation_files_dir) is provided. Please notice this information should be the first column in metadata file;
                                   # 2) phenotype data (numeric, this is the value KOMODO2 uses to rank species when searching for associations)
                                   # 3) normalization data (numeric, this is the value KOMODO2 uses as a denominator to compute annotation term frequencies to remove potential biases caused by, for instance, overannotation of model organisms or large differences in the counts of genomic elements). Please notice KOMODO2 does not require normalization data for GO, as it computes the total number of GO terms per species and uses it as a normalizing factor.

                x.column = , #(mandatory, numeric, do *not* use quotes) - which column in "dataset.info" contains the phenotype data?

                ontology = "", #(mandatory, string, use quotes) - which dictionary data type to use? Possible values are "GO" and "other". For GO, KOMODO2 can compute normalization data. 

                dict.path = "", #(mandatory, string, use quotes) - file for dictionary file (two-column file containing annotation IDs and their descriptions, not needed for GO

                column = "", #(string, mandatory, use quotes) - which column in annotation files should be used (column name)
                denominator.column = , #(numeric, optional for GO, do *not* use quotes) which column contains normalization data
                tree_path = "", #(string, mandatory, use quotes) - path for tree file in either newick or nexus format

                tree_type = "", #(string, mandatory, use quotes) - tree file type (either "nexus" or "newick", case-sensitive)
                cores = , #(numeric, mandatory, do *not* use quotes) - how many cores to use?


                #ADVANCED PARAMETERS (not really...)

                #These parameters are basically to tell KOMODO2 how much graphical output it should produce. We configure it to generate plots only for annotation terms with corrected q-values for phylogenetically independent contrasts smaller than 0.5. 
                
                linear_model_cutoff = 0.5

                )
