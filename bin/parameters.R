# parameters.R - default form of the parameters, possible to expand later.
# For now, just creates the empty list that guides the initial steps, just
# copy and paste the one with the type of analysis to use and fill the values.
# 
# Significance input:
#   type: (char) which comparison module to use:
#                - "significance": compares two groups of genomes within an 
#                                  ontology
#   test.path: (char) path of dir with annotation files of the test group.
#   back.path: (char) path of dir with annotation files of the background one.
#
#   ontology: (char) which ontology to use: "GO" (or "Gene Ontology"),
#                    "KEGG" and "other"
# 
#   dict.path: (char) file with the dictionary (terms and their meaning) 
#                     of an ontology, used if ontology is set as "other"
#   column: (char) name of the column with the functional annotation to be used
#                  in the analysis.
#   bootstrap: (numeric) number of bootstrap samples to run. Set it to
#                        any value below or equal to 1 to skip it.
#   criticalValue: (numeric) cutoff for the significance criteria in the 
#                            bootstrap counting. The terms must show a value 
#                            below it to be significant.
#
#
# Correlation input:
#   type: (char) which comparison module to use
#                - "correlation": establishes how much a variable explains the 
#                                 variations seen in the genomes
#   dataset.info: (char) path to the file with the paths of the dataset files
#                        (genomes and their annotations, y variable)
#                        and the genomes' attributes (x variable) columns.
#   x.column: (integer) number of the column to be used as the x variable.
#                       KOMODO2 assumes the first column has the paths to
#                       the annotation files, don't use 'x.column = 1'.
#   ontology: (char) which ontology to use: "GO" (or "Gene Ontology"),
#                    "KEGG" and "other"
# 
#   dict.path: (char) file with the dictionary (terms and their meaning) 
#                     of an ontology, used if ontology is set as "other"
#   column: (char) name of the column with the functional annotation to be used
#                  in the analysis (y variable).
#


KOMODO2 <- list(type = "correlation",
                output.dir = "../results/encephalization_Family",
                dataset.info = "../projects/encephalization/superfamily_info",
                x.column = 2,
                ontology = "other",
                dict.path = "../projects/encephalization/dict.tab",
                column = "Family ID",
                cores = 10)

KOMODO2 <- list(type = "correlation",
                output.dir = "../results/encephalization_Superfamily",
                dataset.info = "../projects/encephalization/superfamily_info",
                x.column = 2,
                ontology = "other",
                dict.path = "../projects/encephalization/dict.tab",
                column = "Superfamily ID",
                cores = 10)

KOMODO2 <- list(type = "correlation",
                output.dir = "../results/encephalization_InterPro",
                dataset.info = "../projects/encephalization/dataset_info",
                x.column = 2,
                ontology = "other",
                dict.path = "",
                column = "Cross-reference (InterPro)",
                cores = 10)

KOMODO2 <- list(type = "correlation",
                output.dir = "../results/encephalization_KO",
                dataset.info = "../projects/encephalization/dataset_info",
                x.column = 2,
                ontology = "KEGG",
                dict.path = "",
                column = "Cross-reference (KO)",
                cores = 10)


KOMODO2 <- list(type = "correlation",
                output.dir = "../results/encephalization_GO",
                dataset.info = "../projects/encephalization/dataset_info",
                x.column = 2,
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                cores = 10)

# KOMODO2 <- list(type = "correlation",
#                 output.dir = "../results/GO_encephalization",
#                 test.path = "../projects/encephalization/",
#                 back.path = "../projects/e_coli_data/nonEHEC_back",
#                 ontology = "GO",
#                 dict.path = "",
#                 column = "Gene ontology IDs",
#                 bootstrap = 100,
#                 criticalValue = 0.05,
#                 cores = 40)


KOMODO2 <- list(type = "significance",
                output.dir = "../results/GO_nonEHEC",
                test.path = "../projects/e_coli_data/EHEC_test",
                back.path = "../projects/e_coli_data/nonEHEC_back",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)
                
KOMODO2 <- list(type = "significance",
                output.dir = "../results/KO_nonEHEC",
                test.path = "../projects/e_coli_data/EHEC_test",
                back.path = "../projects/e_coli_data/nonEHEC_back",
                ontology = "KEGG",
                dict.path = "",
                column = "Cross-reference (KO)",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/GO_EHEC",
                test.path = "../projects/e_coli_data/EHEC_test",
                back.path = "../projects/e_coli_data/EHEC_back",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.01,
                cores = 40)
                
KOMODO2 <- list(type = "significance",
                output.dir = "../results/KO_EHEC",
                test.path = "../projects/e_coli_data/EHEC_test",
                back.path = "../projects/e_coli_data/EHEC_back",
                ontology = "KEGG",
                dict.path = "",
                column = "Cross-reference (KO)",
                bootstrap = 100,
                criticalValue = 0.01,
                cores = 40)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/plasmids_KO",
                test.path = "../projects/e_coli_data/plasmids_patho",
                back.path = "../projects/e_coli_data/plasmids_null",
                ontology = "KEGG",
                dict.path = "",
                column = "Cross-reference (KO)",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/plasmids_GO",
                test.path = "../projects/e_coli_data/plasmids_patho",
                back.path = "../projects/e_coli_data/plasmids_null",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/plasmids",
                test.path = "../projects/plasmids/patho",
                back.path = "../projects/plasmids/null",
                ontology = "other",
                dict.path = "",
                column = "Plasmid",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)


KOMODO2 <- list(type = "significance",
                output.dir = "../results/GO_ecoli_refactored",
                test.path = "../projects/e_coli_data/renamed_patho",
                back.path = "../projects/e_coli_data/renamed_null",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 2)


KOMODO2 <- list(type = "significance",
                output.dir = "../results/malaria_GO_with_gambiae",
                test.path = "../projects/malaria/a_gambiae_GO_test",
                back.path = "../projects/malaria/a_gambiae_GO_back",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/gambiae_complex_OG",
                test.path = "../projects/malaria/gambiae_OG_test",
                back.path = "../projects/malaria/gambiae_OG_back",
                ontology = "other",
                dict.path = "../projects/malaria/diptera_dict.tsv",
                column = "OG",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 48)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/gambiae_complex_GO",
                test.path = "../projects/malaria/gambiae_GO_test",
                back.path = "../projects/malaria/gambiae_GO_back",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 1,
                criticalValue = 0.05,
                cores = 40)
                
KOMODO2 <- list(type = "significance",
                output.dir = "../results/malaria_GO",
                test.path = "../projects/malaria/malaria_GO_test",
                back.path = "../projects/malaria/malaria_GO_back",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)
                
KOMODO2 <- list(type = "significance",
                output.dir = "../results/malaria_OG",
                test.path = "../projects/malaria/malaria_OG_test",
                back.path = "../projects/malaria/malaria_OG_back",
                ontology = "other",
                dict.path = "../projects/malaria/diptera_dict.tsv",
                column = "OG",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)
                
KOMODO2 <- list(type = "significance",
                output.dir = "../results/malaria_test",
                test.path = "../projects/malaria/anopheles_mRNA_path",
                back.path = "../projects/malaria/drosophila_path",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 40)

KOMODO2 <- list(type = "significance",
                output.dir = "../results/KO_ecoli",
                test.path = "../projects/e_coli_data/renamed_patho",
                back.path = "../projects/e_coli_data/renamed_null",
                ontology = "KEGG",
                dict.path = "",
                column = "Cross-reference (KO)",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 10)
                
KOMODO2 <- list(type = "significance",
                output.dir = "../results/GO_ecoli",
                test.path = "../projects/e_coli_data/renamed_patho",
                back.path = "../projects/e_coli_data/renamed_null",
                ontology = "GO",
                dict.path = "",
                column = "Gene ontology IDs",
                bootstrap = 100,
                criticalValue = 0.05,
                cores = 10)
# 
# KOMODO2 <- list(type = "significance",
#                 test.path = "../projects/e_coli_data/patho",
#                 back.path = "../projects/e_coli_data/null",
#                 ontology = "GO",
#                 dict.path = "",
#                 column = "Gene ontology IDs",
#                 bootstrap = 1,
#                 criticalValue = 0.05,
#                 cores = 40)
# 
# KOMODO2 <- list(type = "significance",
#                 test.path = "../saccharomyces_data/dome_wine",
#                 back.path = "../saccharomyces_data/null_wine",
#                 ontology = "other",
#                 dict.path = "",
#                 column = "Cross-reference (ORTHODB)",
#                 bootstrap = 0,
#                 criticalValue = 0.05)
# 
# KOMODO2 <- list(type = "correlation",
#                 x.path = "../pexpansion/cell_types.tsv",
#                 y.path = "../pexpansion/scop_data.tsv",
#                 ontology = "other",
#                 dict.path = "../pexpansion/dictionary_old_scop.tsv")
# 
# 
# KOMODO2 <- list(type = "significance",
#                 test.path = "../projects/Thiago_Mafra/test",
#                 back.path = "../projects/Thiago_Mafra/back",
#                 ontology = "GO",
#                 dict.path = "",
#                 column = "Cross-reference (ORTHODB)",
#                 bootstrap = 0,
#                 criticalValue = 0.05)
# 
# KOMODO2 <- list(type = "significance",
#                 test.path = "../projects/e_coli_data/EHEC_subgroup",
#                 back.path = "../projects/e_coli_data/null_no_K12",
#                 ontology = "GO",
#                 dict.path = "",
#                 column = "Gene ontology IDs",
#                 bootstrap = 100,
#                 criticalValue = 0.05)
# 
# KOMODO2 <- list(type = "significance",
#                 test.path = "../projects/e_coli_data/urinary_ecoli",
#                 back.path = "../projects/e_coli_data/null_no_K12",
#                 ontology = "GO",
#                 dict.path = "",
#                 column = "Gene ontology IDs",
#                 bootstrap = 0,
#                 criticalValue = 0.05)
# 
# KOMODO2 <- list(type = "significance",
#                 test.path = "../projects/e_coli_data/second_subgroup",
#                 back.path = "../projects/e_coli_data/null_no_K12",
#                 ontology = "GO",
#                 dict.path = "",
#                 column = "Gene ontology IDs",
#                 bootstrap = 100,
#                 criticalValue = 0.05)
