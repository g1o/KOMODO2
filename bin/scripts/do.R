# do.R - last step of the LCFD workflow of KOMODO2.
# Performs the actual analysis of the program and produces all data and tables
#  desired. Sources the "func.R" file to access all necessary functions.
#
# input:
#   test.anno: (list) list of treated data frames, each one with the
#                     annotation of a genome of the test group.
#   back.anno: (list) list of treated data frames, each one with the
#                     annotation of a genome of the background group.
# output:
# 


# Loading functions without filling the workspace
KOMODO2.env <- new.env()
suppressMessages(sys.source("func.R", envir = KOMODO2.env))
suppressMessages(attach(KOMODO2.env))
rm(KOMODO2.env)

# Prepare the Bioconductor and databases
LoadBioconductor()

# Prepare parallelization
if (!is.null(KOMODO2$cores)) {
  if (KOMODO2$cores > 1) {
    KOMODO2$cl <- makeCluster(KOMODO2$cores)
    registerDoParallel(KOMODO2$cl)
  }
}


if (KOMODO2$type == "correlation") {
  tree <- read.tree("/home/chico/projects/KOMODO2/validation/Cetartiodactyla_weight/data/trees/Cetacean_tree_genome_ids.nwk")
  tree<-multi2di(tree)  
if (tolower(KOMODO2$ontology) == "go" | 
      tolower(KOMODO2$ontology) == "gene ontology") {
    KOMODO2$allAncestor <- ListAncestors()
    KOMODO2$allObsolete <- ListObsoletes()
    KOMODO2$allSynonym <- ListSynonyms()
  } else if (tolower(KOMODO2$ontology) == "kegg") {
    KOMODO2$allKOs <- ListKOs()
  } else if (tolower(KOMODO2$ontology) == "other" & KOMODO2$dict.path == "") {
    KOMODO2$dictionary <- CreateDictionary(KOMODO2$y.anno)
  }

  if (is.null(KOMODO2$yElementCount)) {
    KOMODO2$yElementCount <- GroupElementCount(KOMODO2$y.anno)
  }
  
  KOMODO2$y <- AddGenomeVectors(KOMODO2$y.anno, KOMODO2$y.name) 

  print ("Computing sum of annotation elements...")

  KOMODO2$sum <- lapply(KOMODO2$y, sum)

  print("Done")
 
  print ("Computing standard deviation of annotation elements...")

  KOMODO2$sd <- lapply(KOMODO2$y, sd)

  print ("Done")
 
  print("Computing contrasts...") 
  tmp <- FindContrasts(KOMODO2$x, KOMODO2$y,
                                  tree, "pic",
                                  KOMODO2$denominator)  
  print("Done");
  KOMODO2$contrasts <- tmp
#  KOMODO2$contrasts <- tmp$corr_values
#  KOMODO2$contrast_models <- tmp$data

  print("Computing correlations...")

  tmp <- FindCorrelations(KOMODO2$x, KOMODO2$y, 
                          "pearson",
                          KOMODO2$denominator)

  KOMODO2$correlations.pearson <- tmp$cor

  KOMODO2$correlations.pvalue.pearson <- tmp$cor.pvalue

  tmp <- FindCorrelations(KOMODO2$x, KOMODO2$y,
                          "spearman",
                          KOMODO2$denominator)

  KOMODO2$correlations.spearman <- tmp$cor

  KOMODO2$correlations.pvalue.spearman <- tmp$cor.pvalue

  tmp <- FindCorrelations(KOMODO2$x, KOMODO2$y,
                          "kendall",
                          KOMODO2$denominator)

  KOMODO2$correlations.kendall <- tmp$cor

  KOMODO2$correlations.pvalue.kendall <- tmp$cor.pvalue

  print("Done")

#  KOMODO2$correlations.spearman <- FindCorrelations(KOMODO2$x, KOMODO2$y,
#                                                    "spearman",
#                                                    KOMODO2$denominator)
#  KOMODO2$correlations.kendall <- FindCorrelations(KOMODO2$x, KOMODO2$y,
#                                                   "kendall",
#                                                   KOMODO2$denominator)

  KOMODO2$annotation.cor <- AnnotateResults(KOMODO2$correlations.pearson,
                                            KOMODO2$ontology)

  print("Printing results (1)...")

  PrintCResults(KOMODO2$correlations.pearson, KOMODO2$annotation.cor, 
                "p_corr_results.tsv", KOMODO2$ontology)
  PrintCResults(KOMODO2$correlations.spearman, KOMODO2$annotation.cor, 
                "s_corr_results.tsv", KOMODO2$ontology)
  PrintCResults(KOMODO2$correlations.kendall, KOMODO2$annotation.cor, 
                "k_corr_results.tsv", KOMODO2$ontology)

  KOMODO2$results.correlations.pvalue.pearson <- MultipleHypothesisCorrection(KOMODO2$correlations.pvalue.pearson)
  KOMODO2$results.correlations.pvalue.spearman <- MultipleHypothesisCorrection(KOMODO2$correlations.pvalue.spearman)
  KOMODO2$results.correlations.pvalue.kendall <- MultipleHypothesisCorrection(KOMODO2$correlations.pvalue.kendall)
  KOMODO2$contrasts.corrected <- MultipleHypothesisCorrection(KOMODO2$contrasts)

  print("Printing results (2)...")

  KOMODO2$contrasts.cor <- AnnotateResults(KOMODO2$contrasts,
                                            KOMODO2$ontology)

  KOMODO2$sum.cor <- AnnotateResults(KOMODO2$sum, KOMODO2$ontology)

  KOMODO2$sd.cor <- AnnotateResults(KOMODO2$sd, KOMODO2$ontology)

  PrintCResults(KOMODO2$contrasts.corrected, KOMODO2$contrasts.cor,
                "contrasts_corrected.tsv", KOMODO2$ontology)
  
  PrintCResults(KOMODO2$contrasts, KOMODO2$contrasts.cor,
                "contrasts_raw.tsv", KOMODO2$ontology)

 
  PrintCResults(KOMODO2$results.correlations.pvalue.pearson, KOMODO2$annotation.cor,
                "p_corr_qvalues_results.tsv", KOMODO2$ontology)
  PrintCResults(KOMODO2$results.correlations.pvalue.spearman, KOMODO2$annotation.cor,
                "s_corr_qvalues_corr_results.tsv", KOMODO2$ontology)
  PrintCResults(KOMODO2$results.correlations.pvalue.kendall, KOMODO2$annotation.cor,
                "k_corr_qvalues_corr_results.tsv", KOMODO2$ontology)

  PrintCResults(KOMODO2$sum, KOMODO2$sum.cor,
                "sum.tsv", KOMODO2$ontology)
  
  PrintCResults(KOMODO2$sd, KOMODO2$sd.cor,
                "sd.tsv", KOMODO2$ontology)

#  common <- Reduce(intersect, list(names(KOMODO2$sd),names(KOMODO2$correlations.kendall), names(KOMODO2$correlations.pearson) ,names(KOMODO2$contrasts.corrected), names(KOMODO2$correlations.spearman), names(KOMODO2$sum)))
#  KOMODO2$sum <- subset(KOMODO2$sum, names(KOMODO2$sum) %in% common)
#  KOMODO2$sd <- subset(KOMODO2$sd, names(KOMODO2$sd) %in% common)
#  KOMODO2$contrasts.cor <- subset(KOMODO2$contrasts.cor, names(KOMODO2$contrasts.cor) %in% common)
#  KOMODO2$contrasts.corrected <- subset(KOMODO2$contrasts.corrected, names(KOMODO2$contrasts.corrected) %in% common)
#  KOMODO2$correlations.kendall <- subset(KOMODO2$correlations.kendall, names(KOMODO2$correlations.kendall) %in% common)
#  KOMODO2$correlations.spearman <- subset(KOMODO2$correlations.spearman, names(KOMODO2$correlations.spearman) %in% common)
#  KOMODO2$correlations.pearson <- subset(KOMODO2$correlations.pearson, names(KOMODO2$correlations.pearson) %in% common)
#  KOMODO2$correlations.pearson <- subset(KOMODO2$correlations.pearson, names(KOMODO2$correlations.pearson) %in% common)
#  KOMODO2$results.correlations.pvalue.pearson <- subset(KOMODO2$results.correlations.pvalue.pearson, names(KOMODO2$results.correlations.pvalue.pearson) %in% common)
#  KOMODO2$results.correlations.pvalue.spearman <- subset(KOMODO2$results.correlations.pvalue.spearman, names(KOMODO2$results.correlations.pvalue.spearman) %in% common)
#  KOMODO2$results.correlations.pvalue.kendall <- subset(KOMODO2$results.correlations.pvalue.kendall, names(KOMODO2$results.correlations.pvalue.kendall) %in% common)
#  KOMODO2$annotation.cor <- subset(KOMODO2$annotation.cor, names(KOMODO2$annotation.cor) %in% common)

#  KOMODO2$correlations.pearson <- subset(KOMODO2$correlations.pearson, names(KOMODO2$correlations.pearson) %in% common)

#  KOMODO2$correlations.kendall <- subset(KOMODO2$correlations.kendall, names(KOMODO2$correlations.kendall) %in% common)
#  KOMODO2$correlations.spearman <- subset(KOMODO2$correlations.spearman, names(KOMODO2$correlations.spearman) %in% common)

#  KOMODO2$annotation.cor <- subset(KOMODO2$annotation.cor, names(KOMODO2$annotation.cor) %in% common)
#  KOMODO2$annotation.cor <- subset(KOMODO2$annotation.cor, names(KOMODO2$annotation.cor) %in% common)
#  KOMODO2$correlations.pearson
#  KOMODO2$y <- KOMODO2$y[names(KOMODO2$y) %in% common]
cutoff=0.2;
sumY<-sapply(KOMODO2$y,sum) # done as vector, it is a simply sum and it is fast as this, must check how the list type is used before this one.
sumY<-sumY[!sumY==0] # filter out those with 0 counts
Y<-KOMODO2$y[,colSums(KOMODO2$y)!=0]

plotframe<-rbind(KOMODO2$contrasts.corrected[order(names(KOMODO2$contrasts.corrected))],
          KOMODO2$results.correlations.pvalue.pearson[order(names(KOMODO2$results.correlations.pvalue.pearson))],
          sumY[order(names(sumY))],
          KOMODO2$results.correlations.pvalue.spearman[order(names(KOMODO2$results.correlations.pvalue.spearman))],
          KOMODO2$results.correlations.pvalue.kendall[order(names(KOMODO2$results.correlations.pvalue.kendall))],
	  KOMODO2$sd[order(names(KOMODO2$sd))] ,
          Y[,order(colnames(Y))] )  

rownames(plotframe)[1:6]<-c("corrected_contrasts",
                     "PearsonCorrelation",
                     "size",
                     "SpearmanCorrelation",
                     "KendallCorrelation",
                     "sd")
			

plotframe<-as.data.frame(t(plotframe))
description<-unlist(KOMODO2$annotation.cor)
description<-description[order(names(description))]
plotframe$description<-description
plotframe$name<-rownames(plotframe)
df_cutoff<-plotframe[plotframe$corrected_contrasts<cutoff,]
df_cutoff<-df_cutoff[df_cutoff$sd!=0,] #removing trivial cases, constant values. 

wd<-normalizePath(KOMODO2$output.dir);
render("KOMODO2_correlation_report.Rmd",output_file=paste0(wd,'/KOMODO2_report.html') )
print("Done")


} else if (KOMODO2$type == "significance") {
  if (tolower(KOMODO2$ontology) == "go" | 
      tolower(KOMODO2$ontology) == "gene ontology") {
    KOMODO2$allAncestor <- ListAncestors()
    KOMODO2$allObsolete <- ListObsoletes()
    KOMODO2$allSynonym <- ListSynonyms()
  } else if (tolower(KOMODO2$ontology) == "kegg") {
    KOMODO2$allKOs <- ListKOs()
  } else if (tolower(KOMODO2$ontology) == "other" & KOMODO2$dict.path == "") {
    KOMODO2$dictionary <- CreateDictionary(KOMODO2$test.anno,
                                           KOMODO2$back.anno)
  }
  
  KOMODO2$testGV <- AddGenomeVectors(KOMODO2$test.anno, KOMODO2$test.name)
  KOMODO2$backGV <- AddGenomeVectors(KOMODO2$back.anno, KOMODO2$back.name)
  
  
  # Preparing a listing of parameter vectors
  if (is.null(KOMODO2$testElementCount)) {
    KOMODO2$testElementCount <- GroupElementCount(KOMODO2$test.anno,
                                                  KOMODO2$test.name)
  }
  if (is.null(KOMODO2$backElementCount)) {
    KOMODO2$backElementCount <- GroupElementCount(KOMODO2$back.anno,
                                                  KOMODO2$back.name)
  }
  KOMODO2$parameterVectors <- ParameterVectors(KOMODO2$testGV, KOMODO2$backGV,
                                               KOMODO2$testElementCount,
                                               KOMODO2$backElementCount)

  # Statistical test
  KOMODO2$pvalue.f.over <- StatisticalTest(KOMODO2$parameterVectors,
                                           "fisher.over")
  KOMODO2$pvalue.f.under <- StatisticalTest(KOMODO2$parameterVectors,
                                            "fisher.under")
  KOMODO2$pvalue.ks.over <- CompareDistributions(KOMODO2$parameterVectors,
                                                 KOMODO2$testGV,
                                                 KOMODO2$backGV)
  KOMODO2$pvalue.ks.under <- CompareDistributions(KOMODO2$parameterVectors,
                                                  KOMODO2$testGV, 
                                                  KOMODO2$backGV, "ks.under")
  KOMODO2$pvalue.w.over <- CompareDistributions(KOMODO2$parameterVectors,
                                                KOMODO2$testGV, KOMODO2$backGV,
                                                "wilcox.over")
  KOMODO2$pvalue.w.under <- CompareDistributions(KOMODO2$parameterVectors,
                                                KOMODO2$testGV, KOMODO2$backGV,
                                                "wilcox.under")


  # Multiple Hypothesis Correction
  KOMODO2$results.f.over <- MultipleHypothesisCorrection(KOMODO2$pvalue.f.over)
  KOMODO2$results.f.under <- MultipleHypothesisCorrection(
                                                       KOMODO2$pvalue.f.under)
  
  KOMODO2$results.ks.over <- MultipleHypothesisCorrection(
                                                       KOMODO2$pvalue.ks.over)
  KOMODO2$results.ks.under <- MultipleHypothesisCorrection(
                                                       KOMODO2$pvalue.ks.under)
  KOMODO2$results.w.over <- MultipleHypothesisCorrection(
                                                       KOMODO2$pvalue.w.over)
  KOMODO2$results.w.under <- MultipleHypothesisCorrection(
                                                       KOMODO2$pvalue.w.under)


  # Effect size (Phi Coefficient)
  KOMODO2$phi <- GetPhiCoefficient(KOMODO2$parameterVectors)
  KOMODO2$r <- Get.r.Coefficient(KOMODO2$testGV, KOMODO2$backGV,
                                         KOMODO2$pvalue.ks.over)

  # Confidence Intervals
#   KOMODO2$confInt.phi <- ConfidenceInterval(KOMODO2$phi, 1.96, "P",
#                                             KOMODO2$parameterVectors)
  KOMODO2$confInt.r <- ConfidenceInterval(KOMODO2$r, 1.96, "N", 
                                          KOMODO2$test.name,
                                          KOMODO2$back.name)

  # Fold change
  KOMODO2$foldChange.p <- CalculateFoldChange(KOMODO2$parameterVectors)
  KOMODO2$foldChange.n <- CalculateFoldChange(KOMODO2$parameterVectors, "N", 
                                              KOMODO2$test.name,
                                              KOMODO2$back.name)

  # Annotating results
  KOMODO2$annotation <- AnnotateResults(KOMODO2$results.f.over,
                                        KOMODO2$ontology)

  # Bootstrap analysis
  if (KOMODO2$bootstrap > 1) {
    KOMODO2$bootIndices <- BootGenomeVectorsIndices(KOMODO2$testGV,
                                                    KOMODO2$backGV,
                                                    KOMODO2$bootstrap)

    KOMODO2$bParameterVectors <- BootParameterVectors(KOMODO2$testGV,
                                                    KOMODO2$backGV,
                                                    KOMODO2$testElementCount,
                                                    KOMODO2$backElementCount,
                                                    KOMODO2$bootIndices)
                                                  
    KOMODO2$count.f.over <- BootStatisticalTest(KOMODO2$bParameterVectors,
                                                "fisher.over")
    KOMODO2$count.f.under <- BootStatisticalTest(KOMODO2$bParameterVectors,
                                                "fisher.under")
    KOMODO2$count.ks.over <- BootDistributions(KOMODO2$bParameterVectors, 
                                                KOMODO2$bootIndices,
                                                KOMODO2$testGV, KOMODO2$backGV)
    KOMODO2$count.ks.under <- BootDistributions(KOMODO2$bParameterVectors, 
                                                KOMODO2$bootIndices,
                                                KOMODO2$testGV, KOMODO2$backGV,
                                                "ks.under")
    KOMODO2$count.w.over <- BootDistributions(KOMODO2$bParameterVectors, 
                                               KOMODO2$bootIndices,
                                               KOMODO2$testGV, KOMODO2$backGV,
                                               "wilcox.over")
    KOMODO2$count.w.under <- BootDistributions(KOMODO2$bParameterVectors, 
                                                KOMODO2$bootIndices,
                                                KOMODO2$testGV, KOMODO2$backGV,
                                                "wilcox.under")
                                                  
    KOMODO2$hCount.f.over <- Boot_MHCorrection(KOMODO2$count.f.over)
    KOMODO2$hCount.f.under <- Boot_MHCorrection(KOMODO2$count.f.under)
                                               
    KOMODO2$hCount.ks.over <- Boot_MHCorrection(KOMODO2$count.ks.over)
    KOMODO2$hCount.ks.under <- Boot_MHCorrection(KOMODO2$count.ks.under)
    KOMODO2$hCount.w.over <- Boot_MHCorrection(KOMODO2$count.w.over)
    KOMODO2$hCount.w.under <- Boot_MHCorrection(KOMODO2$count.w.under)


    KOMODO2$count.f.over <- SignificantSamples(KOMODO2$hCount.f.over,
                                               KOMODO2$parameterVectors,
                                               KOMODO2$criticalValue)
    KOMODO2$count.f.under <- SignificantSamples(KOMODO2$hCount.f.under,
                                                KOMODO2$parameterVectors,
                                                KOMODO2$criticalValue)
    KOMODO2$count.ks.over <- SignificantSamples(KOMODO2$hCount.ks.over,
                                                KOMODO2$parameterVectors,
                                                KOMODO2$criticalValue)
    KOMODO2$count.ks.under <- SignificantSamples(KOMODO2$hCount.ks.under,
                                                 KOMODO2$parameterVectors,
                                                 KOMODO2$criticalValue)
    KOMODO2$count.w.over <- SignificantSamples(KOMODO2$hCount.w.over,
                                               KOMODO2$parameterVectors,
                                               KOMODO2$criticalValue)
    KOMODO2$count.w.under <- SignificantSamples(KOMODO2$hCount.w.under,
                                                KOMODO2$parameterVectors,
                                                KOMODO2$criticalValue)
  }
  
  
  

  # Printing results into a file
  PrintResults(KOMODO2, "f.over", "fisher_over")
  PrintResults(KOMODO2, "f.under", "fisher_under")

  PrintNPResults(KOMODO2, "ks.over", "ks_over")
  PrintNPResults(KOMODO2, "ks.under", "ks_under")
  PrintNPResults(KOMODO2, "w.over", "wilcox_over")
  PrintNPResults(KOMODO2, "w.under", "wilcox_under")



  PrintVolcanoPlot(KOMODO2$foldChange.p, KOMODO2$results.f.over,
                   outputName = "volcano_fover.png")
  PrintVolcanoPlot(KOMODO2$foldChange.p, KOMODO2$results.f.under,
                   outputName = "volcano_funder.png")
  PrintVolcanoPlot(KOMODO2$foldChange.n, KOMODO2$results.ks.over,
                   outputName = "volcano_ksover.png")
  PrintVolcanoPlot(KOMODO2$foldChange.n, KOMODO2$results.ks.under,
                   outputName = "volcano_ksunder.png")
  PrintVolcanoPlot(KOMODO2$foldChange.n, KOMODO2$results.w.over,
                   outputName = "volcano_wover.png")
  PrintVolcanoPlot(KOMODO2$foldChange.n, KOMODO2$results.w.under,
                   outputName = "volcano_wunder.png")

  PrintSignificantHist(KOMODO2$results.f.over, 0.05, KOMODO2$testGV, 
                       KOMODO2$backGV, KOMODO2$annotation)

  PrintSignTermCor(KOMODO2$ontology, KOMODO2$results.f.over, 
                   KOMODO2$annotation, "pearson", KOMODO2$testGV,
                   KOMODO2$backGV, 0.05)

  # Printing the heatmap with the hierarchical clustering
  HierarchicalClustering(KOMODO2, test = "f.over",
                         outputName = "heatmap_fover.png",
                         criticalValue = 0.05, margins = c(50,30),
                         nResults = 25)
  HierarchicalClustering(KOMODO2, test = "f.under",
                         outputName = "heatmap_funder.png",
                         criticalValue = 0.05, margins = c(50,30),
                         nResults = 25)

  HierarchicalClustering(KOMODO2, test = "ks.over",
                         outputName = "heatmap_ksover.png",
                         criticalValue = 0.05, margins = c(50,30),
                         nResults = 25)
  HierarchicalClustering(KOMODO2, test = "ks.under",
                         outputName = "heatmap_ksunder.png",
                         criticalValue = 0.05, margins = c(50,30),
                         nResults = 25)
  HierarchicalClustering(KOMODO2, test = "w.over",
                         outputName = "heatmap_wover.png",
                         criticalValue = 0.05, margins = c(50,30),
                         nResults = 25)
  HierarchicalClustering(KOMODO2, test = "w.under",
                         outputName = "heatmap_wunder.png",
                         criticalValue = 0.05, margins = c(50,30),
                         nResults = 25)


}

# Close parallel structure
if (!is.null(KOMODO2$cores)) {
  if (KOMODO2$cores > 1) {
    stopCluster(KOMODO2$cl)
  }
}
