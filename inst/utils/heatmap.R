#loading libraries
library(RColorBrewer)
library(gplots)
library(preprocessCore)
library(extrafont)
library(ape)
library(dendextend)
library(heatmaply)
loadfonts()

tmp <- read.table("data_folder/metadata/Pfam_metadata_Pan_proxy.txt", sep="\t")

names <- tmp[, c(1, 7,8)]

names <- as.data.frame(output$y.name, output$groups, output$short.name)

colnames(names) <- c("genomeID", "group", "shortName")

rm(tmp)

tmp <- read.table("../../docs/manuscript/biorx/datasets/domain_IDs.txt")
#tmp <- read.table("../../docs/manuscript/biorx/datasets/gene2GO_IDs.txt")
#tmp <- read.table("../../docs/manuscript/biorx/datasets/domain2GO_IDs.txt")


ids <- as.vector(tmp$V1)
rm(tmp)

#tmp <- KOMODO2$y
#tmp <- KOMODO2$y[KOMODO2$sum > 0]
tmp <- output$y[as.vector(ids)]

norm <- tmp/output$denominator
#norm <- tmp
rm(tmp)

rownames(norm) <- names$shortName[match(rownames(norm), names$genomeID)]
colnames(norm) <- paste0(ids, "\n", output$annotation.cor[ids])

#rownames(norm) <- names[,3]

tree <- output$tree

#getting a good tree

tree$tip.label<-names[[3]][match(tree$tip.label, names[[1]])]

tree$tip.label <- sapply(tree$tip.label, function(x) as.character(x))

#tree$tip.label <- sapply(tree$tip.label, function(x) gsub(" ", "_", x))

#plot(tree)

tree2 <- as.dendrogram(as.hclust.phylo(tree))

clade_order <- order.dendrogram(tree2)

#ord.mat <- norm[clade_order$tip.name,clade_order$tip.name]

max <- ncol(norm)

mat <- as.matrix(norm[,1:max])



ord.mat <- mat[tree$tip.label,]

#z_score <- scale(norm[,1:99], center = TRUE, scale = FALSE)

my_palette <- colorRampPalette(c("white","blue"))(n = 51)


(jColors <-
   with(names,
        data.frame(LABEL = levels(group),
                   COLOR = I(brewer.pal(nlevels(group), name = 'Set1')))))


tmp <- list()

tmp$species <- rownames(ord.mat)

tmp$group <- names$group[match(tmp$species, names$shortName)]

tmp$COLOR <- jColors$COLOR[match(tmp$group, jColors$LABEL)]

species2color <- unique(as.data.frame(tmp)[,c("species","COLOR")])

rm(tmp)


#quantile normalization

tmp <- normalize.quantiles(as.matrix(ord.mat))
colnames(tmp) <- colnames(ord.mat)
rownames(tmp) <- rownames(ord.mat)

ord.norm.mat <- tmp
rm(tmp)

#lala <- list()
#lala$species <- rownames(ord.mat)
#lala$group <-names$V7[match(lala$species, names$V8)]
#lala$COLOR <- species2color$COLOR[match(lala$species, species2color$species)]

#distance1 = dist(as.matrix(norm[,1:max]), method = "euclidean")
#distance2 = dist(as.matrix(t(ord.mat)), method = "euclidean")

#distance1 = dist(as.matrix(z_score), method = "euclidean")
#distance2 = dist(as.matrix(t(z_score)), method = "euclidean")

#cluster1 = hclust(distance1, method = "ward.D2")
#cluster2 = hclust(distance2, method = "ward.D2")


#pdf using phylogeny to cluster data


pdf("teste.pdf", width = 120, height = 40)
distance2 = dist(as.matrix(t(ord.mat)), method = "euclidean")
#cluster2 = hclust(distance2, method = "ward.D2")

cluster2_final = set(as.dendrogram(hclust(distance2, method=c("average"))), "branches_lwd", 2)

tree_final <- set(tree2, "branches_lwd", 2)
#par(oma=c(10,15,10,10))
#lhei=c(2, 2, 2)
#heatmap using phylogeny to cluster genomes
heatmap.2(as.matrix(ord.mat),
          density.info="histogram",
          RowSideColors=species2color$COLOR,
          col = my_palette,
          cexRow=4,
          cexCol=2,
          margins =c(35,12),
          srtCol=45,
          Rowv=(tree_final),
          Colv=(cluster2_final),
          labRow=species2color$species,
          labCol = FALSE,
          trace=c("none"),
          font <- par(family="Courier New"),
#          lmat=rbind(c(5,4), c(3,2), c(0,1)),
#          lhei=c(2,4,0.2)
          )#, labCol = FALSE)
dev.off()


#heatmap using distance to cluster genomes

pdf("teste2.pdf", width = 120, height = 40)

distance1 = dist(as.matrix(ord.mat), method = "euclidean")
distance2 = dist(as.matrix(t(ord.mat)), method = "euclidean")

#distance1 = dist(as.matrix(z_score), method = "euclidean")
#distance2 = dist(as.matrix(t(z_score)), method = "euclidean")

cluster1_final = set(as.dendrogram(hclust(distance1, method=c("average"))), "branches_lwd", 5)
cluster2_final = set(as.dendrogram(hclust(distance2, method=c("average"))), "branches_lwd", 5)
#cluster1 = hclust(distance1, method = "average")
#cluster2 = hclust(distance2, method = "average")



heatmap.2(as.matrix(ord.mat),
          density.info="histogram",
          RowSideColors=species2color$COLOR,
          col = my_palette,
          cexRow=4,
          cexCol=2,
          margins =c(35,12),
          srtCol=45,
          Rowv=(cluster1_final),
          Colv=(cluster2_final),
#          Rowv=as.dendrogram(cluster1),
#          Colv=as.dendrogram(cluster2),
          labRow=species2color$species,
          labCol = FALSE,
          trace=c("none"),
          font <- par(family="Courier New")
          )

dev.off()


pdf("teste2.pdf", width = 20, height = 10)

#par(oma=c(1,15,1,1))

#heatmap using phylogeny to cluster genomes
#heatmap.2(as.matrix(ord.mat), density.info="histogram", RowSideColors=species2color$COLOR, col = my_palette, cexRow=1.5, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(tree2), Colv=as.dendrogram(cluster2), labRow=species2color$species, trace=c("none"), font <- par(family="Courier New"))

#heatmap using distance to cluster genomes
heatmap.2(as.matrix(ord.mat), density.info="histogram", RowSideColors=species2color$COLOR, col = my_palette, cexRow=1.2, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(cluster1), Colv=as.dendrogram(cluster2), labRow=species2color$species, trace=c("none"), font <- par(family="Courier New"))

dev.off()

heatmap.2(as.matrix(ord.mat), density.info="histogram", RowSideColors=species2color$COLOR, col = my_palette, cexRow=1.5, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(cluster1), Colv=as.dendrogram(cluster2), labRow=species2color$species, trace=c("none"), labCol = FALSE)

#heatmap.2(as.matrix(norm[,1:max]), density.info="histogram", RowSideColors=species2color$COLOR, col = my_palette, cexRow=1.6, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(cluster1), Colv=as.dendrogram(cluster2), labRow=species2color$species, trace=c("none"))
#heatmap.2(as.matrix(ord.mat[,1:max]), density.info="histogram", RowSideColors=species2color$COLOR, col = my_palette, cexRow=1.6, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(tree2), Colv=as.dendrogram(cluster2), labRow=species2color$species, trace=c("none"))
#heatmap.2(ord.mat, density.info="histogram", RowSideColors=lala$COLOR, col = my_palette, cexRow=1.6, cexCol=1, margins =c(6,6), Rowv=as.dendrogram(tree), Colv=as.dendrogram(cluster2), labRow=lala$species, trace=c("none"))


par(lend = 1)           # square line ends for the color legend

legend(y=1.2, x=0.8, xpd=TRUE,      # location of the legend on the heatmap plot
       legend = jColors$LABEL, # category labels
       col = jColors$COLOR,# color key
       lty= 1,             # line style
       lwd = 10,           # line width
       cex=1,
       box.lwd = 1, box.col = "white",bg = "transparent"
)

dev.off()
