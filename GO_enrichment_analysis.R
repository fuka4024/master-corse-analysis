if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("clusterProfiler")
BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))
BiocManager::install("GO.db")

GO_gene <- read.csv("C:/Users/monom/Documents/Ortholog_GO_Drosophila.csv")
gene_transcript <- read.csv("C:/Users/monom/Documents/Symbol_BmoriID.csv")
gene_transcript <- gene_transcript[-1,]
gene <- unique(gene_transcript$Gene)
GO <- c()
TR <- c()
for(RI in 1:length(gene)){
  Hit <- which(GO_gene$Gb_p50 == gene[RI])
  Hit_2 <- which(gene_transcript$Gene == gene[RI])
  if(length((Hit))==0){
    next
  }else{
    GO <- c(GO,rep(GO_gene$GO[Hit],length(Hit_2)))
    TR <- c(TR,rep(gene_transcript$Tb[Hit_2],each=length(Hit)))
  }
  print(RI)
}

go2gene.c2 <- as.data.frame(cbind(GO,TR))
names(go2gene.c2) <- c('term', 'gene')
head(go2gene.c2)

library(clusterProfiler)
go2gene.expand <- buildGOmap(go2gene.c2)

library(GO.db)
columns(GO.db)

gonames <- AnnotationDbi::select(GO.db, keys=keys(GO.db), columns=c("TERM", "ONTOLOGY"))

etab.transcript <- read.csv("C:/Users/monom/Documents/DE_7/FB_5_0_vs_5_3_tab.csv")
head(etab.transcript)

upreg <- subset(etab.transcript, PValue < 0.05 & FDR < 0.05 & logFC < 0)
upregGenes <- upreg$X
expScore <- sign(etab.transcript$logFC) * -log(etab.transcript$FDR)
names(expScore) <- etab.transcript$X
expScore <- sort(expScore, decreasing=TRUE)

gonames.bp <- subset(gonames, ONTOLOGY=="BP")
go2gene.expand.bp <- subset(go2gene.expand, GO %in% gonames.bp$GOID)

cprof.up.bp.fisher2 <- enricher(gene=upregGenes, TERM2GENE=go2gene.expand.bp, TERM2NAME=gonames)
significant_ora <- subset(cprof.up.bp.fisher2@result, p.adjust < 0.05)

cprof.bp.gsea2 <- GSEA(geneList=expScore, TERM2GENE=go2gene.expand.bp, TERM2NAME=gonames)
significant_gsea <- subset(cprof.bp.gsea2@result, p.adjust < 0.05)


pdf("C:/Users/monom/Documents/GO_enrichment_analysis/ORA_5_0_vs_5_3_down_MF.pdf", 
    width = 15, height = 50)

barplot(cprof.up.bp.fisher2, showCategory=min(84, nrow(significant_ora)), 
                title="Significant GO terms (p.adjust < 0.05)")

dev.off()


pdf("C:/Users/monom/Documents/GO_enrichment_analysis/GSEA_5_0_vs_5_3_down_BP.pdf", 
    width = 15, height = 50)
gseaplot(cprof.bp.gsea2, geneSetID = 1, title="GSEA for top GO term")

dev.off()
