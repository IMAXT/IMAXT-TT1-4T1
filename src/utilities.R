#__Utilities script for functions and packages__#
suppressMessages(library(dsb))

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scuttle)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyverse)
  library(dsb)
  library(tibble)
  library(cowplot)
  library(pheatmap)
  library(stringr)
  library(readxl)
  library(fgsea)
  library(tibble)
  library(ggplot2)
  library(scran)
  library(scater)
  library(glue)
  library(org.Hs.eg.db)
  library(dsb)
  library(readr)
  library(tidyverse)
  library(UCell)
  library(progeny)
  library(ggrepel)
  library(reshape2)
  library(CellChat)
  library(patchwork)
  library(gridExtra)
})

##__Enrichment Function__##
#refrence for gene sets
#http://bioinf.wehi.edu.au/software/MSigDB/
# this function takes a list of data frames of markers output from DE
run_gene_enrichment <- function(marker_list, pathways_input){
  
  pathways<-readRDS(pathways_input)
  # create a list of data frames for results of enrichmet per cluster/group
  fgseaRes_list <- lapply(marker_list, function(List){ 
    
    # this extracts gene symbol and logFC from every cluster
    List<-List[, grepl("gen", colnames(List)) | grepl("avg_log2FC", colnames(List))]
    List$Gene<-toupper(List$gene)# cabiltalize mouse genes to map to human's
    mm<-select(org.Hs.eg.db, keys = List$Gene, columns = "ENTREZID", keytype = "SYMBOL")
    mm<-mm[!duplicated(mm$SYMBOL),]
    List$ID<-mm$ENTREZID
    List<-List[,c(4,1)]# rearrange so the gene IDs so they become the first column and the fold changes the 2nd, must for deframe function
    
    List_Rank<-deframe(List)
    fgseaRes <- fgsea(pathways,List_Rank,minSize=15, maxSize=500)%>%arrange(-NES)
    
    return(fgseaRes)
    
  }) 
  
  sig_pathways <- lapply(fgseaRes_list, function(x){
    
    x <- x[x$padj < 0.05 , c('pathway', 'pval', 'padj', 'NES', 'leadingEdge','size')] 
    
    return(x)
    
  })
  
  sig_pathways_name <- lapply(sig_pathways, function(x){x$pathway})
  
  output <- list(fgseaRes_list = fgseaRes_list, sig_pathways = sig_pathways, sig_pathways_name = sig_pathways_name)
  
  return(output)
  
}


##__Plotting Function__##


#gene_nerichment_list is the ouptut of the above function
make_ge_plots <- function(gene_enrichment_list, pathway_name, pthresh = 0.05, ntop = 10, y_text_size = 10){
  
  clusters <- names(gene_enrichment_list)
  
  plots <- lapply(clusters, function(cluster, fgsesRes_list, pthresh, ntop, pathway_name, y_text_size){
    
    # select top 10 significant pathways  
    fgsesRes <- fgsesRes_list[[cluster]]
    fgsesRes <- fgsesRes[fgsesRes$padj < pthresh , ]
    fgsesRes <- top_n(fgsesRes, n = ntop, wt = abs(NES))
    
    #fgsesRes$pathway <- as.character(lapply(fgsesRes$pathway, function(pw) str_replace(pw, pattern = glue(pathway_name, "_"), "")))
    
    if(dim(fgsesRes)[1] > 0){
      
      plot <- ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj), width = 0.7) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title= "cluster", cluster, "pathways")+theme_minimal()
      return(plot)
      
    } else{
      
      return(NULL)
      
    }
    
  }, fgsesRes_list = gene_enrichment_list, pthresh = pthresh, ntop = ntop, pathway_name = pathway_name, y_text_size = y_text_size)
  
  plots <- plots[unlist(lapply(plots, function(plot) ! is.null(plot)))]# looks at plots list and keeps those ewntires that are not null
  # unlist
  return(plots)
  
}

#plot_grid(plotlist = plots, ncol = 1, align = 'hv')


## converting mouse to human and vice versa

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

# converting human to mouse
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  no_mouse_genes <- length(x)
  no_human_genes <- length(humanx)
  
  if(no_human_genes != no_mouse_genes){
    print("Some genes could not be translated!")
    genes_not_trans <- setdiff(x,genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  # Print all gene names that could not be translated and the number of genes that were not translated
  
  return(humanx)
}


## quick functions! ##

#this is a visualization function for seurat object to check main QC metrics i care about

preQC<-function(seur){
  seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^mt")
  # Visualize QC metrics as a violin plot
  p1=VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  return(p1)
}


# this function to do qC on seurat based on cvisualization, allows me to dynamically change thresholds!

QC<- function(JZ,nFeature_min,nFeature_max,mt){
  JZ[["percent.mt"]] <- PercentageFeatureSet(JZ, pattern = "^mt")
  JZ=subset(JZ, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt<mt)
  return(JZ)
}

# clustering and dim reuctions
NormaVar<-function(JZ, nfeatures){
  JZ<-NormalizeData(JZ)
  JZ<- FindVariableFeatures(JZ, selection.method = "vst", nfeatures =nfeatures)
  return(JZ)
}


PCA<-function(a){
  a <- ScaleData(a, verbose = T) 
  a <- RunPCA(a,verbose = T)
  return(a)
}

Cluster<-function(a, ndims, resol){
  a <- FindNeighbors(a, dims = ndims)
  a <- FindClusters(a, resolution = resol)
  a<-RunUMAP(a,dims = ndims)
  return(a)
}

plotCombo<-function(cell_meta,x,yy){
  cluster_meta <- cell_meta%>%dplyr::select({{x}},{{yy}})%>% group_by({{x}},{{yy}}) %>% tally(name="NumberOfCells")
  p=ggplot(cluster_meta, aes(x={{x}}, y = NumberOfCells, fill = {{yy}})) +
    geom_bar(stat = 'identity')+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
  return(p)
}


quickDE=function(integrated){
  markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use ="wilcox", verbose = T)
  top10<-markers%>%group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  print(DoHeatmap(integrated, features = top10$gene, size=4, angle=90) + NoLegend())
  return(markers)
}

# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}


# this is how i saved the inpute for cellphoneDB
#write.csv(genesV2,"geneIDs.csv")



