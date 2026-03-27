# library 
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#set working dir####
setwd("D:/PDAC_liver_metastasis")


#set themes ###
general_theme = theme(text = element_text(size = 12, family = "ArialMT"))
# I used the same themes
umap_theme = theme(plot.title = element_text(hjust = 0.5), 
                   legend.position = "none", 
                   axis.text = element_blank(), 
                   axis.title = element_text(size = 12), 
                   text = element_text(size = 12, family = "ArialMT"),
                   axis.ticks = element_blank())
sample_cols <- c(
  GSM6622231_s1 = "#0072B2",  # blue
  GSM6622232_s2 = "#E69F00",  # orange
  GSM6622233_s3 = "#009E73",  # green
  GSM6622234_s4 = "#CC79A7"   # purple
)

# Random seed ----
RandomSeed <- 999

# load data ----
ss<- readxl::read_xlsx("data/ss.xlsx")
data_dir <- "data/samples"

sc_list <- create_seurat_list(ss, data_path = data_dir)

  
# Quality Control ----
# add percent.mt----
for (id in names(sc_list)){
  
  sample <- sc_list[[id]]
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt") 
  sc_list[[id]]<-sample
}
  
# Visualize QC metrics ----
pdf("analysis/QCVlnPlot.pdf", width = 12, height = 6)
for (id in names(sc_list)){
  p <- VlnPlot(sc_list[[id]], 
               features = c("nCount_RNA", "nFeature_RNA","percent.mt"))+ general_theme
  print(p)
}
dev.off()


# Filter -----
min_features <- 200
max_features <- 6000
max_mt_pct <- 10

for (id in names(sc_list)){
  sample <- sc_list[[id]]
  sample <- subset(
    sample,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_mt_pct
  )
  sc_list[[id]]<-sample
  
}

pdf("analysis/QCVlnPlot_PostFiler.pdf", width = 12, height = 6)
for (id in names(sc_list)){
  p <- VlnPlot(sc_list[[id]], 
               features = c("nCount_RNA", "nFeature_RNA","percent.mt"))+ general_theme
  print(p)
}
dev.off()


# Visualize scatter plot ----
pdf("analysis/QCScattePlot_PostFilter.pdf", width = 12, height = 6)
for (id in names(sc_list)){
  p1 <- FeatureScatter(sc_list[[id]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p2 <- FeatureScatter(sc_list[[id]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1+p2)
}
dev.off()


# Merge seurat objects ----
sc_data = merge(x = sc_list[[1]], y = sc_list[-1])
# LogNormalize and Scaling ----
sc_data = NormalizeData(sc_data, scale.factor = 10000)
sc_data = ScaleData(object = sc_data)
sc_data = FindVariableFeatures(object = sc_data, nfeatures = 2000)

# Retained only cells expressing macrophage markers Cd68 or Adgre1----
sc_data[["Adgre1"]] <- PercentageFeatureSet(sc_data, features = "Adgre1")
sc_data[["Cd68"]] <- PercentageFeatureSet(sc_data, features = "Cd68")


# Filter ----
sc_data <- sc_data <- subset(sc_data, subset = Adgre1 > 0 & Cd68 >0)

dim(sc_data)
print("31053   9167") # I have slightly higher cells, as I did not filer Ribo

RandomSeed <-999


# Dim Reduction and Clustering ----
sc_data <- RunPCA(sc_data, seed.use = RandomSeed)
ElbowPlot(sc_data)
dims_to_use <- 1:20

# cluster ----

sc_data <- FindNeighbors(object = sc_data, 
                                   dims = dims_to_use)
sc_data <- FindClusters(object = sc_data, 
                                  resolution = 0.44, 
                                  random.seed = RandomSeed)
sc_data = RunUMAP(sc_data, reduction = "pca", dims = dims_to_use, seed.use = RandomSeed)

png("analysis/umap.png", res=300, width = 6*300, height = 6*300)
DimPlot(sc_data, pt.size = 0.5, order = T, label = T,
        label.size = 6, label.color = "black") +umap_theme
dev.off()

png("analysis/umapBySample.png",res=300, width = 6*300, height = 6*300)
DimPlot(sc_data, pt.size = 0.5, order = F, label = T,
        label.size = 5, label.color = "black", group.by = "group") &umap_theme
dev.off()

#save----
saveRDS(sc_data, "data/sc_data.rds")

# Highlight Timepoint ----

png("analysis/TimepointDist.png", res=300, width = 10*300, height = 4*300)
p1 <- DimPlot(sc_data,
              cells.highlight = colnames(sc_data)[sc_data$timepoint == "Naïve"],
              cols.highlight = "#CC79A7",
              label = T , ncol = 4) &umap_theme

p2 <- DimPlot(sc_data,
              cells.highlight = colnames(sc_data)[sc_data$timepoint == "Early"],
              cols.highlight = "#009E73",
              label = T ) +umap_theme
p3 <- DimPlot(sc_data,
              cells.highlight = colnames(sc_data)[sc_data$timepoint == "Advanced"],
              cols.highlight = "#0072B2",
              label = T ) +umap_theme
p1+p2+p3
dev.off()

# annotation #####
table(sc_data$group)
print ("Tissue resident Kupffer cells (KCs) and monocyte-derived macrophages (MoMs) comprise the macrophage population in the liver. 
        KC markers such as Clec4f, Vsig4, and Timd4. MoM marker Ccr2 ")

mam <- c("Clec4f", "Vsig4",  "Timd4", "Ccr2")


p1 <- FeaturePlot(sc_data, features = mam, label = T, label.size = 5, ncol = 4)&umap_theme 
png("analysis/KcandMomMarkers.png", res=300, width = 12*300, height = 3*300)
p1
dev.off()


# The markers disp across cells ----


png("analysis/VlnKcandMomMarkers.png", res=300, width = 12*300, height = 3*300)
p1 <- VlnPlot(sc_data, features = mam, pt.size = 0, ncol = 4)
p1
dev.off()

# Percentage of each cell type in cluster -----
ss<- readxl::read_xlsx("ss.xlsx", sheet = 2)
## cells

# differential analysis per ---
sc_data<- JoinLayers(sc_data)
marker <- FindAllMarkers(sc_data)

write.csv(marker,"analysis/markers.csv")


topn <- markers %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5)
write.csv(topn,"analysis/Topmarkers.csv")

head(topn)
hm_mtx <- AverageExpression(
  sc_data,
  features = topn$gene,
  group.by = "seurat_clusters"
)$RNA

hm_mtx <- t(scale(t(hm_mtx)))
colnames(hm_mtx) <- c("0", "1", "2" ,"3" ,"4", "5", "6", "7", "8" ,"9")
annotation_df <- data.frame(ss)
rownames(annotation_df)<- annotation_df$seurat_cluster
annotation_df$seurat_cluster <- NULL
topn <- topn %>% arrange(cluster)

hm_mtx <- hm_mtx[topn$gene, ]

# Calculate gaps
gap_rows <- cumsum(rle(as.character(topn$cluster))$lengths)

graphics.off()

png("analysis/hm.png", res=300, width = 10*300, height = 10*300)
p<-pheatmap::pheatmap(hm_mtx, annotation_col = annotation_df, cluster_cols = F,
                      cluster_rows = F,
                      gaps_row = gap_rows)
p
dev.off()


# GO enrichment analysis using BP ----
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

cluster_id <- unique(markers$cluster)

GO_list <- list()

for (c in cluster_id) {
  
  genes <- markers %>%
    filter(cluster == c,
           avg_log2FC > 1,
           p_val_adj < 0.01) %>%
    pull(gene)   
  
  GO_results <- enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
  )
  
  GO_list[[as.character(c)]] <- GO_results
}

# plot 
Plots_file <- "analysis/GO_results"
dir.create(Plots_file, recursive = TRUE, showWarnings = FALSE)

for (id in names(GO_list)) {
  
  p <- barplot(GO_list[[id]], showCategory = 5)
  
  png(file.path(Plots_file, paste0(id, ".png")),
      width = 2000, height = 1600, res = 300)
  
  print(p)
  
  dev.off()
}



# what is the distrubution of these clusters -----

df <- data.frame(table( sc_data$timepoint, sc_data$seurat_clusters))
colnames(df)<- c("classification","cluster","freq")
df<-df %>% group_by(cluster) %>% mutate(sum=sum(freq))
df<-df %>% mutate(pct=round(freq/sum,2)*100)
png("analysis/clustercompositin.png", res=300, width = 6*300, height = 3*300)
ggplot (df, aes(x=cluster, y=pct, fill = classification))+ geom_col()
dev.off()


df <- data.frame(table( sc_data$seurat_clusters, sc_data$timepoint))
colnames(df)<- c("cluster","classification","freq")
df<-df %>% group_by(classification) %>% mutate(sum=sum(freq))
df<-df %>% mutate(pct=round(freq/sum,2)*100)
png("analysis/ClassificationDistribution.png", res=300, width = 6*300, height = 3*300)
ggplot (df, aes(x=classification, y=pct, fill = cluster))+ geom_col() +theme_classic()
dev.off()

# KC vs. MOM in advance vs. early disease ----

png("analysis/MarkerAcrossClass.png", res=300, width = 9*300, height = 3*300)
VlnPlot(sc_data, features =mam, group.by  = "group", pt.size = 0,ncol  = 4 )

dev.off()


# highlight pMAM vs. dMAM---
png("analysis/MAMvsdMAM-.png", res=300, width = 8*300, height = 4*300)

p2 <- DimPlot(sc_data,
              cells.highlight = colnames(sc_data)[sc_data$group == "Advanced Proximal"],
              cols.highlight = "red",
              label = T ) +umap_theme
p3 <- DimPlot(sc_data,
              cells.highlight = colnames(sc_data)[sc_data$group == "Advanced Distal"],
              cols.highlight = "#0072B2",
              label = T ) +umap_theme
p2+p3
dev.off()



# update names ----
ss<- readxl::read_xlsx("ss.xlsx",sheet = "Sheet2")
new.cluster.ids = ss$KC_or_MOM_mapp

names(new.cluster.ids) = levels(sc_data)
sc_data = RenameIdents(sc_data, new.cluster.ids)
sc_data = AddMetaData(sc_data, sc_data@active.ident, col.name = "cell_types")
sc_data = SetIdent(sc_data, value = sc_data$seurat_clusters)



##===============================================
# Part Two: Analyzing pMAM in depth 
#================================================
# check seurat identity 
Idents(sc_data)= sc_data$seurat_clusters
# Visualize Key MoM marker on the Feature plot----------
heatmap_pal2 <- c("#0D0887","#7E03A8","#CC4678","#F0F921")

# define expression min and max to AP and M2- related genes
features <- c("H2-Eb1","H2-Ab1","Cd74","Chil3","Mrc1","Arg1")
expr <- FetchData(sc_data, vars = features)
min_val <- min(expr)
max_val <- max(expr)

#plot umap next to feature plot
p <- DimPlot(sc_data, label = T, label.size = 5, pt.size = 1)+umap_theme
p_ap_m2<- FeaturePlot(
  sc_data,
  features = features,
  min.cutoff = min_val,
  max.cutoff = max_val,
  combine = TRUE
) & 
  scale_color_gradientn(colours = heatmap_pal2, limits = c(min_val, max_val))

# I am collecting all confusing labels into one unified using plot_layout
p_ap_m2 <- p_ap_m2 + plot_layout(guides = "collect") &
  theme(legend.position = "right")

final_plot <- p_ap_m2 + p + plot_layout(guides = "collect")&theme(plot.title = element_text(hjust = 0.5), 
                                                                  axis.text = element_blank(), 
                                                                  axis.title = element_text(size = 12), 
                                                                  text = element_text(size = 12, family = "ArialMT"),
                                                                  axis.ticks = element_blank())
png("analysis/AP_M2_signiture.png", res=300, width = 8*300, height = 8*300) 
final_plot

dev.off()
















# Subset the pMAM clusters 2,3,6---------------
sc_prox <- subset(sc_data, subset = seurat_clusters %in% c(2,3,6) )


# Normalize, scale , dimreduction and recluster----
# LogNormalize and Scaling ----
sc_prox = NormalizeData(sc_prox, scale.factor = 10000)
sc_prox = ScaleData(object = sc_prox)
sc_prox = FindVariableFeatures(object = sc_prox, nfeatures = 2000)
# Dim Reduction and Clustering ----
sc_prox <- RunPCA(sc_prox, seed.use = RandomSeed)
ElbowPlot(sc_prox)
dims_to_use <- 1:20

# cluster ----

sc_prox <- FindNeighbors(object = sc_prox, 
                         dims = dims_to_use)
sc_prox <- FindClusters(object = sc_prox, 
                        resolution = 0.44, 
                        random.seed = RandomSeed)
sc_prox = RunUMAP(sc_prox, reduction = "pca", dims = dims_to_use, seed.use = RandomSeed)




#rename the clusters into alphabets to make it easier to annotate----
ss<- readxl::read_xlsx("ss.xlsx", sheet = "Sheet3")

new.cluster.ids = ss$alpha_name

names(new.cluster.ids) = levels(sc_prox)
sc_prox = RenameIdents(sc_prox, new.cluster.ids)
sc_prox = AddMetaData(sc_prox, sc_prox@active.ident, col.name = "symbols")
Idents(sc_prox)




# Plot ----
p <- DimPlot(sc_prox, label = T, label.size = 5, pt.size = 1)&umap_theme
png("analysis/Mom_umap.png", res=300, width = 5*300, height = 5*300)
p
dev.off()




# Visualize Early and Advanced Tumorr====
cell_to_highlight_advance <-colnames(sc_prox)[sc_prox$timepoint == "Advanced"]
cell_to_highlight_early <-colnames(sc_prox)[sc_prox$timepoint == "Early"]


png("analysis/Early_vs_Advanced_mom.png", res=300, height = 5*300, width = 10*300)
p_advance <- DimPlot(sc_prox,
                     cells.highlight = cell_to_highlight_advance,
                     cols.highlight = "salmon",
                     label = T, pt.size = 2, sizes.highlight = 2) +umap_theme
p_ealry <- DimPlot(sc_prox,
                   cells.highlight =cell_to_highlight_early,
                   cols.highlight = "#0072B2",
                   label = T , pt.size = 2, sizes.highlight = 2) +umap_theme
p_advance+p_ealry
dev.off()




# Visualize which cell is early and which is not -----
p<- dittoSeq::dittoBarPlot(sc_prox, var = sc_prox$timepoint,
                           group.by = Idents(sc_prox),scale =  "count")
?dittoSeq::dittoBarPlot
png("analysis/cluster_composition.png", res=300, height = 5*300, width = 10*300)

p
dev.off()


# comparing new clusters with old clusters
png("analysis/comparing new clusters with old clusters.png", res=300, height = 5*300, width = 5*300)
DimPlot(sc_prox,group.by = "cell_types", pt.size = 1)&theme(plot.title = element_text(hjust = 0.5), 
                                                            axis.text = element_blank(), 
                                                            axis.title = element_text(size = 12), 
                                                            text = element_text(size = 12, family = "ArialMT"),
                                                            axis.ticks = element_blank())
dev.off()

# Markers for Mom subset------------------
unique_clusters <- levels(sc_prox$symbols)
markers<- FindAllMarkers(sc_prox)
# save it as an xlsx object ----

write.csv(markers,"analysis/Mom_markers.csv")


topn_mom <- markers %>%
  filter(avg_log2FC > 1 & p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10)

# Heatmap----



genes <- c("H2-Eb1","H2-Ab1","Cd74",
                       "Chil3","Mrc1","Arg1","Cd300a",
                       "Cd36", "Scarb1", "Anxa2" ,
                       "Gpnmb", "Lipa", "Ctsb", "Ctsd", "Psap", "Grn",
                       "Il6", "Tnf" , "Il1a",
  "Chil3","Mrc1","Arg1","Socs3","Thbs1","Trem2","Mki67",
  "H2-Ab1","H2-Eb1","Cd74","H2-Aa",
  "Cd300a","Cd36","Scarb1","Anxa2","Gpnmb","Fabp5","Abcg1","Lipa","Ctsb","Ctsd","Psap","Grn",
  "Birc5","Tubb5","H2afv","Hmgb1","Mki67","Stmn1",
  "Il6","Nfkbia","Nr4a1","Nfkbiz","Nlrp3","Tnf","Egr1","Il1a","Il10","Ccl6","Ccl9",
  "Cxcl9","Cxcl10","Ccl5","Cd40","Stat1","Socs1","H2-Q7","H2-K1","H2-T23"
)
genes <- unique(genes)
hm_mtx <- AverageExpression(
  sc_prox,
  features = genes,
  group.by = "symbols"
)$RNA

hm_mtx <- t(scale(t(hm_mtx)))
colnames(hm_mtx) <- ss$alpha_name
annotation_df <- data.frame(ss)
rownames(annotation_df)<- annotation_df$seurat_cluster
annotation_df$seurat_cluster <- NULL

png("analysis/Final_HM_mom.png",
    res=300,
    height = 8*300,
    width = 5*300)
pheatmap::pheatmap(hm_mtx, cluster_rows = T, cluster_cols = F)
dev.off()

# Cluster D is new . I think it contains KC MAM ----
p<- VlnPlot(sc_prox, features = c( c("Clec4f", "Vsig4",  "Timd4", "Ccr2")
))
graphics.off()
p
png("analysis/ClusterD_is_KC.png",
    res=300,
    height = 8*300,
    width = 5*300)
p
dev.off()
p<- VlnPlot(sc_prox, features = c( c("Clec4f", "Vsig4",  "Timd4", "Ccr2")))
graphics.off()
p
png("analysis/ClusterD_is_KC.png",
    res=300,
    height = 4*300,
    width = 8*300)
p
dev.off()




#---------------------------------------------
# Pseudotime 
#---------------------------------------------
expr <- GetAssayData(sc_prox, slot = "data")
library(TSCAN)

proc <- preprocess(expr)
order <- TSCANorder(clust)
pseudotime <- rank(order)
sc_prox$pseudotime <- pseudotime[Cells(sc_prox)]
FeaturePlot(sc_prox, features = "pseudotime")
# order from cluster G
order <- TSCANorder(clust="G", flip = TRUE)
VlnPlot(sc_prox, features = "pseudotime", group.by = "seurat_clusters")
#Which genes drive macrophage evolution?
  diff <- difftest(proc, order)
