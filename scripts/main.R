# library 
library(Seurat)
library(ggplot2)
library(patchwork)
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
saveRDS(sc_data, "data/sc_data.rds")
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


# Highlight Timepoint

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
p1 <- VlnPlot(sc_data, features = mam, pt.size = 0, ncol = 4)+ general_theme
p1
dev.off()
