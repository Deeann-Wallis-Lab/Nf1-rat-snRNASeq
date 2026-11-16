data <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g1 <- CreateSeuratObject(counts = data$`Gene Expression`)
g1[["HTO"]] <- CreateAssayObject(counts = data$`Antibody Capture`)

DefaultAssay(g1) <- "RNA"

g1 <- NormalizeData(g1)
g1 <- FindVariableFeatures(g1, selection.method = "mean.var.plot")
g1 <- ScaleData(g1, features = VariableFeatures(g1))

DefaultAssay(g1) <- "HTO"
g1 <- NormalizeData(g1, assay = "HTO", normalization.method = "CLR")
g1 <- HTODemux(g1, assay = "HTO", positive.quantile = 0.99)

table(g1$HTO_classification.global)

Idents(g1) <- "HTO_maxID"
RidgePlot(g1, assay = "HTO", features = rownames(g1[["HTO"]])[1:3], ncol = 3)
Idents(g1) <- "HTO_classification.global"

g1.singlet <-  subset(g1, idents = "Singlet")

Idents(g1.singlet) <- "HTO_classification"

KI_32542 <- subset(g1.singlet, idents = "KI-32542")
KI_32542$orig.ident <- "KI_32542"
KI_32542[["percent.mt"]] <- PercentageFeatureSet(KI_32542, pattern = "^Mt-")
KI_32542 <- subset(KI_32542, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

KI_32550 <- subset(g1.singlet, idents = "KI-32550")
KI_32550$orig.ident <- "KI_32550"
KI_32550[["percent.mt"]] <- PercentageFeatureSet(KI_32550, pattern = "^Mt-")
KI_32550 <- subset(KI_32550, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

KO_33097 <- subset(g1.singlet, idents = "KO-33097")
KO_33097$orig.ident <- "KO_33097"
KO_33097[["percent.mt"]] <- PercentageFeatureSet(KO_33097, pattern = "^Mt-")
KO_33097 <- subset(KO_33097, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)


# next pooled sample

data2 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g2 <- CreateSeuratObject(counts = data2$`Gene Expression`)
g2[["HTO"]] <- CreateAssayObject(counts = data2$`Antibody Capture`)

DefaultAssay(g2) <- "RNA"

g2 <- NormalizeData(g2)
g2 <- FindVariableFeatures(g2, selection.method = "mean.var.plot")
g2 <- ScaleData(g2, features = VariableFeatures(g2))

DefaultAssay(g2) <- "HTO"
g2 <- NormalizeData(g2, assay = "HTO", normalization.method = "CLR")
g2 <- HTODemux(g2, assay = "HTO", positive.quantile = 0.99)

table(g2$HTO_classification.global)

Idents(g2) <- "HTO_maxID"
RidgePlot(g2, assay = "HTO", features = rownames(g2[["HTO"]])[1:4], ncol = 4)
Idents(g2) <- "HTO_classification.global"

g2.singlet <-  subset(g2, idents = "Singlet")

Idents(g2.singlet) <- "HTO_classification"

WT_pool1 <- subset(g2.singlet, idents = "WT-pool1")
WT_pool1$orig.ident <- "WT_pool1"
WT_pool1[["percent.mt"]] <- PercentageFeatureSet(WT_pool1, pattern = "^Mt-")
WT_pool1 <- subset(WT_pool1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)

WT_pool2 <- subset(g2.singlet, idents = "WT-pool2")
WT_pool2$orig.ident <- "WT_pool2"
WT_pool2[["percent.mt"]] <- PercentageFeatureSet(WT_pool2, pattern = "^Mt-")
WT_pool2 <- subset(WT_pool2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)


KO_41092 <- subset(g2.singlet, idents = "KO-41092")
KO_41092$orig.ident <- "KO_41092"
KO_41092[["percent.mt"]] <- PercentageFeatureSet(KO_41092, pattern = "^Mt-")
KO_41092 <- subset(KO_41092, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

### not used since it has a uneven distribution of cells (found in old analysis)
KO_37109 <- subset(g2.singlet, idents = "KO-37109")
KO_37109$orig.ident <- "KO_37109"
KO_37109[["percent.mt"]] <- PercentageFeatureSet(KO_37109, pattern = "^Mt-")
KO_37109 <- subset(KO_37109, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

# next pooled sample

data3 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g3 <- CreateSeuratObject(counts = data3$`Gene Expression`)
g3[["HTO"]] <- CreateAssayObject(counts = data3$`Antibody Capture`)

DefaultAssay(g3) <- "RNA"

g3 <- NormalizeData(g3)
g3 <- FindVariableFeatures(g3, selection.method = "mean.var.plot")
g3 <- ScaleData(g3, features = VariableFeatures(g3))

DefaultAssay(g3) <- "HTO"
g3 <- NormalizeData(g3, assay = "HTO", normalization.method = "CLR")
g3 <- HTODemux(g3, assay = "HTO", positive.quantile = 0.99)

table(g3$HTO_classification.global)

Idents(g3) <- "HTO_maxID"
RidgePlot(g3, assay = "HTO", features = rownames(g3[["HTO"]])[1:4], ncol = 4)
Idents(g3) <- "HTO_classification.global"

g3.singlet <-  subset(g3, idents = "Singlet")

Idents(g3.singlet) <- "HTO_classification"

WT_41042 <- subset(g3.singlet, idents = "WT-41042")
WT_41042$orig.ident <- "WT_41042"
WT_41042[["percent.mt"]] <- PercentageFeatureSet(WT_41042, pattern = "^Mt-")
WT_41042 <- subset(WT_41042, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)

# not enough cells did not use
WT_41014 <- subset(g3.singlet, idents = "WT-41014")
WT_41014$orig.ident <- "WT_41014"
WT_41014[["percent.mt"]] <- PercentageFeatureSet(WT_41014, pattern = "^Mt-")
WT_41014 <- subset(WT_41014, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

KO_32563 <- subset(g3.singlet, idents = "KO-32563")
KO_32563$orig.ident <- "KO_32563"
KO_32563[["percent.mt"]] <- PercentageFeatureSet(KO_32563, pattern = "^Mt-")
KO_32563 <- subset(KO_32563, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

KI_32521 <- subset(g3.singlet, idents = "KI-32521")
KI_32521$orig.ident <- "KI_32521"
KI_32521[["percent.mt"]] <- PercentageFeatureSet(KI_32521, pattern = "^Mt-")
KI_32521 <- subset(KI_32521, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)



# next pooled sample

data4 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g4 <- CreateSeuratObject(counts = data4$`Gene Expression`)
g4[["HTO"]] <- CreateAssayObject(counts = data4$`Antibody Capture`)

DefaultAssay(g4) <- "RNA"

g4 <- NormalizeData(g4)
g4 <- FindVariableFeatures(g4, selection.method = "mean.var.plot")
g4 <- ScaleData(g4, features = VariableFeatures(g4))

DefaultAssay(g4) <- "HTO"
g4 <- NormalizeData(g4, assay = "HTO", normalization.method = "CLR")
g4 <- HTODemux(g4, assay = "HTO", positive.quantile = 0.99)

table(g4$HTO_classification.global)

Idents(g4) <- "HTO_maxID"
RidgePlot(g4, assay = "HTO", features = rownames(g4[["HTO"]])[1:4], ncol = 4)
Idents(g4) <- "HTO_classification.global"

g4.singlet <-  subset(g4, idents = "Singlet")

Idents(g4.singlet) <- "HTO_classification"

WT_pool3 <- subset(g4.singlet, idents = "WT-pool3")
WT_pool3$orig.ident <- "WT_pool3"
WT_pool3[["percent.mt"]] <- PercentageFeatureSet(WT_pool3, pattern = "^Mt-")
WT_pool3 <- subset(WT_pool3, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)

KI_42280 <- subset(g4.singlet, idents = "KI-42280")
KI_42280$orig.ident <- "KI_42280"
KI_42280[["percent.mt"]] <- PercentageFeatureSet(KI_42280, pattern = "^Mt-")
KI_42280 <- subset(KI_42280, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)

KO_42270 <- subset(g4.singlet, idents = "KO-42270")
KO_42270$orig.ident <- "KO_42270"
KO_42270[["percent.mt"]] <- PercentageFeatureSet(KO_42270, pattern = "^Mt-")
KO_42270 <- subset(KO_42270, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

# next pooled sample

data5 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g5 <- CreateSeuratObject(counts = data5$`Gene Expression`)
g5[["HTO"]] <- CreateAssayObject(counts = data5$`Antibody Capture`)

DefaultAssay(g5) <- "RNA"

g5 <- NormalizeData(g5)
g5 <- FindVariableFeatures(g5, selection.method = "mean.var.plot")
g5 <- ScaleData(g5, features = VariableFeatures(g5))

DefaultAssay(g5) <- "HTO"
g5 <- NormalizeData(g5, assay = "HTO", normalization.method = "CLR")
g5 <- HTODemux(g5, assay = "HTO", positive.quantile = 0.99)

table(g5$HTO_classification.global)

Idents(g5) <- "HTO_maxID"
RidgePlot(g5, assay = "HTO", features = rownames(g5[["HTO"]])[1:4], ncol = 4)
Idents(g5) <- "HTO_classification.global"

g5.singlet <-  subset(g5, idents = "Singlet")

Idents(g5.singlet) <- "HTO_classification"

KI_42268 <- subset(g5.singlet, idents = "KI-42268")
KI_42268$orig.ident <- "KI_42268"
KI_42268[["percent.mt"]] <- PercentageFeatureSet(KI_42268, pattern = "^Mt-")
KI_42268 <- subset(KI_42268, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

KI_42278 <- subset(g5.singlet, idents = "KI-42278")
KI_42278$orig.ident <- "KI_42278"
KI_42278[["percent.mt"]] <- PercentageFeatureSet(KI_42278, pattern = "^Mt-")
KI_42278 <- subset(KI_42278, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

KO_43880 <- subset(g5.singlet, idents = "KO-43880")
KO_43880$orig.ident <- "KO_43880"
KO_43880[["percent.mt"]] <- PercentageFeatureSet(KO_43880, pattern = "^Mt-")
KO_43880 <- subset(KO_43880, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)


# next pooled sample

data6 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g6 <- CreateSeuratObject(counts = data6$`Gene Expression`)
g6[["HTO"]] <- CreateAssayObject(counts = data6$`Antibody Capture`)

DefaultAssay(g6) <- "RNA"

g6 <- NormalizeData(g6)
g6 <- FindVariableFeatures(g6, selection.method = "mean.var.plot")
g6 <- ScaleData(g6, features = VariableFeatures(g6))

DefaultAssay(g6) <- "HTO"
g6 <- NormalizeData(g6, assay = "HTO", normalization.method = "CLR")
g6 <- HTODemux(g6, assay = "HTO", positive.quantile = 0.99)

table(g6$HTO_classification.global)

Idents(g6) <- "HTO_maxID"
RidgePlot(g6, assay = "HTO", features = rownames(g6[["HTO"]])[1:4], ncol = 4)
Idents(g6) <- "HTO_classification.global"

g6.singlet <-  subset(g6, idents = "Singlet")

Idents(g6.singlet) <- "HTO_classification"

WT_pool4 <- subset(g6.singlet, idents = "WT-pool4")
WT_pool4$orig.ident <- "WT_pool4"
WT_pool4[["percent.mt"]] <- PercentageFeatureSet(WT_pool4, pattern = "^Mt-")
WT_pool4 <- subset(WT_pool4, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)

KI_41716 <- subset(g6.singlet, idents = "KI-41716")
KI_41716$orig.ident <- "KI_41716"
KI_41716[["percent.mt"]] <- PercentageFeatureSet(KI_41716, pattern = "^Mt-")
KI_41716 <- subset(KI_41716, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

# did not use in final dataset too few cells
KO_44538 <- subset(g6.singlet, idents = "KO-44538")
KO_44538$orig.ident <- "KO_44538"
KO_44538[["percent.mt"]] <- PercentageFeatureSet(KO_44538, pattern = "^Mt-")
KO_44538 <- subset(KO_44538, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)

# next pooled sample

data7 <- Read10X_h5(filename =  "filtered_feature_bc_matrix.h5")

g7 <- CreateSeuratObject(counts = data7$`Gene Expression`)
g7[["HTO"]] <- CreateAssayObject(counts = data7$`Antibody Capture`)

DefaultAssay(g7) <- "RNA"

g7 <- NormalizeData(g7)
g7 <- FindVariableFeatures(g7, selection.method = "mean.var.plot")
g7 <- ScaleData(g7, features = VariableFeatures(g7))

DefaultAssay(g7) <- "HTO"
g7 <- NormalizeData(g7, assay = "HTO", normalization.method = "CLR")
g7 <- HTODemux(g7, assay = "HTO", positive.quantile = 0.99)

table(g7$HTO_classification.global)

Idents(g7) <- "HTO_maxID"
RidgePlot(g7, assay = "HTO", features = rownames(g7[["HTO"]])[1:4], ncol = 4)
Idents(g7) <- "HTO_classification.global"

g7.singlet <-  subset(g7, idents = "Singlet")

Idents(g7.singlet) <- "HTO_classification"

# did not use in final dataset too few cells
WT_pool5 <- subset(g7.singlet, idents = "WT-pool5")
WT_pool5$orig.ident <- "WT_pool5"
WT_pool5[["percent.mt"]] <- PercentageFeatureSet(WT_pool5, pattern = "^Mt-")
WT_pool5 <- subset(WT_pool5, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)
# did not use in final dataset too few cells
KO_1098 <- subset(g7.singlet, idents = "KO-1098")
KO_1098$orig.ident <- "KO_1098"
KO_1098[["percent.mt"]] <- PercentageFeatureSet(KO_1098, pattern = "^Mt-")
KO_1098 <- subset(KO_1098, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
# did not use in final dataset too few cells
KO_1024 <- subset(g7.singlet, idents = "KO-1024")
KO_1024$orig.ident <- "KO_1024"
KO_1024[["percent.mt"]] <- PercentageFeatureSet(KO_1024, pattern = "^Mt-")
KO_1024 <- subset(KO_1024, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)


# put groups together

KI <- merge(x = KI_32521, y = c(KI_32542,KI_32550,KI_41716,KI_42268,KI_42278,KI_42280), add.cell.ids = c("KI_32521","KI_32542","KI_32550", "KI_41716" ,"KI_42268" ,"KI_42278" ,"KI_42280"), project = "KI")

KI$orig.ident <- "KI"

KO <- merge(x = KO_33097, y = c(KO_41092,KO_32563,KO_42270, KO_43880 ), add.cell.ids = c("KO_33097", "KO_41092", "KO_32563", "KO_42270" , "KO_43880" ), project = "KO")
KO$orig.ident <- "KO"

WT <- merge(x = WT_41042, y = c(WT_pool1,WT_pool2,WT_pool3, WT_pool4), add.cell.ids = c("WT_41042", "WT_pool1", "WT_pool2", "WT_pool3", "WT_pool4"), project = "WT")
WT$orig.ident <- "WT"

all_samples <- merge(x = KO, y = c(KI,WT), add.cell.ids = c("KO", "KI", "WT"), project = "all")


DefaultAssay(all_samples) <- "RNA"

# Quality control

all_samples[["percent.mt"]] <- PercentageFeatureSet(all_samples, pattern = "^Mt-")

# Visualize QC metrics as a violin plot
g <- VlnPlot(all_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 5, pt.size = 0)
g1 <- VlnPlot(all_samples, features = "percent.mt", ncol = 5, pt.size = 0, y.max = 25)
g2 <- VlnPlot(all_samples, features = "nFeature_RNA", ncol = 5, pt.size = 0, y.max = 5000)

ggsave(plot = g,
       filename = "nFeature_count_PercentMito.png",
       path = "Figures/",
       width = 20,
       height = 5,
       dpi = 300)

ggsave(plot = g1,
       filename = "PercentMito.png",
       path = "Figures/",
       width = 20,
       height = 5,
       dpi = 300)

ggsave(plot = g2,
       filename = "nFeature_RNA.png",
       path = "Figures/",
       width = 20,
       height = 5,
       dpi = 300)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# continue onto intagrating the data

Test.list <- SplitObject(all_samples, split.by = "orig.ident")

all_samples.list <- Test.list[c("KO", "KI", "WT")]

# Normalize the data for Intergration step

for (i in 1:length(all_samples.list)) {
  all_samples.list[[i]] <- NormalizeData(all_samples.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  all_samples.list[[i]] <- FindVariableFeatures(all_samples.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
}
# Create a list of Anchor points for intergration

reference.list <- all_samples.list[c("KO", "KI", "WT")]

all_samples.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)

# intergrate the data
all_samples.integrated <- IntegrateData(anchorset = all_samples.anchors, dims = 1:50)

DefaultAssay(all_samples.integrated) <- "integrated"

all_samples.integrated <- ScaleData(all_samples.integrated, verbose = T)

all_samples.integrated <- RunPCA(all_samples.integrated, npcs = 50, verbose = T)

# Visualize Variable features
top15 <- head(VariableFeatures(all_samples.integrated), 15)

top15

# Determine the dimentionality

ElbowPlot(all_samples.integrated, ndims = 50)

data <- all_samples.integrated

data <- RunHarmony(data, "HTO_classification")

DimPlot(data, reduction = "harmony", pt.size = 1, group.by = "HTO_classification")

DimPlot(data, reduction = "harmony", pt.size = 1, group.by = "HTO_classification", split.by = "HTO_classification")

data <- data %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 0.4) %>% identity()

g5 <- DimPlot(data, reduction = "umap", pt.size = 1, split.by = "HTO_classification")

ggsave(plot = g5,
       filename = "umap dim 40 res0.4 sample split.png",
       path = "Figures/",
       width = 40,
       height = 5,
       dpi = 300)

g6 <- DimPlot(data, reduction = "umap", pt.size = 0.5, group.by = "HTO_classification")

ggsave(plot = g6,
       filename = "umap dim 40 res0.4 sample.png",
       path = "Figures/",
       width = 8,
       height = 5,
       dpi = 300)

g7 <- DimPlot(data, reduction = "umap", pt.size = 0.5)


ggsave(plot = g7,
       filename = "umap dim 40 res0.4.png",
       path = "Figures/",
       width = 8,
       height = 5,
       dpi = 300)


saveRDS(data, file = "RDS/harmony_integrated.RDS")

data <- harmony_all_samples_integrated_Mut_WT

DimPlot(data)

data <- ScaleData(data, assay = "RNA")

saveRDS(data, file = "RDS/harmony_integrated_Mut_WT_ScaleRNA.RDS")