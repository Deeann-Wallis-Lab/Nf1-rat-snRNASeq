all_samples.integrated <- harmony_all_samples_integrated_cell_type
# Macrophage
Macrophage_features <- c("Mrc1", "Cd163", "Itgam", "Rbpj", "Ms4a4a", "Ms4a6a", "Lyz", "C1qa", "Msr1","Fcer1g", "Cd4")

g <- FeaturePlot(all_samples.integrated, features = Macrophage_features)

ggsave(plot = g,
       filename = "Macrophage_features.png",
       path = ".../Figures/",
       width = 20,
       height = 15,
       dpi = 300)

# epithelial 
epithelial_features <- c("Krt8", "Krt18", "Epcam", "Cd24", "Krt7", "Krt17")

g <- FeaturePlot(all_samples.integrated, features = epithelial_features )

ggsave(plot = g,
       filename = "epithelial_features.png",
       path = ".../Figures/",
       width = 15,
       height = 15,
       dpi = 300)

# fibrocyte

fibrocyte_features <- c("Ptprc",	"Cd163",	"S100a8",	"S100a9",	"Cd90",	"Col1a1",	"Il6",	"Ccl3",	"Ccl4")

g <- FeaturePlot(all_samples.integrated, features = fibrocyte_features)

ggsave(plot = g,
       filename = "fibrocyte_features.png",
       path = ".../Figures/",
       width = 15,
       height = 15,
       dpi = 300)

# Bcell nuc
bcell_features <- c("Jchain",	"Ighg1",	"Ighg2", "Ighg3",	"Mzb1",	"Ms4a1")

g <- FeaturePlot(all_samples.integrated, features = bcell_features)

ggsave(plot = g,
       filename = "bcell_features.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

# mesenchyme_fibroblasts

mesenchyme_fibroblasts <- c("Fsp1",	"Acta2","Fgf10",	"Tbx4",	"Vim",	"Pdgfra",	"Des")
g <- FeaturePlot(all_samples.integrated, features = mesenchyme_fibroblasts)

ggsave(plot = g,
       filename = "mesenchyme_fibroblasts.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

# Mast cells

Mast_cells <- c("Tpsab1",	"Tpsb2",	"Cpa3",	"Kit",	"Fcer1g")
g <- FeaturePlot(all_samples.integrated, features = Mast_cells)

ggsave(plot = g,
       filename = "Mast_cells2.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

# endothelial cells

endothelial <- c("Pecam1"	,"Cdh5",	"Tie2",	"Foxf1")

g <- FeaturePlot(all_samples.integrated, features = endothelial)

ggsave(plot = g,
       filename = "endothelial.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

# basal cells

basal <- c("Krt5",	"Ngfr",	"P63",	"Krt14",	"Krt15"	,"Itgb4")
g <- FeaturePlot(all_samples.integrated, features = basal)

ggsave(plot = g,
       filename = "basal.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

# cancer stem cell

cancer_stem_cell <- c("Prom1",	"Kit",	"Abcg2",	"Nes", "B4galnt1",	"Fzd6")

g <- FeaturePlot(all_samples.integrated, features = cancer_stem_cell)

ggsave(plot = g,
       filename = "cancer_stem_cell.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

# secretory_cells

secretory_cell<- c("Sftpb",	"Scgb1a2",	"Scgb1a1")

g <-  FeaturePlot(all_samples.integrated, features = secretory_cell)

ggsave(plot = g,
       filename = "secretory_cell.png",
       path = ".../Figures/",
       width = 8,
       height = 5,
       dpi = 300)

# melanocyte
melanocyte <- c("Csh3", "Mlana"	,"Gja3",	"Mmp17",	"Oca2",	"Clec11a"	,"Sox10",	"Gpr137b",	"Tyr",	"Tyrp1",	"Kcnab2",	"Pir",	"Tnfrsf14",	"Farp2",	"Qpct",	"Slc45a2",	"Wipi1",	"Mcoln1",	"Wfdc1",	"Aph1b",	"Uap1l1",	"Tbc1d16")

g <- FeaturePlot(all_samples.integrated, features = melanocyte)

ggsave(plot = g,
       filename = "melanocyte.png",
       path = ".../Figures/",
       width = 15,
       height = 15,
       dpi = 300)


myoepithelial_cell <- c("Mylk", "Myl9", "Krt5", "Actg2", 	"Krt14")

g <- FeaturePlot(all_samples.integrated, features = myoepithelial_cell)

ggsave(plot = g,
       filename = "myoepithelial_cell.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)


Mast_Cell <- c('Kit', "Ltc4s", "Tpsab1", "Il1rl1", "Tpsb2", " Cma1", "Fcer2")


g <- FeaturePlot(all_samples.integrated, features = Mast_Cell)

ggsave(plot = g,
       filename = "Mast_Cell2.png",
       path = ".../Figures/",
       width = 8,
       height = 5,
       dpi = 300)

monocyte <- c("Cd40","Csf3r", "S100a4", "Cxcl10", "S100a8", "Fcgr3b", "Cxcr4", "Fcgr1a", "Mrc1", "Ly6c2", "Ly6g", "Cd7")

g <- FeaturePlot(all_samples.integrated, features = monocyte)

ggsave(plot = g,
       filename = "monocyte.png",
       path = ".../Figures/",
       width = 10,
       height = 10,
       dpi = 300)

Fibroblast <- c("Col6a2", "Lum", " Aspn", "Lox", "Postn", "Cela1", "Ccr2", "Fbn1", "Fbln2", "Igf1")

g <- FeaturePlot(all_samples.integrated, features = Fibroblast)

ggsave(plot = g,
       filename = "Fibroblast.png",
       path = ".../Figures/",
       width = 10,
       height = 10,
       dpi = 300)

Adipocyte <- c("Adipoq", "Plin1", "Pnpla2", "Lrp1", "Fabp4")

g <- FeaturePlot(all_samples.integrated, features = Adipocyte)

ggsave(plot = g,
       filename = "Adipocyte.png",
       path = ".../Figures/",
       width = 8,
       height = 8,
       dpi = 300)

Tcell_features <- c("Trac",	"Themis",	"Il7r",	"Cd247", "Cd96", "Cd8a","Cd8b", "Cd4", "Ilr2a", "Itgae",  "Tbx21","Cd3d", " Cd3E", "Cd3g", "Gzmk", "Ltb", "Trbc1", "Cd6" , "Prf1", "Gzmp")

g <- FeaturePlot(all_samples.integrated, features = Tcell_features)

ggsave(plot = g,
       filename = "Tcell_features.png",
       path = ".../Figures/",
       width = 15,
       height = 15,
       dpi = 300)


#Unidentifeid clusters PangioDB was used for 

data <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(data)

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T)

UK1_markers <- FindMarkers(data, only.pos = F, min.pct = 0.05, logfc.threshold = 0.5, assay = "RNA", ident.1 = "2")

UK1_markers$genes <- row.names(UK1_markers)

UK1_marker_list <- PanglaoDB$`official gene symbol`

#export and capitalize data

writexl::write_xlsx(UK1_markers, path = ".../Unknown 1 markers.xlsx")

UK1_marker_use <- Unknown_1_markers %>%  filter(gene_Upper %in% UK1_marker_list)

filter_terms <- UK1_marker_use$gene_Upper

Pangio_UK <- PanglaoDB %>%  filter(PanglaoDB$`official gene symbol` %in% filter_terms)

Pangio_UK <- as.data.frame(Pangio_UK)
UK1_marker_use <- as.data.frame(UK1_marker_use)


test <- merge(x = Pangio_UK, y = UK1_marker_use, by.x ="official gene symbol", by.y = "gene_Upper")




writexl::write_xlsx(test, path = ".../Unknown 1 cell type possibilities.xlsx")
### used this

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T, split.by = "orig.ident")

DefaultAssay(data) <- "RNA"

Myocytes <- subset(Unknown_1_cell_type_possibilities, Unknown_1_cell_type_possibilities$`cell type` == "Myocytes")

Myocytes <- subset(Myocytes, Myocytes$avg_log2FC > 1.5)

Myocytes <- Myocytes$genes

g1 <- FeaturePlot(data, features = Myocytes)

ggsave(plot = g1,
       filename = "Myocytes.png",
       path = ".../Figures/Cell Type/",
       width = 15,
       height = 15,
       dpi = 300, limitsize = F)

# UK2 

#data <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(data)

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T)

UK2_markers <- FindMarkers(data, only.pos = F, min.pct = 0.05, logfc.threshold = 1.0, assay = "RNA", ident.1 = "4")

UK2_markers$genes <- row.names(UK2_markers)

UK2_marker_list <- PanglaoDB$`official gene symbol`

#export and capitalize data

writexl::write_xlsx(UK2_markers, path = ".../Unknown 2 markers.xlsx")

UK2_marker_use <- Unknown_2_markers %>%  filter(gene_upper %in% UK2_marker_list)

filter_terms <- UK2_marker_use$gene_upper

Pangio_UK <- PanglaoDB %>%  filter(PanglaoDB$`official gene symbol` %in% filter_terms)

Pangio_UK <- as.data.frame(Pangio_UK)
UK2_marker_use <- as.data.frame(UK2_marker_use)


test <- merge(x = Pangio_UK, y = UK2_marker_use, by.x ="official gene symbol", by.y = "gene_upper")




writexl::write_xlsx(test, path = ".../Unknown 2 cell type possibilities.xlsx")
### used this

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T, split.by = "orig.ident")

DefaultAssay(data) <- "RNA"

Myoepithelial_cells <- subset(Unknown_2_cell_type_possibilities, Unknown_2_cell_type_possibilities$`cell type` == "Myoepithelial cells")

Myoepithelial_cells <- subset(Myoepithelial_cells, Myoepithelial_cells$avg_log2FC >= 1)

Myoepithelial_cells <- Myoepithelial_cells$genes

g1 <- FeaturePlot(data, features = Myoepithelial_cells)

ggsave(plot = g1,
       filename = "Myoepithelial_cells.png",
       path = ".../Figures/Cell Type/",
       width = 15,
       height = 15,
       dpi = 300, limitsize = F)
#UK4

#data <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(data)

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T)

UK4_markers <- FindMarkers(data, only.pos = F, min.pct = 0.40, logfc.threshold = 1.0, assay = "RNA", ident.1 = "9")

UK4_markers$genes <- row.names(UK4_markers)

UK4_marker_list <- PanglaoDB$`official gene symbol`


#export and capitalize data

writexl::write_xlsx(UK4_markers, path = ".../Unknown 4 markers.xlsx")

UK4_marker_use <- Unknown_4_markers %>%  filter(gene_upper %in% UK4_marker_list)

filter_terms <- UK4_marker_use$gene_upper

Pangio_UK <- PanglaoDB %>%  filter(PanglaoDB$`official gene symbol` %in% filter_terms)

Pangio_UK <- as.data.frame(Pangio_UK)
UK4_marker_use <- as.data.frame(UK4_marker_use)


test <- merge(x = Pangio_UK, y = UK4_marker_use, by.x ="official gene symbol", by.y = "gene_upper")




writexl::write_xlsx(test, path = ".../Unknown 4 cell type possibilities.xlsx")
### used this

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T, split.by = "orig.ident")

DefaultAssay(data) <- "RNA"

Endothelial_cells <- subset(Unknown_4_cell_type_possibilities, Unknown_4_cell_type_possibilities$`cell type` == "Endothelial cells")

Endothelial_cells <- subset(Endothelial_cells, Endothelial_cells$avg_log2FC >= 1)

Endothelial_cells <- Endothelial_cells$genes

g1 <- FeaturePlot(data, features = Endothelial_cells)

ggsave(plot = g1,
       filename = "Endothelial_cells.png",
       path = ".../Figures/Cell Type/",
       width = 15,
       height = 15,
       dpi = 300, limitsize = F)
# UK5

#data <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(data)

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T)

UK5_markers <- FindMarkers(data, only.pos = F, min.pct = 0.40, logfc.threshold = 1.0, assay = "RNA", ident.1 = "17")

UK5_markers$genes <- row.names(UK5_markers)

UK5_marker_list <- PanglaoDB$`official gene symbol`


#export and capitalize data

writexl::write_xlsx(UK5_markers, path = ".../Unknown 5 markers.xlsx")

UK5_marker_use <- Unknown_5_markers %>%  filter(gene_upper %in% UK5_marker_list)

filter_terms <- UK5_marker_use$gene_upper

Pangio_UK <- PanglaoDB %>%  filter(PanglaoDB$`official gene symbol` %in% filter_terms)

Pangio_UK <- as.data.frame(Pangio_UK)
UK5_marker_use <- as.data.frame(UK5_marker_use)


test <- merge(x = Pangio_UK, y = UK5_marker_use, by.x ="official gene symbol", by.y = "gene_upper")




writexl::write_xlsx(test, path = ".../Unknown 5 cell type possibilities.xlsx")
### used this

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T, split.by = "orig.ident")

DefaultAssay(data) <- "RNA"

B_cells <- subset(Unknown_5_cell_type_possibilities, Unknown_5_cell_type_possibilities$`cell type` == c("B cells", "B cells memory", "B cells naive"))

B_cells <- subset(B_cells, B_cells$avg_log2FC >= 1)

B_cells <- B_cells$genes

g1 <- FeaturePlot(data, features = B_cells)

ggsave(plot = g1,
       filename = "B_cells.png",
       path = ".../Figures/Cell Type/",
       width = 15,
       height = 15,
       dpi = 300, limitsize = F)


#data <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(data)

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T)

UK6_markers <- FindMarkers(data, only.pos = F, min.pct = 0.2, logfc.threshold = 0.5, assay = "RNA", ident.1 = "16")

UK6_markers$genes <- row.names(UK6_markers)

UK6_marker_list <- PanglaoDB$`official gene symbol`


#export and capitalize data

writexl::write_xlsx(UK6_markers, path = ".../Unknown 6 markers.xlsx")

UK6_marker_use <- Unknown_6_markers %>%  filter(gene_upper %in% UK6_marker_list)

filter_terms <- UK6_marker_use$gene_upper

Pangio_UK <- PanglaoDB %>%  filter(PanglaoDB$`official gene symbol` %in% filter_terms)

Pangio_UK <- as.data.frame(Pangio_UK)
UK6_marker_use <- as.data.frame(UK6_marker_use)


test <- merge(x = Pangio_UK, y = UK6_marker_use, by.x ="official gene symbol", by.y = "gene_upper")




writexl::write_xlsx(test, path = ".../Unknown 6 cell type possibilities.xlsx")
### used this

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T, split.by = "orig.ident")

DefaultAssay(data) <- "RNA"

Myoblasts <- subset(Unknown_6_cell_type_possibilities, Unknown_6_cell_type_possibilities$`cell type` == c("Myoblasts"))

Myoblasts<- subset(Myoblasts, Myoblasts$avg_log2FC >= 1)

Myoblasts <- Myoblasts$genes

g1 <- FeaturePlot(data, features = Myoblasts)

ggsave(plot = g1,
       filename = "Myoblasts.png",
       path = ".../Figures/Cell Type/",
       width = 15,
       height = 15,
       dpi = 300, limitsize = F)
#UK7

#data <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(data)

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T)

UK7_markers <- FindMarkers(data, only.pos = F, min.pct = 0.2, logfc.threshold = 1, assay = "RNA", ident.1 = "18")

UK7_markers$genes <- row.names(UK7_markers)

UK7_marker_list <- PanglaoDB$`official gene symbol`


#export and capitalize data

writexl::write_xlsx(UK7_markers, path = ".../Unknown 7 markers.xlsx")

UK7_marker_use <- Unknown_7_markers %>%  filter(gene_upper %in% UK7_marker_list)

filter_terms <- UK7_marker_use$gene_upper

Pangio_UK <- PanglaoDB %>%  filter(PanglaoDB$`official gene symbol` %in% filter_terms)

Pangio_UK <- as.data.frame(Pangio_UK)
UK7_marker_use <- as.data.frame(UK7_marker_use)


test <- merge(x = Pangio_UK, y = UK7_marker_use, by.x ="official gene symbol", by.y = "gene_upper")




writexl::write_xlsx(test, path = ".../Unknown 7 cell type possibilities.xlsx")
### used this

Idents(data) <- "seurat_clusters"

DimPlot(data, label = T, split.by = "orig.ident")

DefaultAssay(data) <- "RNA"

Adipocytes<- subset(Unknown_7_cell_type_possibilities, Unknown_7_cell_type_possibilities$`cell type` == c("Adipocytes"))

Adipocytes<- subset(Adipocytes, Adipocytes$avg_log2FC >= 1.5)

Adipocytes <- Adipocytes$genes

g1 <- FeaturePlot(data, features = Adipocytes)

ggsave(plot = g1,
       filename = "Adipocytes.png",
       path = ".../Figures/Cell Type/",
       width = 15,
       height = 15,
       dpi = 300, limitsize = F)


# After determining all cell types the clusters were relabeled in the object


all_samples.integrated <- harmony_integrated_Mut_WT_ScaleRNA

DimPlot(all_samples.integrated)

Idents(all_samples.integrated) <- "seurat_clusters"

DimPlot(all_samples.integrated)

new.cluster.ids <- c("1", "15", "11", "5", "16", "17", "4", "9", "12", "18", "14", "2", "13", "3", "6", "7", "19", "8", "10")

names(new.cluster.ids) <- levels(all_samples.integrated)

all_samples.integrated <- RenameIdents(all_samples.integrated, new.cluster.ids)

levels(x=all_samples.integrated)<- c("1","2","3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")

DimPlot(all_samples.integrated, label = T)

new.cluster.ids <- c("Epithelial cells 1", "Epithelial cells 1", "Epithelial cells 1", "T cells", "Myeloid cells", "Myeloid cells", "Myeloid cells", "B cells", "Adipocyte", "Adipocyte", "Myocyte", "Epithelial cells 1", "Epithelial cells 2", "Epithelial cells 2", "Epithelial cells 2", "Myoepithelial cells", "Fibroblast", "Endothelial cells", "Myoblast")

names(new.cluster.ids) <- levels(all_samples.integrated)

all_samples.integrated <- RenameIdents(all_samples.integrated, new.cluster.ids)

DimPlot(all_samples.integrated, label = T)

my_cols <- c("Epithelial cells 1" = "#663c8f", "T cells" = "#9cb241", "Myeloid cells" = "#7080d9", "B cells" = "#c99135", "Adipocyte" = "#c868b1", "Myocyte" = "#5bbc72", "Epithelial cells 2" =  "#b54361", "Myoepithelial cells" = "#43c9b0", "Fibroblast" = "#b74a32", "Endothelial cells" = "#6b8b3a", "Myoblast" = "#b67c46")

g1 <- DimPlot(all_samples.integrated, label = F, cols = my_cols)

ggsave(plot = g1,
       filename = "Umap all cell.png",
       path = ".../Figures/",
       width = 8,
       height = 5,
       dpi = 300, limitsize = F)

all_samples.integrated@meta.data[["cell type"]] <- Idents(all_samples.integrated)

saveRDS(all_samples.integrated, file = ".../RDS/harmony_all_samples_integrated_cell_type_Mut_WT_Scale.rds")


n_cells <- FetchData(all_samples.integrated, 
                     vars = c("ident", "HTO_classification")) %>%
  dplyr::count(ident, HTO_classification) %>%
  tidyr::spread(ident, n)

write_xlsx(n_cells,path = ".../harmony cell count by cell type all.xlsx")


myeloid_cells <- subset(all_samples.integrated, idents = "Myeloid cells")

saveRDS(myeloid_cells, file = ".../RDS/harmony_myeloid_cells.rds")

Epithelial_cell_1 <- subset(all_samples.integrated, idents = "Epithelial cells 1")

saveRDS(Epithelial_cell_1, file = ".../RDS/harmony_Epithelial_cell_1.rds")

Epithelial_cell_2 <- subset(all_samples.integrated, idents = "Epithelial cells 2")

saveRDS(Epithelial_cell_2, file = ".../RDS/harmony_Epithelial_cell_2.rds")

TCell <- subset(all_samples.integrated, idents = "T cells")

saveRDS(TCell, file = ".../RDS/harmony_TCell.rds")
