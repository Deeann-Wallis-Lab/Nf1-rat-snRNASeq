all_samples.integrated <- harmony_all_samples_integrated_cell_type_Mut_WT

all_samples.integrated <- all_samples.integrated

Idents(all_samples.integrated) <- "orig.ident"

all_samples.integrated_DE_KO_TO_KI <- FindMarkers(all_samples.integrated, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "KI")

all_samples.integrated_DE_KO_TO_KI$genes <- row.names(all_samples.integrated_DE_KO_TO_KI)

all_samples.integrated_DE_KO_TO_WT <- FindMarkers(all_samples.integrated, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "WT")

all_samples.integrated_DE_KO_TO_WT$genes <- row.names(all_samples.integrated_DE_KO_TO_WT)

all_samples.integrated_DE_KI_TO_WT <- FindMarkers(all_samples.integrated, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KI", ident.2 = "WT")

all_samples.integrated_DE_KI_TO_WT$genes <- row.names(all_samples.integrated_DE_KI_TO_WT)

# Epithelial 1

Epithelial_cell_1 <- subset(harmony_all_samples_integrated_cell_type_Mut_WT, idents = "Epithelial cells 1")

Idents(Epithelial_cell_1) <- "orig.ident"

Epithelial_cell_1_DE_KO_TO_KI <- FindMarkers(Epithelial_cell_1, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "KI")

Epithelial_cell_1_DE_KO_TO_KI$genes <- row.names(Epithelial_cell_1_DE_KO_TO_KI)

Epithelial_cell_1_DE_KO_TO_WT <- FindMarkers(Epithelial_cell_1, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "WT")

Epithelial_cell_1_DE_KO_TO_WT$genes <- row.names(Epithelial_cell_1_DE_KO_TO_WT)

Epithelial_cell_1_DE_KI_TO_WT <- FindMarkers(Epithelial_cell_1, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KI", ident.2 = "WT")

Epithelial_cell_1_DE_KI_TO_WT$genes <- row.names(Epithelial_cell_1_DE_KI_TO_WT)

# Epithelial 2

Epithelial_cell_2 <- subset(harmony_all_samples_integrated_cell_type_Mut_WT, idents = "Epithelial cells 2")

Idents(Epithelial_cell_2) <- "orig.ident"

Epithelial_cell_2_DE_KO_TO_KI <- FindMarkers(Epithelial_cell_2, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "KI")

Epithelial_cell_2_DE_KO_TO_KI$genes <- row.names(Epithelial_cell_2_DE_KO_TO_KI)

Epithelial_cell_2_DE_KO_TO_WT <- FindMarkers(Epithelial_cell_2, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "WT")

Epithelial_cell_2_DE_KO_TO_WT$genes <- row.names(Epithelial_cell_2_DE_KO_TO_WT)

Epithelial_cell_2_DE_KI_TO_WT <- FindMarkers(Epithelial_cell_2, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KI", ident.2 = "WT")

Epithelial_cell_2_DE_KI_TO_WT$genes <- row.names(Epithelial_cell_2_DE_KI_TO_WT)

# Myeloid cells

myeloid_cells <- subset(harmony_all_samples_integrated_cell_type_Mut_WT, idents = "Myeloid cells")

Idents(myeloid_cells) <- "orig.ident"

myeloid_cells_DE_KO_TO_KI <- FindMarkers(myeloid_cells, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "KI")

myeloid_cells_DE_KO_TO_KI$genes <- row.names(myeloid_cells_DE_KO_TO_KI)

myeloid_cells_DE_KO_TO_WT <- FindMarkers(myeloid_cells, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "WT")

myeloid_cells_DE_KO_TO_WT$genes <- row.names(myeloid_cells_DE_KO_TO_WT)

myeloid_cells_DE_KI_TO_WT <- FindMarkers(myeloid_cells, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KI", ident.2 = "WT")

myeloid_cells_DE_KI_TO_WT$genes <- row.names(myeloid_cells_DE_KI_TO_WT)

# T cells

TCell <- subset(harmony_all_samples_integrated_cell_type_Mut_WT, idents = "T cells")

Idents(TCell) <- "orig.ident"

TCell_DE_KO_TO_KI <- FindMarkers(TCell, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "KI")

TCell_DE_KO_TO_KI$genes <- row.names(TCell_DE_KO_TO_KI)


TCell_DE_KO_TO_WT <- FindMarkers(TCell, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KO", ident.2 = "WT")

TCell_DE_KO_TO_WT$genes <- row.names(TCell_DE_KO_TO_WT)

TCell_DE_KI_TO_WT <- FindMarkers(TCell, only.pos = F, min.pct = 0.05, logfc.threshold = 0.0, assay = "RNA", ident.1 = "KI", ident.2 = "WT")

TCell_DE_KI_TO_WT$genes <- row.names(TCell_DE_KI_TO_WT)

# make excel or CSV for each dataframe