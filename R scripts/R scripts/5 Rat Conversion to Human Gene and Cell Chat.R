##########################################
# Ortholog gene name conversion in seurat
##########################################

# Add the path to rds file

rat1 <- `Final_ harmony_all_samples_integrated_cell_type_Mut_WT`
DefaultAssay(rat1)<-"RNA"

#read in file created (mart_human_rat.csv) downloaded from ensemble v110 accessed on 7/5/2023

# the mart file has column 6= human gene name, col 5= rat gene names
mart1<-mart_export[,c(5,6)]

# create a column with Human genename- rat gene name to remove any duplicates (keep unique pairs)

mart_u <- na.omit(mart1)
mart_u <- unique(mart_u)

ratgenes<-as.vector(unlist(mart_u$`Rat gene name`))
humangenes<-as.vector(unlist(mart1$`Gene name`))

#convert gene names for rat data- NOW in the DATA slot of rds file
for (i in ratgenes){
  pn <-paste0("^",i,"$")
  p=which(grepl(pn,mart_u$`Rat gene name`))
  p1<-p[1]
  f<-mart_u$`Gene name`[p1]
  pn2<-i
  pn3<-paste0("^",pn2,"$")
  k=grep(pn3,rat1@assays$RNA@data@Dimnames[[1]])
  k1<-k[1]
  rat1@assays$RNA@data@Dimnames[[1]][k1]<-as.character(f)
}

#save the rds with gene names changed
saveRDS(rat1,file="Human_genes_data_slot.RDS")

## run CellChat
#### All samples combined
options(stringsAsFactors = FALSE)

DimPlot(Human_genes_data_slot)

Idents(Human_genes_data_slot) <- "cell type"
DimPlot(Human_genes_data_slot)
Human_genes_data_slot_select <- subset(Human_genes_data_slot, idents = c("Epithelial cells 1", "T cells", "Myeloid cells", "Epithelial cells 2", "Myoepithelial cells", "Fibroblast", "Endothelial cells"))


cellchat <- createCellChat(object = Human_genes_data_slot_select, group.by = "ident", datatype = "RNA", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 6)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
test <- netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

setwd(".../Cell Chat/selected populations/")

png(filename = paste("all interactions.png"), width = 8, height = 6, units = "in", res = 300)
test 
dev.off()


setwd(".../Cell Chat/selected populations/all combined/circle/")

# for circle plots
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
}

setwd(".../selected populations/all combined/heatmap/")

# for heatmaps
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
}

saveRDS(cellchat, file = ".../Cell Chat/selected populations/all_geno_seurat_cell_chat_selected.rds")

#### KI

options(stringsAsFactors = FALSE)

DimPlot(Human_genes_data_slot_select)
Idents(Human_genes_data_slot_select) <- "orig.ident"
data <- subset(Human_genes_data_slot_select, idents = "KI")
Idents(data) <- "cell type"
DimPlot(data)

cellchat <- createCellChat(object = data, group.by = "ident", datatype = "RNA", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 6)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
test <- netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

setwd(".../Cell Chat/selected populations/KI/")

png(filename = paste("all interactions.png"), width = 8, height = 6, units = "in", res = 300)
test 
dev.off()


# for circle plots
setwd(".../Cell Chat/selected populations/KI/circle/")
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
}

# for heatmaps

setwd(".../Cell Chat/selected populations/KI/heatmaps/")
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
}

saveRDS(cellchat, file = ".../Cell Chat/selected populations/KI_seurat_cell_chat_selected.rds")

#### KO

options(stringsAsFactors = FALSE)

DimPlot(Human_genes_data_slot_select)
Idents(Human_genes_data_slot_select) <- "orig.ident"
data <- subset(Human_genes_data_slot_select, idents = "KO")
Idents(data) <- "cell type"
DimPlot(data)

cellchat <- createCellChat(object = data, group.by = "ident", datatype = "RNA", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 6)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
test <- netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

setwd(".../Cell Chat/selected populations/KO/")

png(filename = paste("all interactions.png"), width = 8, height = 6, units = "in", res = 300)
test 
dev.off()

# for circle plots
setwd(".../Cell Chat/selected populations/KO/circle/")
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
}

# for heatmaps

setwd(".../Cell Chat/selected populations/KO/heatmaps/")
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
}

saveRDS(cellchat, file = ".../Cell Chat/selected populations/KO_seurat_cell_chat_selected.rds")

#### WT

options(stringsAsFactors = FALSE)

DimPlot(Human_genes_data_slot_select)
Idents(Human_genes_data_slot_select) <- "orig.ident"
data <- subset(Human_genes_data_slot_select, idents = "WT")
Idents(data) <- "cell type"
DimPlot(data)

cellchat <- createCellChat(object = data, group.by = "ident", datatype = "RNA", assay = "RNA")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 6)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
test <-  netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

setwd(".../Cell Chat/selected populations/WT/")

png(filename = paste("all interactions.png"), width = 8, height = 6, units = "in", res = 300)
test 
dev.off()



# for circle plots
setwd(".../Cell Chat/selected populations/WT/circle/")
pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  
}

# for heatmaps

setwd(".../Cell Chat/selected populations/WT/heatmaps/")

#cellchat <- WT_seurat_cell_chat_selected

pathways <- cellchat@netP$pathways
pathways <- as.data.frame(pathways)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (i in 1:nrow(pathways)){
  
  pathways.show <-pathways[i,]
  par(mfrow=c(1,1))
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)
  test <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
  
}

saveRDS(cellchat, file = ".../Cell Chat/selected populations/WT_seurat_cell_chat_selected.rds")


# comparison of all KI to KO

# Define the cell labels to lift up
group.new = levels(KI_seurat_cell_chat_selected@idents)
KI_KO <- liftCellChat(KO_seurat_cell_chat_selected, group.new)


object.list <- list(KI = KI_seurat_cell_chat_selected, KO = KO_seurat_cell_chat_selected)

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

KI_KO@netP$pathways

pathways_KI <- KI_seurat_cell_chat_selected@netP$pathways
pathways_KI <- as.data.frame(pathways_KI)
colnames(pathways_KI)[1] <- "pathways"

pathways_KO <- KO_seurat_cell_chat_selected@netP$pathways
pathways_KO <- as.data.frame(pathways_KO)
colnames(pathways_KO)[1] <- "pathways"

pathways <- pathways_KI %>% filter(pathways %in% pathways_KO$pathways)

setwd(".../Cell Chat/selected populations/Comparison/KI to KO/")



for (i in 1:nrow(pathways)){
  
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)  
  
  pathways.show <- pathways[i,]
  
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
  
  vertex.receiver = seq(1,10) 
  
  par(mfrow = c(1,2), xpd=TRUE)
  
  
  
  (for (a in 1:length(object.list)) {
    
    test2 <- netVisual_aggregate(object.list[[a]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[a]))
    
  })
  
  
  test2
  dev.off()
}
#dev.off()

# comparison of all KI WT

# Define the cell labels to lift up
group.new = levels(WT_seurat_cell_chat_selected@idents)
KI_KO <- liftCellChat(KO_seurat_cell_chat_selected, group.new)


object.list <- list(KI = KI_seurat_cell_chat_selected, WT = WT_seurat_cell_chat_selected)

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


pathways_KI <- KI_seurat_cell_chat_selected@netP$pathways
pathways_KI <- as.data.frame(pathways_KI)
colnames(pathways_KI)[1] <- "pathways"

pathways_WT <- WT_seurat_cell_chat_selected@netP$pathways
pathways_WT <- as.data.frame(pathways_WT)
colnames(pathways_WT)[1] <- "pathways"

pathways <- pathways_KI %>% filter(pathways %in% pathways_WT$pathways)

setwd(".../Cell Chat/selected populations/Comparison/KI to WT/")

for (i in 1:nrow(pathways)){
  
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)  
  
  pathways.show <- pathways[i,]
  
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
  
  vertex.receiver = seq(1,10) 
  
  par(mfrow = c(1,2), xpd=TRUE)
  
  for (a in 1:length(object.list)) {
    
    test2 <- netVisual_aggregate(object.list[[a]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[a]))
    
  }
  
  
  test2
  dev.off()
}
#dev.off()


# comparison of all KO WT

# Define the cell labels to lift up
group.new = levels(WT_seurat_cell_chat_selected@idents)
KI_KO <- liftCellChat(KO_seurat_cell_chat_selected, group.new)

# 
object.list <- list(KO = KO_seurat_cell_chat_selected, WT = WT_seurat_cell_chat_selected)

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


pathways_KO <- KO_seurat_cell_chat_selected@netP$pathways
pathways_KO <- as.data.frame(pathways_KO)
colnames(pathways_KO)[1] <- "pathways"

pathways_WT <- WT_seurat_cell_chat_selected@netP$pathways
pathways_WT <- as.data.frame(pathways_WT)
colnames(pathways_WT)[1] <- "pathways"

pathways <- pathways_KO %>% filter(pathways %in% pathways_WT$pathways)

setwd(".../Cell Chat/selected populations/Comparison/KO to WT/")

for (i in 1:nrow(pathways)){
  
  png(filename = paste(pathways[i,],".png"), width = 8, height = 6, units = "in", res = 300)  
  
  pathways.show <- pathways[i,]
  
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
  
  vertex.receiver = seq(1,10) 
  
  par(mfrow = c(1,2), xpd=TRUE)
  
  for (a in 1:length(object.list)) {
    
    test2 <- netVisual_aggregate(object.list[[a]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[a]))
    
  }
  
  
  test2
  dev.off()
}
#dev.off()

