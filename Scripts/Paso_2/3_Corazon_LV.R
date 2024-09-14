# Carga de paquetes necesarios para el análisis
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(CCCbank)
devtools::load_all("C:/Users/USUARIO/Downloads/NICHES")


# **Procesamiento**
# Creación de la matriz de conteo con metadata
matriz_c <- read.csv("C:/Users/USUARIO/Desktop/Elías/VIU/TFM/Data/Dataset/Corazón/GSE109816_normal_heart_umi_matrix.csv", header = TRUE, row.names = 1)
metadata=read.delim("C:/Users/USUARIO/Desktop/Elías/VIU/TFM/Data/Dataset/Corazón/GSE109816_normal_heart_cell_cluster_info.txt", header = TRUE, row.names = 1)
genes.comunes=intersect(rownames(metadata),colnames(matriz_c))
matriz_fl=matriz_c[genes.comunes]
rm(matriz_c)
rm(genes.comunes)


# Filtrar la metadata y la matriz de conteo para la condición LV
barcodes.LV=metadata[grepl("LV", metadata$Condition), ]
genes.LV=intersect(rownames(barcodes.LV), colnames(matriz_fl))
metadata.LV <- barcodes.LV[rownames(barcodes.LV) %in% genes.LV, ]
matriz_LV=matriz_fl[genes.LV]
rm(barcodes.LV)
rm(genes.LV)
rm(metadata)

# ** Flujo de trabajo para LV ** 
# Crear el objeto Seurat y agregar metadata
LV_seurat <- CreateSeuratObject(counts = matriz_LV)
LV_seurat <- AddMetaData(LV_seurat, metadata = metadata.LV)
# Normalización y escalado usando SCTransform
LV_seurat <- SCTransform(LV_seurat)
# Reducción dimensional mediante PCA y UMAP
LV_seurat <- RunPCA(LV_seurat, npcs = 50, verbose = FALSE)
ElbowPlot(LV_seurat, ndims = 50)
LV_seurat <- FindNeighbors(LV_seurat, dims = 1:10)
LV_seurat <- RunUMAP(LV_seurat, dims = 1:10)
rm(matriz_LV)
rm(metadata.LV)

# Modificación del nombre de metadata
LV_seurat@meta.data$CellType[grepl("MP", LV_seurat@meta.data$CellType)] <- "Macrophage"
LV_seurat@meta.data$CellType[grepl("CM", LV_seurat@meta.data$CellType)] <- "Cardiac Myocite"
LV_seurat@meta.data$CellType[grepl("SMC", LV_seurat@meta.data$CellType)] <- "Smooth Muscle Cell"
LV_seurat@meta.data$CellType[grepl("FB", LV_seurat@meta.data$CellType)] <- "Cardiac Fibroblast"
LV_seurat@meta.data$CellType[grepl("EC", LV_seurat@meta.data$CellType)] <- "Endothelial Cell"
Idents(LV_seurat) <- LV_seurat@meta.data$CellType
saveRDS(LV_seurat, file =
          "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LV/Flujo de trabajo/seurat.rds")

# **NICHES**
# Inferencia de la comunicación celular usando NICHES
NICHES <- RunNICHES(LV_seurat,
                    assay = 'SCT',  # Usar datos normalizados con SCTransform
                    species = 'human',
                    LR.database = 'fantom5', # Base de datos de ligandos-receptores
                    cell_types = 'CellType', # Tipos celulares 
                    CellToCell = T) # Matriz comunicación celular-celular
NICHES_seurat <- NICHES$CellToCell  # Extracción de la matriz de comunicación celular a celular
NICHES_seurat <- ScaleData(NICHES_seurat)

# Búsqueda de los ligandos y receptores diferenciados en todos los clústeres  
marcadoresNICHES <- FindAllMarkers(
  object = NICHES_seurat,
  test.use = "wilcox",      
  min.pct = 0.01,             
  only.pos = TRUE      
)
rm(NICHES_seurat)

# **Creación de la Matriz de Interacción (MI)**
resultadoNICHES= marcadoresNICHES %>%
  separate(col= cluster, sep="—", into=c("Sender","Receiver")) %>%
  separate(col=gene, sep="—", into=c("Ligand","Receptor")) %>%
  select(c(Sender, Receiver,Ligand,Receptor)) %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_"))
resultadoNICHES <- resultadoNICHES[,c('Ligand', 'Receptor', 'Sender', 'Receiver', 'all')]

# Guardar la matriz de interacción
saveRDS(resultadoNICHES,
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LV/Flujo de trabajo/LR_NICHES_LV.rds")
rm(resultadoNICHES)

# **CellChat**

# Preparar el objeto Seurat para análisis con CellChat
LV_seurat@meta.data$celltype=LV_seurat@meta.data$CellType # Asignar las anotaciones a celltype
Idents(LV_seurat)=LV_seurat@meta.data$CellType # Definir los identificadores celulares
LV_seurat@assays[["RNA"]]=LV_seurat@assays[["SCT"]] # Usar datos SCT como RNA
DefaultAssay(LV_seurat) <- "RNA" # Establecer el conjunto de datos por defecto

# Inferencia de la comunicación celular usando CellChat
result <- CCCbank(LV_seurat, method = 'CellChat', species = 'human', 
                   database = 'FANTOM5', extension =FALSE)
rm(LV_seurat)

# **Creación de la Matriz de Interacción (MI)**
result= result %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_")) %>%
  select(c(Sender,Receiver,Ligand,Receptor, all))

# Guardar la matriz de interacción
saveRDS(result, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LV/Flujo de trabajo/LR_CellChat_LV.rds")
rm(result)
