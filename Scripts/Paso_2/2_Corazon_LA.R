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
metadata <- read.delim("C:/Users/USUARIO/Desktop/Elías/VIU/TFM/Data/Dataset/Corazón/GSE109816_normal_heart_cell_cluster_info.txt", header = TRUE, row.names = 1)
genes.comunes <- intersect(rownames(metadata),colnames(matriz_c))
matriz_fl <- matriz_c[genes.comunes]
rm(matriz_c)
rm(genes.comunes)


# Filtrar la metadata y la matriz de conteo para la condición LA
barcodes.LA <- metadata[grepl("LA", metadata$Condition), ]
genes.LA <- intersect(rownames(barcodes.LA), colnames(matriz_fl))
metadata.LA <- barcodes.LA[rownames(barcodes.LA) %in% genes.LA, ]
matriz_LA <- matriz_fl[genes.LA]
rm(barcodes.LA)
rm(genes.LA)
rm(metadata)

# ** Flujo de trabajo para LA **
# Crear el objeto Seurat y agregar metadata
LA_seurat <- CreateSeuratObject(counts = matriz_LA)
LA_seurat <- AddMetaData(LA_seurat, metadata = metadata.LA)
# Normalización y escalado usando SCTransform
LA_seurat <- SCTransform(LA_seurat)
# Reducción dimensional mediante PCA y UMAP
LA_seurat <- RunPCA(LA_seurat, npcs = 50, verbose = FALSE)
LA_seurat <- RunUMAP(LA_seurat, dims = 1:30)
rm(matriz_LA)
rm(metadata.LA)

# Modificación del nombre de metadata
LA_seurat@meta.data$CellType[grepl("MP", LA_seurat@meta.data$CellType)] <- "Macrophage"
LA_seurat@meta.data$CellType[grepl("CM", LA_seurat@meta.data$CellType)] <- "Cardiac Myocite"
LA_seurat@meta.data$CellType[grepl("SMC", LA_seurat@meta.data$CellType)] <- "Smooth Muscle Cell"
LA_seurat@meta.data$CellType[grepl("FB", LA_seurat@meta.data$CellType)] <- "Cardiac Fibroblast"
LA_seurat@meta.data$CellType[grepl("EC", LA_seurat@meta.data$CellType)] <- "Endothelial Cell"
Idents(LA_seurat) <- LA_seurat@meta.data$CellType
saveRDS(LA_seurat, file =
          "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LA/Flujo de trabajo/seurat.rds")

# **NICHES**
# Inferencia de la comunicación celular usando NICHES 
NICHES <- RunNICHES(LA_seurat,
                    assay = 'SCT', # Usar datos normalizados con SCTransform
                    species = 'human',
                    LR.database = 'fantom5', # Base de datos de ligandos-receptores
                    cell_types = 'CellType', # Tipos celulares 
                    CellToCell = T) # Matriz comunicación celular-celular
NICHES_seurat <- NICHES$CellToCell # Extracción de la matriz de comunicación celular a celular
NICHES_seurat <- ScaleData(NICHES_seurat)
NICHES_seurat <- RunPCA(NICHES_seurat,features = rownames(NICHES_seurat))
# Búsqueda de los ligandos y receptores diferenciados en todos los clústeres 
marcadoresNICHES <- FindAllMarkers(
  object = NICHES_seurat,
  test.use = "wilcox",      
  min.pct = 0.01,             
  only.pos = TRUE      
)
rm(NICHES_seurat)

# **Creación de la Matriz de Interacción (MI)**
resultadoNICHES <- marcadoresNICHES %>%
  separate(col= cluster, sep="—", into=c("Sender","Receiver")) %>%
  separate(col=gene, sep="—", into=c("Ligand","Receptor")) %>%
  select(c(Sender, Receiver,Ligand,Receptor)) %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_"))
resultadoNICHES <- resultadoNICHES[,c('Ligand', 'Receptor', 'Sender', 'Receiver', 'all')]

# Guardar la matriz de interacción
saveRDS(resultadoNICHES,
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LA/Flujo de trabajo/LR_NICHES_LA.rds")
rm(resultadoNICHES)


# **CellChat**

# Preparar el objeto Seurat para análisis con CellChat
LA_seurat@meta.data$celltype <- LA_seurat@meta.data$CellType # Asignar las anotaciones a celltype
Idents(LA_seurat) <- LA_seurat@meta.data$CellType # Definir los identificadores celulares
LA_seurat@assays[["RNA"]] <- LA_seurat@assays[["SCT"]] # Usar datos SCT como RNA
DefaultAssay(LA_seurat) <- "RNA" # Establecer el conjunto de datos por defecto

# Inferencia de la comunicación celular usando CellChat
result <- CCCbank(LA_seurat, method = 'CellChat', species = 'human', 
                   database = 'FANTOM5', extension =FALSE)
rm(LA_seurat)

# **Creación de la Matriz de Interacción (MI)**
result <- result %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_")) %>%
  select(c(Sender,Receiver,Ligand,Receptor, all))

# Guardar la matriz de interacción
saveRDS(result, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LA/Flujo de trabajo/LR_CellChat_LA.rds")
rm(result)