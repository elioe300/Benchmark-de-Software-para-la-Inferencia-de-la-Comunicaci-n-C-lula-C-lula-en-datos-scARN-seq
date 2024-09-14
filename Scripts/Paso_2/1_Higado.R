# Carga de paquetes necesarios para el análisis
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(HumanLiver)
devtools::load_all("C:/Users/USUARIO/Downloads/NICHES")
library(CCCbank)

#Carga del dataset (Al ejecutar la función se abrirá una página web, se deberá
#cerrar y aparecerá el objeto HumanLiverSeurat)

viewHumanLiver()
rm(sCVdL)

# **Procesamiento**
# Creación de la matriz de conteo con metadata
HumanLiverSeurat1 <- HumanLiverSeurat
metadata <- read.delim("C:/Users/USUARIO/Desktop/Elías/VIU/TFM/Data/Dataset/Hígado/GSE115469_CellClusterType.txt", header = TRUE, row.names = 1)
HumanLiverSeurat1 <- AddMetaData(HumanLiverSeurat1, metadata = metadata)

# Normalización, escalado de datos, y reducción dimensional
HumanLiverSeurat1 <- SCTransform(HumanLiverSeurat1)
HumanLiverSeurat1 <- RunUMAP(HumanLiverSeurat1, dims = 1:20)
Idents(HumanLiverSeurat1) <- HumanLiverSeurat1@meta.data$CellType

# Guardar el objeto Seurat
saveRDS(HumanLiverSeurat1, file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Higado/Flujo de trabajo/seuratHIGADO.rds")

# **NICHES**
# Inferencia de la comunicación celular usando NICHES 
NICHES <- RunNICHES(HumanLiverSeurat1,
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
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Higado/Flujo de trabajo/LR_NICHES_Higado.rds")
rm(resultadoNICHES)
rm(marcadoresNICHES)


# **CellChat**

# Preparar el objeto Seurat para análisis con CellChat
HumanLiverSeurat1@meta.data$celltype <- HumanLiverSeurat1@meta.data$CellType # Asignar las anotaciones a celltype
Idents(HumanLiverSeurat1) <- HumanLiverSeurat1@meta.data$CellType # Definir los identificadores celulares
HumanLiverSeurat1@assays[["RNA"]] <- HumanLiverSeurat1@assays[["SCT"]] # Usar datos SCT como RNA
DefaultAssay(HumanLiverSeurat1) <- "RNA"  # Establecer el conjunto de datos por defecto
 
# Inferencia de la comunicación celular usando CellChat
result <- CCCbank(HumanLiverSeurat1, method = 'CellChat', species = 'human', 
                   database = 'FANTOM5', extension =FALSE)
rm(HumanLiverSeurat1)
rm(HumanLiverSeurat)


# **Creación de la Matriz de Interacción (MI)**
result <- as.data.frame(result) 
result <- result %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_")) %>%
  select(c(Sender,Receiver,Ligand,Receptor, all))

# Guardar la matriz de interacción
saveRDS(result, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Higado/Flujo de trabajo/LR_CellChat_Higado.rds")
rm(result)