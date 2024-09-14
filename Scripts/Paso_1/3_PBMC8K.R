# Carga de paquetes necesarios para el análisis
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(tidyr)
library(dplyr)
library(ROCR)
library(tidyverse)
library(SingleR)
library(CCCbank)
library(ggplot2)
library(data.table)
devtools::load_all("C:/Users/USUARIO/Downloads/NICHES")
source("C:/Users/USUARIO/Desktop/Elías/VIU/TFM/Data/Script/Paso_1/Función_para_calculo_CAGE_y_proteomica.R") # Carga de funciones adicionales

# **Procesamiento**
# Procesamiento de PBMC8K con Seurat

# Carga de PBMC8K a un objeto Seurat
query=Read10X(data.dir = "C:/Users/USUARIO/Downloads/filtered_gene_bc_matrices/GRCh38")
query=CreateSeuratObject(query)

# Anotación de SingleR usando como referencia PBMC3k de Seurat
pbmc3k1=pbmc3k.final
pbmc3k1=UpdateSeuratObject(pbmc3k1)
reference_sce <- as.SingleCellExperiment(pbmc3k1)
query_sce <- as.SingleCellExperiment(query)


# Extracción de las anotaciones del dataset de referencia y anotación
reference_labels <- pbmc3k1$seurat_annotations 
pred <- SingleR(test = query_sce, ref = reference_sce, labels = reference_labels)
query$SingleR.labels <- pred$labels

# Eliminación de células con alto porcentaje de ARN mitocondrial
query[['percent.mt']] <- PercentageFeatureSet(query, pattern = '^MT-')
query <- subset(query, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
rm(reference_sce)
rm(query_sce)
rm(reference_labels)
rm(pred)
rm(pbmc3k1)


# Imputación de datos y reducción dimensional
Idents(query)=query@meta.data$SingleR.labels
imputed <- SeuratWrappers::RunALRA(query)  # Imputación de datos faltantes
imputed <- ScaleData(imputed)
imputed <- RunPCA(imputed,features = rownames(imputed))
imputed <- RunUMAP(imputed,dims = 1:8)

# Guardar el objeto procesado
saveRDS(imputed, file =
          "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/seurat.rds")
rm(query)

# **NICHES**
# Inferencia de la comunicación celular usando NICHES  
result <- RunNICHES(imputed,
                    assay = 'alra', # Uso del conjunto de datos imputado
                    species = 'human',
                    LR.database = 'fantom5', # Base de datos de ligandos-receptores
                    cell_types = 'SingleR.labels', # Tipos celulares 
                    CellToCell = T) # Matriz comunicación celular-celular
result <- result$CellToCell # Extracción de la matriz de comunicación celular a celular
result <- ScaleData(result)
result <- RunPCA(result,features = rownames(result))

# Búsqueda de los ligandos y receptores diferenciados en todos los clústeres  
marcadoresNICHES <- FindAllMarkers(
  object = result,
  test.use = "wilcox",      
  min.pct = 0.01,             
  only.pos = TRUE      
)

# **Creación de la Matriz de Interacción (MI) para análisis posterior**
# Preparar la matriz de interacción
resultadoNICHES= marcadoresNICHES %>%
  separate(col= cluster, sep="—", into=c("Sender","Receiver")) %>%
  separate(col=gene, sep="—", into=c("Ligand","Receptor")) %>%
  select(c(Sender, Receiver,Ligand,Receptor)) %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_"))
rm(marcadoresNICHES)

# Guardar la matriz de interacción
saveRDS(resultadoNICHES, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/LR_NICHES_PBMC8K.rds")


# **Preparación del dataframe para la comparación con los datos de CAGE FANTOM5 y datos proteómicos**
desired_cell_types <- c('B', "NK", 'CD14+ Mono',
                        'CD8 T', "Memory CD4 T", "Naive CD4 T", "DC") # Tipos celulares de interés
imputed1 <- subset(imputed, idents = desired_cell_types) # Subconjunto de tipos celulares deseados
data <- GetAssayData(imputed1, 'data', 'RNA') %>% 
  as.matrix(.) %>% t(.) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(., 'barcodes') # Conversión de datos en formato adecuado
imputed1@meta.data[["celltype"]]=imputed1@active.ident
meta <- imputed1@meta.data %>% .[, 'celltype', drop =FALSE] %>%
  tibble::rownames_to_column(., 'barcodes')
data <- merge(data, meta, by = 'barcodes')
data$barcodes <- NULL
CellType=imputed1@meta.data[["celltype"]]
data <- aggregate(.~celltype, data, mean)
data <- pivot_longer(data, -celltype,    # Cálculo de la media de expresión por tipo celular      
                     names_to = "genes")
data$ct_genes <- paste(data$celltype, data$genes, sep = '_')
data <- data[ -c(1:2)]
rm(imputed1)
rm(CellType)
rm(meta)
rm(desired_cell_types);gc()

# Guardar datos de expresión
saveRDS(data,
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/CAGE_FANTOM5/CAGE_data_PBMC8k.rds")

# **Obtención del LRscore a partir de la matriz de NICHES y el dataframe de expresión**

# Calcular el valor de ligand_value y receptor_value usando los datos de expresión
resultadoNICHES$sl <- paste(resultadoNICHES$Sender, resultadoNICHES$Ligand, sep = '_')  # Crear identificadores únicos para ligandos
resultadoNICHES <- merge(resultadoNICHES, data, by.x = 'sl', by.y = "ct_genes")  # Merge con los datos de expresión
colnames(resultadoNICHES)[7] <- 'ligand_value'  # Renombrar columna
resultadoNICHES$rr <- paste(resultadoNICHES$Receiver, resultadoNICHES$Receptor, sep = '_')  # Crear identificadores únicos para receptores
resultadoNICHES <- merge(resultadoNICHES, data, by.x = 'rr', by.y = "ct_genes")  # Merge con los datos de expresión
colnames(resultadoNICHES)[9] <- 'receptor_value'  # Renombrar columna
resultadoNICHES$LRscore <- resultadoNICHES$ligand_value * resultadoNICHES$receptor_value  # Calcular el score de ligando-receptor
resultadoNICHES <- resultadoNICHES[, c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore', 'all')]  # Seleccionar columnas relevantes


#Carga de CAGEdata
cage_genes=readRDS(
  file = "C:/Users/USUARIO/Desktop/Elías/TFM/CAGE_FANTOM5/CAGE_PBMC.rds")

# **Evaluación de PBMC 8k mediante AUPRC comparando resultados de NICHES con datos de CAGE**

# Obtención de AUPRC CAGE 
result=resultadoNICHES
result$receptor_ct <- lapply(1:dim(result)[[1]], function(i){
  receptor <- result$Receptor[i]
  if(grepl('&', receptor)){
    receptor <- unlist(strsplit(receptor, '&'))
  } 
  receiver <- result$Receiver[i]
  tmp <- paste(receptor, receiver, sep = '_')
  label <- ifelse(all(tmp %in% cage_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(receptor %in% cage_genes$all), FALSE, NA)
  }
  label
}) %>% unlist() # Etiqueta en una nueva columna con TRUE or FALSE si se encuentra el receptor en CAGE data

result$ligand_ct <- lapply(1:dim(result)[[1]], function(i){
  ligand <- result$Ligand[i]
  if(grepl('&', ligand)){
    ligand <- unlist(strsplit(ligand, '&'))
  }
  sender <- result$Sender[i]
  tmp <- paste(ligand, sender, sep = '_')
  label <- ifelse(all(tmp %in% cage_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(ligand %in% cage_genes$all), FALSE, NA)
  }
  label
}) %>% unlist() # Etiqueta en una nueva columna con TRUE or FALSE si se encuentra el ligando en CAGE data

result <- result[which(!is.na(result$receptor_ct)), ] # Filtrar resultados que tengan etiquetas no NA
result <- result[which(!is.na(result$ligand_ct)),] # Filtrar resultados que tengan etiquetas no NA

result$label <- lapply(1:dim(result)[[1]], function(i){ # Etiqueta con TRUE or FALSE
  all(result$ligand_ct[i] & result$receptor_ct[i]) # si ligando y receptor se encuentra en CAGE
}) %>% unlist()

result$sr <- paste(result$Sender, result$Receiver, sep = '_')
result <- result[, c('LRscore', 'label', 'sr')]
result <- split(result, result$sr)

res_index <- lapply(result, function(res){ # Cálculo de AUPRC
  index <- get_evaluate_metrics(as.numeric(res$LRscore), res$label)
  index <- index$perf_metrics
})
res_index <- do.call(rbind, res_index) #Presenta los resultados en una lista
res_index <- as.data.frame(res_index)
res_index=na.omit(res_index) #Elimina todos los resultados que contenga NA
saveRDS(res_index, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/CAGE_NICHES_PBMC8k.rds")

 #Carga de la tabla proteomica
pep_genes= readRDS(file = "C:/Users/USUARIO/Desktop/Elías/TFM/Datos proteómicos/tabla_PROTE_PBMC.rds")

# **Evaluación de PBMC 8k mediante AUPRC comparando resultados de NICHES con el dataset proteómico **

# Obtención de AUPRC proteómico
result=resultadoNICHES
result$receptor_ct <- lapply(1:dim(result)[[1]], function(i){
  receptor <- result$Receptor[i]
  if(grepl('&', receptor)){
    receptor <- unlist(strsplit(receptor, '&'))
  }
  receiver <- result$Receiver[i]
  tmp <- paste(receptor, receiver, sep = '_')
  label <- ifelse(all(tmp %in% pep_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(receptor %in% pep_genes$all), FALSE, NA)
  }
  label
}) %>% unlist() # Etiqueta en una nueva columna con TRUE or FALSE si se encuentra el receptor en datos proteómicos

result$ligand_ct <- lapply(1:dim(result)[[1]], function(i){
  ligand <- result$Ligand[i]
  if(grepl('&', ligand)){
    ligand <- unlist(strsplit(ligand, '&'))
  }
  sender <- result$Sender[i]
  tmp <- paste(ligand, sender, sep = '_')
  label <- ifelse(all(tmp %in% pep_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(ligand %in% pep_genes$all), FALSE, NA)
  }
  label 
}) %>% unlist() # Etiqueta en una nueva columna con TRUE or FALSE si se encuentra el ligando en datos proteómicos

result <- result[which(!is.na(result$receptor_ct)), ] # Etiqueta con TRUE or FALSE
result <- result[which(!is.na(result$ligand_ct)),] # si ligando y receptor se encuentra en datos proteómicos

result$label <- lapply(1:dim(result)[[1]], function(i){
  all(result$ligand_ct[i] & result$receptor_ct[i])
}) %>% unlist()

result$sr <- paste(result$Sender, result$Receiver, sep = '_')
result <- result[, c('LRscore', 'label', 'sr')]
result <- split(result, result$sr)

res_index <- lapply(result, function(res){
  index <- get_evaluate_metrics(as.numeric(res$LRscore), res$label)
  index <- index$perf_metrics
})
res_index <- do.call(rbind, res_index)
res_index <- as.data.frame(res_index)
res_index <- na.omit(res_index)
saveRDS(res_index, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/PROTE_NICHES_PBMC8k.rds")

# **CellChat**

# Preparación de los datos PBMC para CellChat
imputed@meta.data$celltype=imputed@meta.data$SingleR.labels # Asignar las anotaciones a celltype
Idents(imputed)=imputed@meta.data$SingleR.labels # Definir los identificadores celulares
imputed@assays[["RNA"]]=imputed@assays[["alra"]] # Utilizar el conjunto de datos imputado
DefaultAssay(imputed) <- "RNA" # Establecer el conjunto de datos por defecto

# Inferencia de la comunicación celular usando CellChat
result1 <- CCCbank(imputed, method = 'CellChat', species = 'human', 
                   database = 'FANTOM5', extension =FALSE)
rm(imputed)


# **Creación de la Matriz de Interacción (MI) para análisis posterior**

# Preparar la matriz de interacción
result=result1
result= result %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_")) %>%
  select(c(Sender,Receiver,Ligand,Receptor, all))

# Guardar la matriz de interacción
saveRDS(result, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/LR_CellChat_PBMC8K.rds")

# **Obtención del LRscore a partir de la matriz de CellChat y el dataframe de expresión**

# Calcular el valor de ligand_value y receptor_value usando los datos de expresión
result$sl <- paste(result$Sender, result$Ligand, sep = '_')
result <- merge(result, data, by.x = 'sl', by.y = "ct_genes")
colnames(result)[7] <- 'ligand_value'
result$rr <- paste(result$Receiver, result$Receptor, sep = '_')
result <- merge(result, data, by.x = 'rr', by.y = "ct_genes")
colnames(result)[9] <- 'receptor_value'
result$LRscore <- result$ligand_value*result$receptor_value
result <- result[,c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore', 'all')]

# **Evaluación de PBMC 8k mediante AUPRC comparando resultados de CellChat con datos de CAGE**

# Obtención de AUPRC CAGE
result$receptor_ct <- lapply(1:dim(result)[[1]], function(i){
  receptor <- result$Receptor[i]
  if(grepl('&', receptor)){
    receptor <- unlist(strsplit(receptor, '&'))
  } 
  receiver <- result$Receiver[i]
  tmp <- paste(receptor, receiver, sep = '_')
  label <- ifelse(all(tmp %in% cage_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(receptor %in% cage_genes$all), FALSE, NA)
  }
  label
}) %>% unlist()

result$ligand_ct <- lapply(1:dim(result)[[1]], function(i){
  ligand <- result$Ligand[i]
  if(grepl('&', ligand)){
    ligand <- unlist(strsplit(ligand, '&'))
  }
  sender <- result$Sender[i]
  tmp <- paste(ligand, sender, sep = '_')
  label <- ifelse(all(tmp %in% cage_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(ligand %in% cage_genes$all), FALSE, NA)
  }
  label
}) %>% unlist()

result <- result[which(!is.na(result$receptor_ct)), ]
result <- result[which(!is.na(result$ligand_ct)),]

result$label <- lapply(1:dim(result)[[1]], function(i){
  all(result$ligand_ct[i] & result$receptor_ct[i])
}) %>% unlist()

result$sr <- paste(result$Sender, result$Receiver, sep = '_')
result <- result[, c('LRscore', 'label', 'sr')]
result <- split(result, result$sr)

res_index <- lapply(result, function(res){
  index <- get_evaluate_metrics(as.numeric(res$LRscore), res$label)
  index <- index$perf_metrics
})
res_index <- do.call(rbind, res_index)
res_index <- as.data.frame(res_index)
res_index=na.omit(res_index)
rm(cage_genes)
saveRDS(res_index, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/CAGE_CellChat_PBMC8k.rds")

# **Preparación y obtención de métricas de comparación con datos de proteómica**

# Creación de la MI (Matrix de Interacción) para los posteriores cálculos
result=result1
result= result %>%
  mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_")) %>%
  select(c(Sender,Receiver,Ligand,Receptor, all))

# Calcular el valor de ligand_value y receptor_value usando los datos de expresión
result$sl <- paste(result$Sender, result$Ligand, sep = '_')
result <- merge(result, data, by.x = 'sl', by.y = "ct_genes")
colnames(result)[7] <- 'ligand_value'
result$rr <- paste(result$Receiver, result$Receptor, sep = '_')
result <- merge(result, data, by.x = 'rr', by.y = "ct_genes")
colnames(result)[9] <- 'receptor_value'
result$LRscore <- result$ligand_value*result$receptor_value
result <- result[,c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore', 'all')]
rm(data)



# Obtención de AUPRC proteómico
result$receptor_ct <- lapply(1:dim(result)[[1]], function(i){
  receptor <- result$Receptor[i]
  if(grepl('&', receptor)){
    receptor <- unlist(strsplit(receptor, '&'))
  }
  receiver <- result$Receiver[i]
  tmp <- paste(receptor, receiver, sep = '_')
  label <- ifelse(all(tmp %in% pep_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(receptor %in% pep_genes$all), FALSE, NA)
  }
  label
}) %>% unlist()

result$ligand_ct <- lapply(1:dim(result)[[1]], function(i){
  ligand <- result$Ligand[i]
  if(grepl('&', ligand)){
    ligand <- unlist(strsplit(ligand, '&'))
  }
  sender <- result$Sender[i]
  tmp <- paste(ligand, sender, sep = '_')
  label <- ifelse(all(tmp %in% pep_genes$genes_ct), TRUE, FALSE)
  if(!label){
    label <- ifelse(all(ligand %in% pep_genes$all), FALSE, NA)
  }
  label
}) %>% unlist()

result <- result[which(!is.na(result$receptor_ct)), ]
result <- result[which(!is.na(result$ligand_ct)),]

result$label <- lapply(1:dim(result)[[1]], function(i){
  all(result$ligand_ct[i] & result$receptor_ct[i])
}) %>% unlist()

result$sr <- paste(result$Sender, result$Receiver, sep = '_')
result <- result[, c('LRscore', 'label', 'sr')]
result <- split(result, result$sr)

res_index <- lapply(result, function(res){
  index <- get_evaluate_metrics(as.numeric(res$LRscore), res$label)
  index <- index$perf_metrics
})
res_index <- do.call(rbind, res_index)
res_index <- as.data.frame(res_index)
res_index <- na.omit(res_index)
rm(result)
rm(result1)
rm(pep_genes)
saveRDS(res_index, 
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/PROTE_CellChat_PBMC8k.rds")
