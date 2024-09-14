# Cargar los paquetes necesarios
devtools::load_all("C:/Users/USUARIO/Downloads/NICHES")
library(CCCbank)
library(profvis)
library(Seurat)
library(dplyr)
library(tidyr)


# Carga de la información y el listado para las rutas de las carpetas
datasets <- c("Corazon","Higado", "PBMC")
numero <- c("3k","6k","8k")
corazon <- c("LA","LV")
metodos <- c("NICHES")
ratios <- c(0.9,0.8,0.7,0.6,0.5)
paths <- list()

#Extracción de las rutas de los datasets elegidos

for (data in datasets) {
  if (data == "PBMC") {
    for (num in numero) {
      ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/", num,"/Flujo de trabajo/seurat.rds", sep = "")
      paths[[data]][[num]] <- ruta
    }
  }
  if (data == "Corazon") {
    for (parte in corazon) {
      ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/", parte,"/Flujo de trabajo/seurat.rds", sep = "")
      paths[[data]][[parte]] <- ruta
    }
  }
  if (data == "Higado") {
    ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/Flujo de trabajo/seurat.rds", sep = "")
    paths[[data]] <- ruta
  }
}


# Función que nombra los datasets analizados por su tipo tisular, clasificación y ratio

identificar_dataset <- function(file_path) {
  dataset_name <- if (grepl("PBMC", file_path)) {
    "PBMC"
  } else if (grepl("Corazon", file_path)) {
    "Corazon"
  } else if (grepl("Higado", file_path)) {
    "Higado"
  }
  
  subcategoria <- if (dataset_name == "PBMC") {
    if (grepl("3k", file_path)) {
      "3k"
    } else if (grepl("6k", file_path)) {
      "6k"
    } else if (grepl("8k", file_path)) {
      "8k"
    }
  } else if (dataset_name == "Corazon") {
    if (grepl("LA", file_path)) {
      "LA"
    } else if (grepl("LV", file_path)) {
      "LV"
    }
  } else {
    ""
  }
  
  identifier <- if (subcategoria != "") {
    paste(dataset_name, subcategoria, ratio, sep = "_")
  } else {
    paste(dataset_name, ratio, sep = "_")
  }
  return(identifier)
}


# ** Sampleo y análisis de los datasets con NICHES **

for (file_path in unlist(paths)) {
  for (ratio in ratios) {
    
    identifier <- identificar_dataset(file_path)
    
    # Lee el archivo de resultados
    resultado_data <- readRDS(file_path)

    
    # Preparar el dataset según el ratio a calcular
    if (ratio < 1) {
      total_cells <- ncol(resultado_data)
      num_cells_to_keep <- ceiling(total_cells * ratio)
      set.seed(123)
      cells_to_keep <- sample(colnames(resultado_data), num_cells_to_keep)
      resultado_data <- subset(resultado_data, cells = cells_to_keep)
    }
    
    profvis_output <- profvis({ #Función que mide la memoria consumida durante el análisis
        # Determinar el tipo de dataset y correr RunNICHES
        pbmc6k_8k <- grepl("6k|8k", file_path, ignore.case = TRUE) #Lee la ruta en búsqueda de PBMC 6k, y 8k
        pbmc_3k <- grepl("3k", file_path, ignore.case = TRUE) #Lee la ruta en búsqueda de PBMC3k
        
        if (pbmc6k_8k) { #Al usar diferentes tipos celulares como anotaciones se crea estos 
          execution_time <- system.time({ # condicionales para que se puedan analizar correctamente
          NICHES <- RunNICHES( # Este condicional es para sólo PBMC 6k y 8k
            resultado_data,
            assay = 'alra',
            species = 'human',
            LR.database = 'fantom5',
            cell_types =  'SingleR.labels',
            CellToCell = TRUE
          )
          # Procesar los resultados finales para NICHES
          NICHES <- NICHES$CellToCell
          NICHES <- ScaleData(NICHES)
          NICHES <- RunPCA(NICHES,features = rownames(NICHES))
          resultado_data <- FindAllMarkers(
            object = NICHES,
            test.use = "wilcox",
            min.pct = 0.01,
            only.pos = TRUE)
          resultado_data <- resultado_data %>%
              separate(col = cluster, sep = "—", into = c("Sender", "Receiver")) %>%
              separate(col = gene, sep = "—", into = c("Ligand", "Receptor")) %>%
              select(c(Sender, Receiver, Ligand, Receptor)) %>%
              mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_"))
        })} else if (pbmc_3k) { # Este condicional es para sólo PBMC 3k
          execution_time <- system.time({ 
          NICHES <- RunNICHES(
            resultado_data,
            assay = 'alra',
            species = 'human',
            LR.database = 'fantom5',
            cell_types = 'seurat_annotations',
            CellToCell = TRUE
          )
          # Procesar los resultados finales para NICHES
          NICHES <- NICHES$CellToCell
          NICHES <- ScaleData(NICHES)
          NICHES <- RunPCA(NICHES,features = rownames(NICHES))
          resultado_data <- FindAllMarkers(
            object = NICHES,
            test.use = "wilcox",
            min.pct = 0.01,
            only.pos = TRUE)
          resultado_data <- resultado_data %>%
              separate(col = cluster, sep = "—", into = c("Sender", "Receiver")) %>%
              separate(col = gene, sep = "—", into = c("Ligand", "Receptor")) %>%
              select(c(Sender, Receiver, Ligand, Receptor)) %>%
              mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_"))
          })} else { # Este condicional es para el resto de datasets 
          execution_time <- system.time({ 
          NICHES <- RunNICHES(
            resultado_data,
            assay = 'SCT',
            species = 'human',
            LR.database = 'fantom5',
            cell_types = 'CellType',
            CellToCell = TRUE
          )
          # Procesar los resultados finales para NICHES
          NICHES <- NICHES$CellToCell
          NICHES <- ScaleData(NICHES)
          NICHES <- RunPCA(NICHES,features = rownames(NICHES))
          resultado_data <- FindAllMarkers(
            object = NICHES,
            test.use = "wilcox",
            min.pct = 0.01,
            only.pos = TRUE)
          resultado_data <- resultado_data %>%
              separate(col = cluster, sep = "—", into = c("Sender", "Receiver")) %>%
              separate(col = gene, sep = "—", into = c("Ligand", "Receptor")) %>%
              select(c(Sender, Receiver, Ligand, Receptor)) %>%
              mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_"))
        })}
    })
      # Crea los directorios y guarda los resultados 
      dir.resultados <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data//Resultados/Consistencia/Sampling/NICHES/Interacciones/", sep = "")
      if(!dir.exists(dir.resultados)){
        dir.create(dir.resultados, recursive = TRUE)
      }
      dir.tiempo <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/NICHES/Informacion/Tiempo/", sep = "")
      if(!dir.exists(dir.tiempo)){
        dir.create(dir.tiempo, recursive = TRUE)
      }
      dir.memoria <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/NICHES/Informacion/Memoria/", sep = "")
      if(!dir.exists(dir.memoria)){
        dir.create(dir.memoria, recursive = TRUE)
      }
      profvis_output <- mean(profvis_output[["x"]][["message"]][["prof"]][["memalloc"]])
      # Guardar los datos y los resultados
      saveRDS(profvis_output,
              file = paste(dir.memoria,identifier,".rds", sep = ""))
      saveRDS(resultado_data,
              file = paste(dir.resultados,identifier,".rds", sep = ""))
      saveRDS(execution_time, 
              file = paste(dir.tiempo,identifier,".rds", sep = ""))
    }
}

# ** Sampleo y análisis de los datasets con CellChat **

for (file_path in unlist(paths)) {
  for (ratio in ratios) {
    
    identifier <- identificar_dataset(file_path)
    
    # Lee el archivo de resultados
    resultado_data <- readRDS(file_path)
    
    
    # Preparar el dataset según el ratio a calcular
    if (ratio < 1) {
      total_cells <- ncol(resultado_data)
      num_cells_to_keep <- ceiling(total_cells * ratio)
      set.seed(123)
      cells_to_keep <- sample(colnames(resultado_data), num_cells_to_keep)
      resultado_data <- subset(resultado_data, cells = cells_to_keep)
    }
    
    # Medir el tiempo de ejecución y ejecutar el método
    profvis_output <- profvis({ 
    execution_time <- system.time({
      # Determinar el tipo de dataset y correr RunNICHES
      pbmc6k_8k <- grepl("6k|8k", file_path, ignore.case = TRUE) # Sigue la misma lógica que NICHES
      pbmc_3k <- grepl("3k", file_path, ignore.case = TRUE)
        if (pbmc6k_8k) { 
          resultado_data@meta.data$celltype = resultado_data@meta.data$SingleR.labels
          resultado_data@meta.data$celltype = factor(resultado_data@meta.data$celltype)
          Idents(resultado_data) = droplevels(resultado_data@meta.data$celltype)
          resultado_data@assays[["RNA"]] = resultado_data@assays[["alra"]]
          DefaultAssay(resultado_data) <- "RNA"
          result <- CCCbank(resultado_data, method = 'CellChat', species = 'human', database = 'FANTOM5', extension = FALSE)
        } else if (pbmc_3k) {
          resultado_data@meta.data$celltype = resultado_data@meta.data$seurat_annotations
          resultado_data@meta.data$celltype = factor(resultado_data@meta.data$celltype)
          Idents(resultado_data) = droplevels(resultado_data@meta.data$celltype)
          resultado_data@assays[["RNA"]] = resultado_data@assays[["alra"]]
          DefaultAssay(resultado_data) <- "RNA"
          result <- CCCbank(resultado_data, method = 'CellChat', species = 'human', database = 'FANTOM5', extension = FALSE)
        } else {
          resultado_data@meta.data$celltype = resultado_data@meta.data$CellType
          resultado_data@meta.data$celltype = factor(resultado_data@meta.data$celltype)
          Idents(resultado_data) = droplevels(resultado_data@meta.data$celltype)
          resultado_data@assays[["RNA"]] = resultado_data@assays[["SCT"]]
          DefaultAssay(resultado_data) = "RNA"
          result <- CCCbank(resultado_data, method = 'CellChat', species = 'human', database = 'FANTOM5', extension = FALSE)
        }
      # Procesar los resultados finales para CellChat
        result <- result %>%
          mutate(all = paste(Sender, Ligand, Receiver, Receptor, sep = "_")) %>%
          select(c(Sender, Receiver, Ligand, Receptor, all))
      })
    }) 
    
    # Crea los directorios y guarda los resultados 
    dir.resultados <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/CellChat/Interacciones/", sep = "")
    if(!dir.exists(dir.resultados)){
      dir.create(dir.resultados, recursive = TRUE)
    }
    
    dir.tiempo <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/CellChat/Informacion/Tiempo/", sep = "")
    if(!dir.exists(dir.tiempo)){
      dir.create(dir.tiempo, recursive = TRUE)
    }
    
    dir.memoria <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/CellChat/Informacion/Memoria/", sep = "")
    if(!dir.exists(dir.memoria)){
      dir.create(dir.memoria, recursive = TRUE)
    }
    profvis_output <- mean(profvis_output[["x"]][["message"]][["prof"]][["memalloc"]])
    # Guardar los datos y los resultados
    saveRDS(profvis_output,
            file = paste(dir.memoria,identifier,".rds", sep = ""))
    saveRDS(result,
            file = paste(dir.resultados,identifier,".rds", sep = ""))
    saveRDS(execution_time, 
            file = paste(dir.tiempo,identifier,".rds", sep = ""))
  }
}

