# Cargar los paquetes necesarios
library(dplyr)

# Cargar la función para calcular el índice de Jaccard
CalJaccardIndex <- function(a, b) {
  a <- as.character(a$all)
  b <- as.character(b$all)
  intersection <- length(intersect(unique(a), unique(b)))
  union <- length(unique(c(a, b)))  # Unión de los dos conjuntos
  return(intersection / union)
}


# Cargar la función para calcular índices de Jaccard por método
CalcularIndicesJaccard <- function(ruta_a, ruta_b) {
  
  resultados_jaccard <- list() # Guarda en una lista todos los resultados
  
  for (archivo_a in ruta_a) {
    nombre_a <- basename(archivo_a)
    data_a <- readRDS(archivo_a)
    
    for (archivo_b in ruta_b) {
      nombre_b <- basename(archivo_b)
      nombre_b <- gsub("\\.rds$", "", nombre_b)  # Quitar la extensión .rds
      data_b <- readRDS(archivo_b)
      
      indice_jaccard <- CalJaccardIndex(data_a, data_b)
      resultados_jaccard[[paste0(nombre_a, "_vs_", nombre_b)]] <- indice_jaccard
    }
  }
  
  return(resultados_jaccard)
}

# Carga de la información y el listado para el análisis del índice de Jaccard 
datasets <- c("PBMC","Corazon","Higado")
numero <- c("3k","6k","8k")
dato <- c("LA","LV")
metodos <- c("CellChat","NICHES")
ratios <- c(0.9,0.8,0.7,0.6,0.5)
resultados_robustez <- list()

# ** Análisis del índice de Jaccard entre el análisis sin subsamplear y el sampleado por el ratio **
for (met in metodos) {
for (data in datasets) {
  if (data == "PBMC") { #Condicional para acceder a los directorios PBMC y clasificarlos en la lista
    for (i in numero) {
      for (ratio in ratios) {
      LR <- paste0("^LR_",met) 
      rute <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/", i, "/Flujo de trabajo/", sep = "")
      ruta_a <- list.files(rute, pattern = LR, full.names = TRUE)
      ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Interacciones/", sep = "")
      nombre <- paste0(data,"_",i,"_",ratio,".rds")
      ruta_b <- list.files(ruta, pattern = nombre, full.names = TRUE)
      resultados <- CalcularIndicesJaccard(ruta_a, ruta_b)
      resultados_robustez[[met]][[data]][[i]][[as.character(ratio)]] <- resultados
      rm(resultados)
      }
    }
  }
  if (data == "Corazon") { #Condicional para acceder a los directorios Corazon y y clasificarlos en la lista 
    for (i in dato) {
      for (ratio in ratios) {
      LR <- paste0("^LR_",met)
      ruta_a <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/", i, "/Flujo de trabajo/", sep = "")
      ruta_a <- list.files(ruta_a, pattern = LR, full.names = TRUE)
      ruta_b <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Interacciones/", sep = "")
      nombre <- paste0(data,"_",i,"_",ratio,".rds")
      ruta_b <- list.files(ruta_b, pattern = nombre, full.names = TRUE)
      resultados <- CalcularIndicesJaccard(ruta_a, ruta_b)
      resultados_robustez[[met]][[data]][[i]][[as.character(ratio)]] <- resultados
      rm(resultados)
    }
    }
  }
  if (data == "Higado") { #Condicional para acceder al directorio Higado
    for (ratio in ratios) {
    LR <- paste0("^LR_", met)
    ruta_a <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/Flujo de trabajo/", sep = "")
    ruta_a <- list.files(ruta_a, pattern = LR, full.names = TRUE)
    ruta_b <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Interacciones/", sep = "")
    nombre <- paste0(data,"_",ratio,".rds")
    ruta_b <- list.files(ruta_b, pattern = nombre, full.names = TRUE)
    resultados <- CalcularIndicesJaccard(ruta_a, ruta_b)
    resultados_robustez[[met]][[data]][[as.character(ratio)]] <- resultados
    rm(resultados)
  }
}
}
}
# Crea los directorios y guarda los resultados 
IndiceJaccard <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/")
if(!dir.exists(IndiceJaccard)){
  dir.create(IndiceJaccard, recursive = TRUE)
}
rm(IndiceJaccard)
saveRDS(resultados_robustez,
        file="C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/IndiceJaccard.rds")
