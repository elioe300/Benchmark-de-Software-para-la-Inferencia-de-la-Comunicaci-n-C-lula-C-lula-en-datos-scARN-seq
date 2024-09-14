# Carga de las rutas que contienen las interacciones predecidas de NICHES y CellChat
Higado <- list.files(path = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Higado/Flujo de trabajo/", pattern = "LR_", full.names = TRUE)
CorazonLA <- list.files(path = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LA/Flujo de trabajo/", pattern = "LR_", full.names = TRUE)
CorazonLV <-list.files(path = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Corazon/LV/Flujo de trabajo/", pattern = "LR_", full.names = TRUE)
PBMC3K <- list.files(path = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/3k/Flujo de trabajo/", pattern = "LR_", full.names = TRUE)
PBMC6K <- list.files(path = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/6k/Flujo de trabajo/", pattern = "LR_", full.names = TRUE)
PBMC8K <- list.files(path = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/", pattern = "LR_", full.names = TRUE)

# Listado de rutas de las carpetas
paths <- list(Higado, CorazonLA, CorazonLV, PBMC3K, PBMC6K, PBMC8K)

# Crear una lista para guardar los resultados
SI_results <- list()

# Función para calcular el índice de similitud
CalSI_1 <- function(a, b) {
  a <- as.character(a$all)
  b <- as.character(b$all)
  intersection <- length(intersect(unique(a), unique(b)))
  union <- min(length(unique(a)), length(unique(b)))
  return(intersection / union)
}

# ** Índice de Similitud **
# Bucle para calcular el índice de similitud entre los datasets
for (path in paths) {
  for (i in 1:length(path)) {
    resultado <- path[i]
    resultado_data <- readRDS(resultado)
    resultado_data$LRscore <- NULL
    rownames(resultado_data) <- 1:nrow(resultado_data)
    resultado_data <- resultado_data[order(resultado_data$Ligand), ]
    nombre_resultado <- tools::file_path_sans_ext(basename(resultado))
    
    for (j in (i+1):length(path)) {
      if (j <= length(path)) {
        resultado1 <- path[j]
        resultado1_data <- readRDS(resultado1)
        resultado1_data$LRscore <- NULL
        rownames(resultado1_data) <- 1:nrow(resultado1_data)
        resultado1_data <- resultado1_data[order(resultado1_data$Ligand), ]
        nombre_resultado1 <- tools::file_path_sans_ext(basename(resultado1))
        
        SI_res <- CalSI_1(resultado_data, resultado1_data)
        
        result_name <- paste0(nombre_resultado, "_vs_", nombre_resultado1)
        SI_results[[result_name]] <- SI_res
      }
    }
  }
}
resultado_path <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Concordancia/", sep = "")
if(!dir.exists(resultado_path)){
  dir.create(resultado_path, recursive = TRUE)
}
# Guardar los resultados en un archivo .RDS
saveRDS(SI_results, file = paste0(resultado_path, "/Indice_Similitud.rds"))


# ** Cálculo de ligandos y receptores totales por dataset **

#Elección de datasets
metodos <- c("NICHES", "CellChat")
archivos <- c("PBMC", "Higado", "Corazon")
numero <- c("3k","6k","8k")
Dato <- c("LA", "LV")
LR <- list()

# Bucle para calcular los ligandos y receptores totales entre los datasets
for (archivo in archivos) {
  LR[[archivo]] <- list()
  if (archivo == "PBMC"){
    for (num in numero) {
    LR[[archivo]][[num]] <- list()
      for (met in metodos) {
        LR[[archivo]][[num]][[met]] <- list()
        nombre <- paste0("^LR_",met) 
        ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", archivo, "/", num, "/Flujo de trabajo/", sep = "")
        ruta_ <- list.files(ruta, pattern = nombre, full.names = TRUE)
        matriz <- readRDS(ruta_)
        ligand <- as.vector(matriz$Ligand)
        ligand <- unique(ligand)
        ligand <- sort(ligand)
        receptor <- as.vector(matriz$Receptor)
        receptor <- sort(receptor)
        receptor <- unique(receptor)
        union <- c(ligand,receptor)
        LR[[archivo]][[num]][[met]] <- union
      }
    }
  }
    if (archivo == "Corazon") {
      for (parte in Dato) {
        for (met in metodos) {
          LR[[archivo]][[parte]][[met]] <- list()
          nombre <- paste0("^LR_",met) 
          ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", archivo, "/", parte, "/Flujo de trabajo/", sep = "")
          ruta_ <- list.files(ruta, pattern = nombre, full.names = TRUE)
          matriz <- readRDS(ruta_)
          ligand <- as.vector(matriz$Ligand)
          ligand <- unique(ligand)
          ligand <- sort(ligand)
          receptor <- as.vector(matriz$Receptor)
          receptor <- sort(receptor)
          receptor <- unique(receptor)
          union <- c(ligand,receptor)
          LR[[archivo]][[parte]][[met]] <- union
        }
      }
    }
    if (archivo == "Higado") {
      for (met in metodos) {
      LR[[archivo]][[met]] <- list()
      nombre <- paste0("^LR_",met) 
      ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", archivo, "/Flujo de trabajo/", sep = "")
      ruta_ <- list.files(ruta, pattern = nombre, full.names = TRUE)
      matriz <- readRDS(ruta_)
      ligand <- as.vector(matriz$Ligand)
      ligand <- unique(ligand)
      ligand <- sort(ligand)
      receptor <- as.vector(matriz$Receptor)
      receptor <- sort(receptor)
      receptor <- unique(receptor)
      union <- c(ligand,receptor)
      LR[[archivo]][[met]] <- union
    }
  }
}
resultado_path <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Concordancia/", sep = "")
if(!dir.exists(resultado_path)){
  dir.create(resultado_path, recursive = TRUE)
} 
# Guardar los resultados en un archivo .RDS
saveRDS(LR, file = paste0(resultado_path, "/LigandoReceptorConcordancia.rds"))
