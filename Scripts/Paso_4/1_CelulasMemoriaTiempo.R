#Gráfica
celulas_ratio <- list()
datasets <- c("PBMC","Corazon","Higado")
numero <- c("3k","6k","8k")
corazon <- c("LA","LV")
metodos <- c("NICHES","CellChat")
ratios <- c(0.9,0.8,0.7,0.6,0.5)

for (data in datasets) {
  celulas_ratio[[data]] <- list()
  if (data == "PBMC") {
    for (i in numero) {
      celulas_ratio[[data]][[i]] <- list()
      ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/", i, "/Flujo de trabajo/", sep = "")
      ruta <- list.files(ruta,pattern = "seurat" ,full.names = TRUE)
      seurat <- readRDS(ruta)
      celulas_totales <- length(colnames(seurat))
      for (ratio in ratios) {
        celulas <- ceiling(celulas_totales * ratio)
        celulas_ratio[[data]][[i]][[as.character(ratio)]] <- celulas
      }
    }
  }
  if (data == "Corazon") {
    for (i in corazon) {
      celulas_ratio[[data]][[i]] <- list()
      ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/", i, "/Flujo de trabajo/", sep = "")
      ruta <- list.files(ruta,pattern = "seurat" ,full.names = TRUE)
      seurat <- readRDS(ruta)
      celulas_totales <- length(colnames(seurat))
      
      for (ratio in ratios) {
        celulas <- ceiling(celulas_totales * ratio)
        celulas_ratio[[data]][[i]][[as.character(ratio)]] <- celulas
      }
    }
  }
  if (data == "Higado") {
    ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/", data, "/Flujo de trabajo/", sep = "")
    ruta <- list.files(ruta,pattern = "seurat" ,full.names = TRUE)
    seurat <- readRDS(ruta)
    celulas_totales <- length(colnames(seurat))
    
    for (ratio in ratios) {
      celulas <- ceiling(celulas_totales * ratio)
      celulas_ratio[[data]][[as.character(ratio)]] <- celulas
    }
  }
}
dir.celulas <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/")
if(!dir.exists(dir.celulas)){
  dir.create(dir.celulas, recursive = TRUE)
}
rm(dir.celulas)
saveRDS(celulas_ratio,
        file="C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/celulas.rds")


# Crear una lista vacía para almacenar los resultados
resultados_tiempo <- list()

# Función para ajustar el tiempo en minutos
ajustar_tiempo <- function(tiempo_segundos) {
  tiempo_minutos <- tiempo_segundos / 60
  minutos <- floor(tiempo_minutos)
  segundos <- (tiempo_minutos - minutos) * 60
  
  # Ajustar si los segundos son >= 60
  if (segundos >= 59.5) {
    minutos <- minutos + 1
    segundos <- 0
  }
  
  print(minutos)
}

for (met in metodos) {
  for (data in datasets) {
    if (data == "PBMC") {
      for (i in numero) {
        for (ratio in ratios) {
          ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Informacion/Tiempo/", sep = "")
          nombre <- paste0(data,"_",i,"_",ratio)
          ruta_a <- list.files(ruta, pattern = nombre, full.names = TRUE)
          tiempo <- readRDS(ruta_a)
          tiempo <- tiempo[["elapsed"]]
          tiempo <- ajustar_tiempo(tiempo)
          resultados_tiempo[[data]][[i]][[met]][[as.character(ratio)]] <- tiempo
        }
      }
    }
  if (data == "Corazon") {
    for (dato in corazon) {
      for (ratio in ratios) {
        ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Informacion/Tiempo/", sep = "")
        nombre <- paste0(data,"_",dato,"_",ratio)
        ruta_a <- list.files(ruta, pattern = nombre, full.names = TRUE)
        tiempo <- readRDS(ruta_a)
        tiempo <- tiempo[["elapsed"]]
        tiempo <- ajustar_tiempo(tiempo)
        resultados_tiempo[[data]][[dato]][[met]][[as.character(ratio)]] <- tiempo
      }
    }
  } 
    if (data == "Higado") {
      for (ratio in ratios) {
        ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Informacion/Tiempo/", sep = "")
        nombre <- paste0(data,"_",ratio)
        ruta_a <- list.files(ruta, pattern = nombre, full.names = TRUE)
        tiempo <- readRDS(ruta_a)
        tiempo <- tiempo[["elapsed"]]
        tiempo <- ajustar_tiempo(tiempo)
        resultados_tiempo[[data]][[met]][[as.character(ratio)]] <- tiempo
      }
    }
  }
}

dir.tiempo <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/")
if(!dir.exists(dir.tiempo)){
  dir.create(dir.tiempo, recursive = TRUE)
}
rm(dir.tiempo)
saveRDS(resultados_tiempo,
        file="C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/tiempo.rds")

#Memoria
memoria_list <- list()
for (met in metodos) {
  for (data in datasets) {
    if (data == "PBMC") {
      for (i in numero) {
        for (ratio in ratios) {
          ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Informacion/Memoria/", sep = "")
          nombre <- paste0(data,"_",i,"_",ratio,".rds")
          ruta <- list.files(ruta, pattern = nombre, full.names = TRUE)
          memoria <- readRDS(ruta)
          memoria_list[[data]][[i]][[met]][[as.character(ratio)]] <- memoria
        }
      }
    }
  if (data == "Corazon") {
    for (dato in corazon) {
      for (ratio in ratios) {
        ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Informacion/Memoria/", sep = "")
        nombre <- paste0(data,"_",dato,"_",ratio,".rds")
        ruta <- list.files(ruta, pattern = nombre, full.names = TRUE)
        memoria <- readRDS(ruta)
        memoria_list[[data]][[dato]][[met]][[as.character(ratio)]] <- memoria
      }
    }
  } 
    if (data == "Higado") {
      for (ratio in ratios) {
        ruta <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Sampling/",met,"/Informacion/Memoria/", sep = "")
        nombre <- paste0(data,"_",ratio,".rds")
        ruta <- list.files(ruta, pattern = nombre, full.names = TRUE)
        memoria <- readRDS(ruta)
        memoria_list[[data]][[met]][[as.character(ratio)]] <- memoria
      }
    }
  } 
}
dir.memoria <- paste("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/")
if(!dir.exists(dir.memoria)){
  dir.create(dir.memoria, recursive = TRUE)
}
rm(dir.memoria)
saveRDS(memoria_list,
        file="C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/memoria.rds")
