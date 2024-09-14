# Carga de los datasets para graficar 
metodo <- c("CellChat", "NICHES")
PBMC <- c("3k","6k", "8k")
dataset <- c("PBMC", "Corazon", "Higado")
Corazon <- c("LA", "LV")

# Carga de las rutas donde se encuentra la información a graficar
celulas <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/celulas.rds")
memoria <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/memoria.rds")
tiempo <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/tiempo.rds")


# Bucle para obtener la información en un dataframe para graficar
for (met in metodo) { 
  if (met == "CellChat") {
    df_list <- list()
    for (data in dataset) {
      # Crear un dataframe vacío para cada dataset
      df_total <- data.frame()
      
      if (data == "PBMC") {
        for (numero in PBMC) {
          cell <- unlist(celulas[[data]][[numero]])
          memoria_data <- unlist(memoria[[data]][[numero]][[met]])
          tiempo_data <- unlist(tiempo[[data]][[numero]][[met]])
          
          # Crear un dataframe para los resultados actuales
          df <- data.frame(
            Celulas = cell,
            Tiempo = tiempo_data,
            Memoria = memoria_data
          )
          
          # Combinar el dataframe actual con el dataframe total
          df_total <- rbind(df_total, df)
        }
      } 
      if (data == "Corazon") {
        for (numero in Corazon) {
          cell <- unlist(celulas[[data]][[numero]])
          memoria_data <- unlist(memoria[[data]][[numero]][[met]])
          tiempo_data <- unlist(tiempo[[data]][[numero]][[met]])
          
          # Crear un dataframe para los resultados actuales
          df <- data.frame(
            Celulas = cell,
            Tiempo = tiempo_data,
            Memoria = memoria_data
          )
          
          # Combinar el dataframe actual con el dataframe total
          df_total <- rbind(df_total, df)
        }
      }
      if (data == "Higado") {
        cell <- unlist(celulas[[data]])
        memoria_data <- unlist(memoria[[data]][[met]])
        tiempo_data <- unlist(tiempo[[data]][[met]])
        
        # Crear un dataframe para los resultados actuales
        df <- data.frame(
          Celulas = cell,
          Tiempo = tiempo_data,
          Memoria = memoria_data
        )
        # Combinar el dataframe actual con el dataframe total
        df_total <- rbind(df_total, df)
      }
      # Almacenar el dataframe total en la lista por dataset
      df_list[[data]][[met]] <- df_total
    }
    
    CellChat <- do.call(rbind, unlist(df_list, recursive = FALSE)) #Contendrá toda la información CellChat
  }
  if (met == "NICHES") {
    df_list <- list()
    for (data in dataset) {
      # Crear un dataframe vacío para cada dataset
      df_total <- data.frame()
      
      if (data == "PBMC") {
        for (numero in PBMC) {
          cell <- unlist(celulas[[data]][[numero]])
          memoria_data <- unlist(memoria[[data]][[numero]][[met]])
          tiempo_data <- unlist(tiempo[[data]][[numero]][[met]])
          
          # Crear un dataframe para los resultados actuales
          df <- data.frame(
            Celulas = cell,
            Tiempo = tiempo_data,
            Memoria = memoria_data
          )
          
          # Combinar el dataframe actual con el dataframe total
          df_total <- rbind(df_total, df)
        }
      } 
      if (data == "Corazon") {
        for (numero in Corazon) {
          cell <- unlist(celulas[[data]][[numero]])
          memoria_data <- unlist(memoria[[data]][[numero]][[met]])
          tiempo_data <- unlist(tiempo[[data]][[numero]][[met]])
          
          # Crear un dataframe para los resultados actuales
          df <- data.frame(
            Celulas = cell,
            Tiempo = tiempo_data,
            Memoria = memoria_data
          )
          
          # Combinar el dataframe actual con el dataframe total
          df_total <- rbind(df_total, df)
        }
      }
      if (data == "Higado") {
        cell <- unlist(celulas[[data]])
        memoria_data <- unlist(memoria[[data]][[met]])
        tiempo_data <- unlist(tiempo[[data]][[met]])
        
        # Crear un dataframe para los resultados actuales
        df <- data.frame(
          Celulas = cell,
          Tiempo = tiempo_data,
          Memoria = memoria_data
        )
        # Combinar el dataframe actual con el dataframe total
        df_total <- rbind(df_total, df)
      }
      # Almacenar el dataframe total en la lista por dataset
      df_list[[data]][[met]] <- df_total
    }
    
    NICHES <- do.call(rbind, unlist(df_list, recursive = FALSE)) #Contendrá toda la información NICHES
  }
}
  

# Gráfico con puntos y líneas
ggplot(CellChat, aes(x = Celulas, y = Memoria)) +
  geom_point(color = "blue", size = 3) +  # Añadir puntos al gráfico
  labs(title = "A) Consumo de Memoria durante la Ejecución de CellChat",
       x = "Células",
       y = "Mb") +
  theme_minimal()

ggplot(NICHES, aes(x = Celulas, y = Memoria)) +
  geom_point(color = "blue", size = 3) +  # Añadir puntos al gráfico
  labs(title = "B) Consumo de Memoria durante la Ejecución de NICHES",
       x = "Células",
       y = "Mb") +
  theme_minimal()

ggplot(CellChat, aes(x = Celulas, y = Tiempo)) +
  geom_point(color = "blue", size = 3) +  # Añadir puntos al gráfico
  labs(title = "A) Tiempo de Ejecución de CellChat",
       x = "Células",
       y = "Minutos") +
  theme_minimal()

ggplot(NICHES, aes(x = Celulas, y = Tiempo)) +
  geom_point(color = "blue", size = 3) +  # Añadir puntos al gráfico
  labs(title = "B) Tiempo de Ejecución de NICHES",
       x = "Células",
       y = "Minutos") +
  theme_minimal()

NICHES$dataset <- rownames(NICHES)
NICHES <- NICHES %>%
  filter(!grepl("Higado", dataset))
                
ggplot(NICHES, aes(x = Celulas, y = Tiempo)) +
  geom_point(color = "blue", size = 3) +  # Añadir puntos al gráfico
  labs(title = "Tiempo de Ejecución de NICHES sin el dataset Higado",
       x = "Células",
       y = "Minutos") +
  theme_minimal()

  

