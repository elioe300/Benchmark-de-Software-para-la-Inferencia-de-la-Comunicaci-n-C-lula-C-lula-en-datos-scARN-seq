# Cargar los paquetes necesarios
library(dplyr)
library(ggplot2)

# Carga de los datos
resultados <- readRDS(file = "C:/Users/USUARIO/Desktop/Elías/TFM/Data/Resultados/Consistencia/Datos/IndiceJaccard.rds")

# Inicializar el dataframe que contendrá los resultados de los sampleos 
df_total <- data.frame()
# Listar los datasets de interés
metodos <- c("NICHES", "CellChat")

# Extracción de la información de los datasets Corazon
ratios=c(0.9,0.8,0.7,0.6,0.5)
for(met in metodos) {
  dataset = resultados[[met]]
  for (data in names(dataset)){
    if (data == "Corazon"){
      numero = dataset[[data]]
      for (num in names(numero)) {
        lista = numero[[num]]
        for (ratio in ratios) {
          resultado = lista[[as.character(ratio)]]
          resultado = unlist(resultado)
          # Crear un data.frame común para ambos métodos
          df <- data.frame(
            Dataset = paste0(data,"_",num),  # Columna para indicar el dataset
            Metodo = met,     # Columna para indicar el método
            Jaccard = resultado, # Valores Jaccard
            ratio = as.character(ratio) # Columna para indicar el ratio
          )
          # Acumular resultados en df_total usando bind_rows
          df_total <- bind_rows(df_total, df)
        }
      }
    }
  }
}
# Extracción de la información de los datasets PBMC
for(met in metodos) {
  dataset = resultados[[met]]
  for (data in names(dataset)){
    if (data == "PBMC"){
      numero = dataset[[data]]
      for (num in names(numero)) {
        lista = numero[[num]]
        for (ratio in ratios) {
          resultado = lista[[as.character(ratio)]]
          resultado = unlist(resultado)
          # Crear un data.frame común para ambos métodos
          df <- data.frame(
            Dataset = paste0(data,"_",num), # Columna para indicar el dataset
            Metodo = met,    # Columna para indicar el método
            Jaccard = resultado, # Valores Jaccard
            ratio = as.character(ratio) # Columna para indicar el ratio
          )
          # Acumular resultados en df_total usando bind_rows
          df_total <- bind_rows(df_total, df)
        }
      }
    }
  }
}
# Extracción de la información del dataset Higado
for(met in metodos) {
  dataset = resultados[[met]]
  for (data in names(dataset)){
    if (data == "Higado"){
      numero = dataset[[data]]
        for (ratio in ratios) {
          resultado = numero[[as.character(ratio)]]
          resultado = unlist(resultado)
          # Crear un data.frame común para ambos métodos
          df <- data.frame(
            Dataset = paste0(data), # Columna para indicar el dataset
            Metodo = met,     # Columna para indicar el método
            Jaccard = resultado, # Valores Jaccard
            ratio = as.character(ratio) # Columna para indicar el ratio
          )
          # Acumular resultados en df_total usando bind_rows
          df_total <- bind_rows(df_total, df)
        }
      }
    }
}

#CellChat
CellChat <- df_total[which(df_total$Metodo == "CellChat"),]

p1 <- ggplot(CellChat, aes(x=" ",y = Jaccard, fill = ratio)) +
  geom_boxplot() +
  labs(title = "A) Media por ratio CellChat",
       y = "Indice Jaccard") +
  theme_minimal()    +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16)    
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    hjust = 0.5,   # Ajustar la alineación horizontal
    position = position_dodge(width = 0.75),  # Ajuste para evitar solapamientos
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.25),  # Calcular el primer cuartil (Q1)
    geom = "text", 
    aes(label = paste("Q1:", round(..y.., 3))), 
    vjust = 1.5, 
    hjust = 0.5, 
    position = position_dodge(width = 0.75),
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.75),  # Calcular el tercer cuartil (Q3)
    geom = "text", 
    aes(label = paste("Q3:", round(..y.., 3))), 
    vjust = -1.5, 
    hjust = 0.5, 
    position = position_dodge(width = 0.75),
    color = "black"
  )

CellChat <- ggplot_build(p1)
CellChat <- CellChat$data[[1]]
CellChat <- CellChat[,c(2:7)]


#NICHES
NICHES <- df_total[which(df_total$Metodo == "NICHES"),]

p2 <- ggplot(NICHES, aes(x="",y = Jaccard, fill = ratio)) +
  geom_boxplot() +
  labs(title = "B) Media por ratio NICHES",
       y = "Indice Jaccard") +
  theme_minimal()    +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16)    
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    hjust = 0.5,   # Ajustar la alineación horizontal
    position = position_dodge(width = 0.75),  # Ajuste para evitar solapamientos
    color = "black"
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    hjust = 0.5,   # Ajustar la alineación horizontal
    position = position_dodge(width = 0.75),  # Ajuste para evitar solapamientos
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.25),  # Calcular el primer cuartil (Q1)
    geom = "text", 
    aes(label = paste("Q1:", round(..y.., 3))), 
    vjust = 1.5, 
    hjust = 0.5, 
    position = position_dodge(width = 0.75),
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.75),  # Calcular el tercer cuartil (Q3)
    geom = "text", 
    aes(label = paste("Q3:", round(..y.., 3))), 
    vjust = -1.5, 
    hjust = 0.5, 
    position = position_dodge(width = 0.75),
    color = "black"
  )

NICHES <- ggplot_build(p2)
NICHES <- NICHES$data[[1]]
NICHES <- NICHES[,c(2:7)]