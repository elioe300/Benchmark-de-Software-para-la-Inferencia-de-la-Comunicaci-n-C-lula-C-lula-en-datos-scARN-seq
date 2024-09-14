# Carga de paquetes y datos para obtener las gráficas de AUPRC proteómico
library(ggplot2)
PROTE_CellChat_PBMC3k = readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/3k/Flujo de trabajo/PROTE_CellChat_PBMC3k.rds")
PROTE_NICHES_PBMC3k <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/3k/Flujo de trabajo/PROTE_NICHES_PBMC3k.rds")
PROTE_CellChat_PBMC6k <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/6k/Flujo de trabajo/PROTE_CellChat_PBMC6k.rds")
PROTE_NICHES_PBMC6k <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/6k/Flujo de trabajo/PROTE_NICHES_PBMC6k.rds")
PROTE_CellChat_PBMC8k <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/PROTE_CellChat_PBMC8k.rds")
PROTE_NICHES_PBMC8k <- readRDS("C:/Users/USUARIO/Desktop/Elías/TFM/Data/PBMC/8k/Flujo de trabajo/PROTE_NICHES_PBMC8k.rds")

#Flujo de trabajo para obtener AUPRC CAGE NICHES
# Crear la lista de resultados de NICHES
datasets <- list(
  NICHES_3k = PROTE_NICHES_PBMC3k$AUPRC,
  NICHES_6k = PROTE_NICHES_PBMC6k$AUPRC,
  NICHES_8k = PROTE_NICHES_PBMC8k$AUPRC
)

# Convertir la lista en un dataframe adecuado para ggplot
result_df <- do.call(rbind, lapply(names(datasets), function(name) {
  data.frame(
    Dataset = name,
    AUPRC = as.numeric(datasets[[name]])
  )
}))

# Crear el boxplot con ggplot2 y colores diferenciados por dataset
p1 <- ggplot(result_df, aes(x = Dataset, y = AUPRC, fill = Dataset)) +
  geom_boxplot() +
  labs(title = "B) Boxplot AURPC Proteómico por Dataset",
       x = "Dataset",
       y = "AUPRC") +
  scale_fill_manual(values = c("NICHES_3k" = "lightblue", "NICHES_6k" = "lightgreen", "NICHES_8k" = "lightcoral")) +
theme_minimal()  +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16)    
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.25),  # Calcular el primer cuartil (Q1)
    geom = "text", 
    aes(label = paste("Q1:", round(..y.., 3))), 
    vjust = 1.5,
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.75),  # Calcular el tercer cuartil (Q3)
    geom = "text", 
    aes(label = paste("Q3:", round(..y.., 3))), 
    vjust = -1.5,
    color = "black"
  )

# Obtención de datos a partir de la gráfica Boxplot AUPRC proteómico por Dataset
datos_NICHES_prote <- ggplot_build(p1)
valores_atipicos_NICHES <- datos_NICHES_prote[["data"]][[1]][["outliers"]]
datos_NICHES_prote <- datos_NICHES_prote$data[[1]]
datos_NICHES_prote <- datos_NICHES_prote[,c(2:7)]


# Crear el boxplot con ggplot2 para graficar la media de AUPRC proteómico NICHES
p2 <- ggplot(result_df, aes(x="", y = AUPRC)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "B) Media AUPRC proteómico en NICHES",
       y = "AUPRC") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16)    
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.25),  # Calcular el primer cuartil (Q1)
    geom = "text", 
    aes(label = paste("Q1:", round(..y.., 3))), 
    vjust = 1.5,
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.75),  # Calcular el tercer cuartil (Q3)
    geom = "text", 
    aes(label = paste("Q3:", round(..y.., 3))), 
    vjust = -1.5,
    color = "black"
  )

# Obtención de datos a partir de la gráfica Boxplot media de AUPRC proteómico en NICHES
media_Proteómico_NICHES <- ggplot_build(p2)
media_Proteómico_NICHES <- media_Proteómico_NICHES$data[[1]]
media_Proteómico_NICHES <- media_Proteómico_NICHES[,c(2:6)]


#Flujo de trabajo para obtener AUPRC Proteómico CellChat
# Crear la lista de resultados de CellChat
datasets <- list(
  CellChat_3K = PROTE_CellChat_PBMC3k$AUPRC,
  CellChat_6K = PROTE_CellChat_PBMC6k$AUPRC,
  CellChat_8K = PROTE_CellChat_PBMC8k$AUPRC
)

# Convertir la lista en un dataframe adecuado para ggplot
result_df <- do.call(rbind, lapply(names(datasets), function(name) {
  data.frame(
    Dataset = name,
    AUPRC = as.numeric(datasets[[name]])
  )
}))

# Crear el boxplot con ggplot2 y colores diferenciados por dataset
p3 <- ggplot(result_df, aes(x = Dataset, y = AUPRC, fill = Dataset)) +
  geom_boxplot() +
  labs(title = "A) Boxplot AURPC Proteómico por Dataset",
       x = "Dataset",
       y = "AUPRC") +
  scale_fill_manual(values = c("CellChat_3K" = "lightblue", "CellChat_6K" = "lightgreen", "CellChat_8K" = "lightcoral")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16)    
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.25),  # Calcular el primer cuartil (Q1)
    geom = "text", 
    aes(label = paste("Q1:", round(..y.., 3))), 
    vjust = 1.5,
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.75),  # Calcular el tercer cuartil (Q3)
    geom = "text", 
    aes(label = paste("Q3:", round(..y.., 3))), 
    vjust = -1.5,
    color = "black"
  )

# Obtención de datos a partir de la gráfica Boxplot AUPRC Proteómico por Dataset
datos_CellChat_prote <- ggplot_build(p3)
valores_atipicos_CellChat <- datos_CellChat_prote[["data"]][[1]][["outliers"]]
datos_CellChat_prote <- datos_CellChat_prote$data[[1]]
datos_CellChat_prote <- datos_CellChat_prote[,c(2:7)]

# Crear el boxplot con ggplot2 para graficar la media de AUPRC Proteómico CellChat
p4 <- ggplot(result_df, aes(x = "",y = AUPRC)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "A) Media AUPRC proteómico en CellChat",
       y = "AUPRC") +
  theme_minimal()  +
  theme(
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    plot.title = element_text(size = 16)    
  ) + stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 3)),  # Mostrar la mediana con 3 decimales
    vjust = -0.5, 
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.25),  # Calcular el primer cuartil (Q1)
    geom = "text", 
    aes(label = paste("Q1:", round(..y.., 3))), 
    vjust = 1.5,
    color = "black"
  ) +
  stat_summary(
    fun = function(y) quantile(y, 0.75),  # Calcular el tercer cuartil (Q3)
    geom = "text", 
    aes(label = paste("Q3:", round(..y.., 3))), 
    vjust = -1.5,
    color = "black"
  )

# Obtención de datos a partir de la gráfica Boxplot media de AUPRC Proteómico en CellChat
media_Proteómico_CellChat <- ggplot_build(p4)
media_Proteómico_CellChat <- media_Proteómico_CellChat$data[[1]]
media_Proteómico_CellChat <- media_Proteómico_CellChat[,c(2:6)]