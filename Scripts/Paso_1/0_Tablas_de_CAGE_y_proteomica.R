# Carga de paquetes necesarios para el análisis
library(data.table)  
library(tibble)     
library(tidyr)       
library(dplyr)      
source("C:/Users/USUARIO/Desktop/Elías/VIU/TFM/Data/Script/Paso_1/Función_para_calculo_CAGE_y_proteomica.R")

# Función para cargar el script que contiene la función de combinación ponderada
celltype <- c('CD14+ Monocytes', 'CD8+ T cells',
              'CD4+CD25-CD45RA- memory conventional T cells',
              'CD4+CD25-CD45RA+ naive conventional T cells',
              'Dendritic Monocyte Immature derived',
              'NK cells','CD19+ B cells'
)

# Carga de la tabla de datos CAGE (Cap Analysis of Gene Expression)
cage_all <- data.table::fread("C:/Users/USUARIO/Desktop/Elías/TFM/CAGE_FANTOM5/ExpressionGenes.txt", header = TRUE, sep = '\t')
cage_all <- as.data.frame(cage_all)
cage_all <- cage_all[ , c(1, which(colnames(cage_all) %in% celltype))]
cage_all <- tibble::column_to_rownames(cage_all, 'ApprovedSymbol')

# Renombrar las columnas con nombres de tipos celulares
colnames(cage_all) <- c('CD14+ Mono', 'B', 'Naive CD4 T',
                        'Memory CD4 T','CD8 T','DC','NK'
)

# Para cada tipo celular, seleccionamos los genes cuya expresión es mayor a 10
# y los almacenamos en una lista
cage_all <- as.matrix(cage_all)
cage <- lapply(1:dim(cage_all)[[2]], function(i){
  x <- cage_all[,i]
  genes <- rownames(cage_all)[which(x>10)]
  paste(genes, colnames(cage_all)[i], sep = '_')
})
cage <- unlist(cage)
cage_genes <- list(all = rownames(cage_all), genes_ct = cage)
rm(cage_all)
rm(celltype)
saveRDS(cage_genes,
        file = "C:/Users/USUARIO/Desktop/Elías/TFM/CAGE_FANTOM5/CAGE_PBMC.rds")


# Carga de la tabla proteómica
proteinGroup <- fread("C:/Users/USUARIO/Desktop/Elías/TFM/Datos proteómicos/proteinGroups.txt")
pep <- proteinGroup[, c(7,549:816)]
rm(proteinGroup)
pep <- pep[which(pep$Gene.names != ''), ]
  colnames(pep) <- gsub('Peptides.', '', colnames(pep))
  pep <- as.data.frame(pep)
  pep <- pep[, c(1, which(grepl('steady-state', colnames(pep))))]
  pep <- pep[, which(!grepl('Library_single', colnames(pep)))]
  colnames(pep) <- gsub(' \\(steady-state\\)', '', colnames(pep))
  
  colnames(pep)[1] <- 'genes'
  pep <- tidyr::separate_rows(pep, genes, sep = ';')
  pep <- aggregate(.~genes, data = pep, FUN=mean) %>%
    tibble::column_to_rownames(., 'genes') %>% 
    as.matrix() %>% t() %>% as.data.frame()
  
  pep$celltype <- gsub('_[0-9]+', '', rownames(pep))
  pep <- aggregate(.~celltype, data = pep, FUN = mean) %>%
    tibble::column_to_rownames(., 'celltype') %>%
    as.matrix() %>% t() %>% as.data.frame()
  
  celltype <- c('Unique.peptides.B.memory', 'Unique.peptides.B.naive', 
                'Unique.peptides.B.plasma', 'Unique.peptides.T8.naive',
                'Unique.peptides.T8.EMRA', 'Unique.peptides.T8.EM',
                'Unique.peptides.T8.CM', 'Unique.peptides.NK.dim',
                'Unique.peptides.NK.brigh', 'Unique.peptides.MO.nonclassical',
                'Unique.peptides.MO.intermediate', 'Unique.peptides.MO.classical',
                'Unique.peptides.T4.naive',
                'Unique.peptides.T4.EMRA', 'Unique.peptides.T4.EM', 
                'Unique.peptides.T4.CM', 'Unique.peptides.pDC','Unique.peptides.mDC'
  )
  pep <- pep[ , which(colnames(pep) %in% celltype)]
  
  pep <- as.matrix(pep) %>% t() %>% as.data.frame()
  pep$celltype <- c('B', 'B', 'B', 'DC', 'CD14+ Mono',
                    'CD14+ Mono', 'CD14+ Mono', 'NK', 'DC', 'Memory CD4 T',
                    'Memory CD4 T','Memory CD4 T', 'Naive CD4 T', 'CD8T', 'CD8T', 'CD8T', 'CD8 T')
  pep <- aggregate(.~celltype, data = pep, FUN = mean) %>%
    tibble::column_to_rownames(., 'celltype') %>%
    as.matrix() %>% t() %>% as.data.frame()
  
  pep_genes <- lapply(1:dim(pep)[[2]], function(i){
    x <- pep[,i]
    genes <- rownames(pep)[which(x>=2)]
    paste(genes, colnames(pep)[i], sep = '_')
  })
  pep_genes <- unlist(pep_genes)
  pep_genes <- list(all = rownames(pep), genes_ct = pep_genes)
  rm(pep)
  rm(celltype)
  saveRDS(pep_genes, 
          file = "C:/Users/USUARIO/Desktop/Elías/TFM/Datos proteómicos/tabla_PROTE_PBMC.rds")