library(dplyr)

# 数据导入
setwd('/Users/dengshuyang/Desktop/Bioinfo2019/Breast_project/')
Protein_Breast <- read.csv('Protein_Breast_2.csv', stringsAsFactors = FALSE)
RNA_Breast <- read.csv('RNA_Breast_2.csv', stringsAsFactors = FALSE)
pam50_genes <- read.csv('pam50_genes.csv', stringsAsFactors = FALSE)

# 数据筛选
Protein_Breast <- merge(pam50_genes, Protein_Breast , by.y = 'X', by.x = 'gene')
RNA_Breast <- merge(pam50_genes, RNA_Breast, by.y = 'X', by.x = 'gene', all.x = TRUE)

# 数据转换
row.names(Protein_Breast) <- Protein_Breast$gene
Protein_Breast <- select(Protein_Breast, -gene) 
Protein_Breast <- as.matrix(Protein_Breast)

row.names(RNA_Breast) <- RNA_Breast$gene
RNA_Breast <- select(RNA_Breast, -gene)
RNA_Breast <- as.matrix(RNA_Breast)

# 画图
rc.PB <- rainbow(nrow(Protein_Breast), start = 0, end = .3)
cc.PB <- rainbow(ncol(Protein_Breast), start = 0, end = .3)
hv.PB <- heatmap(
  Protein_Breast, col = cm.colors(256), scale = "column",
  RowSideColors = rc.PB, ColSideColors = cc.PB, margins = c(5,10),
  xlab = " ", ylab =  " ",
  main = "Protein_Breast"
)

rc.RB <- rainbow(nrow(RNA_Breast), start = 0, end = .3)
cc.RB <- rainbow(ncol(RNA_Breast), start = 0, end = .3)
hv.RB <- heatmap(
  RNA_Breast, col = cm.colors(256), scale = "column",
  RowSideColors = rc.RB, ColSideColors = cc.RB, margins = c(5,10),
  xlab = " ", ylab =  " ",
  main = "RNA_Breast"
)

# 第三幅图
# =====================
Protein_Breast <- read.csv('Protein_Breast_2.csv', stringsAsFactors = FALSE)
RNA_Breast <- read.csv('RNA_Breast_2.csv', stringsAsFactors = FALSE)
pam50_genes <- read.csv('pam50_genes.csv', stringsAsFactors = FALSE)

colnames(Protein_Breast) <- paste0(colnames(Protein_Breast), '_pro')
colnames(RNA_Breast) <- paste0(colnames(RNA_Breast), '_rna')

data <- merge(pam50_genes, Protein_Breast, by.x = 'gene', by.y = 'X_pro', all.x = TRUE)
data <- merge(data, RNA_Breast, by.x = 'gene', by.y = 'X_rna', all.x = TRUE)

row.names(data) <- data$gene
data <- select(data, -gene)
data <- as.matrix(data)

rc <- rainbow(nrow(data), start = 0, end = .3)
cc <- rainbow(ncol(data), start = 0, end = .3)
hv <- heatmap(
  data, col = cm.colors(256), scale = "column",
  RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
  xlab = " ", ylab =  " ",
  main = "data"
)