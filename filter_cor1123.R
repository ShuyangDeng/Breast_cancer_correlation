library(tidyr)
library(dplyr)

## 读取数据
setwd('/Users/dengshuyang/Desktop/Bioinfo2019/Breast_project/')
P <- read.csv('Protein_Breast_2.csv', stringsAsFactors = FALSE, header = TRUE)
R <- read.csv('RNA_Breast_2.csv', stringsAsFactors = FALSE, header = TRUE)

## 进行sample相同的筛选
sample.same <- merge(
  data.frame(P = colnames(P), stringsAsFactors = FALSE), 
  data.frame(R = colnames(R), stringsAsFactors = FALSE), 
  by.x = 'P', by.y = 'R'
)

## 对原始数据筛选
P <- P[, sample.same$P]
R <- R[, sample.same$P]

## 转置
row.names(P) <- P$X
P <- select(P, -X)
P <- t(as.matrix(P))
P <- as.data.frame(P)
P$sample <- row.names(P)

row.names(R) <- R$X
R <- select(R, -X)
R <- t(as.matrix(R))
R <- as.data.frame(R)
R$sample <- row.names(R)

## 进行gene相同的筛选
gene.same <- merge(
  data.frame(P = colnames(P), stringsAsFactors = FALSE), 
  data.frame(R = colnames(R), stringsAsFactors = FALSE), 
  by.x = 'P', by.y = 'R'
)

## gene筛选
P <- P[, gene.same$P]
R <- R[, gene.same$P]

#################
## 先计算相同名称的相关系数-做筛选条件, 减少数据量
gene.same$SAME_CORR <- NA
for (col in gene.same$P) {
  if (col == 'sample') {
    gene.same[gene.same$P == col, ]$SAME_CORR <- 1
  } else {
    gene.same[gene.same$P == col, ]$SAME_CORR <- cor(
      x = P[, col], y = R[, col], use = "complete.obs")
  }
}

## 筛选同样基因之间相关系数 < 0.3 的数据
gene.same <- gene.same %>% filter(abs(SAME_CORR) < 0.3)
#################

## 满足基本的条件之后的数据
P.fiter.col <- P[, gene.same$P]
R.fiter.col <- R[, gene.same$P]

## 计算数据
result <- data.frame(stringsAsFactors = FALSE)
gene.num <- nrow(gene.same)
gene.num <- 100 #try data

for (gene_P in c(1:gene.num)) {
  for (gene_R in c(1:gene.num)) {
    corr <- cor(x = P.fiter.col[, gene_P], y = R.fiter.col[, gene_R], use = "complete.obs")
    ##  判别是否满足条件
    if (gene_P == gene_R) {
      result <- rbind(
        result,
        data.frame(
          gene_P = gene.same$P[gene_P], gene_R = gene.same$P[gene_R], 
          corr = corr, cor = '<0.3', stringsAsFactors = FALSE
        )
      )
    } else if ((gene_P != gene_R) & abs(corr) > 0.6) {
      result <- rbind(
        result,
        data.frame(
          gene_P = gene.same$P[gene_P], gene_R = gene.same$P[gene_R], 
          corr = corr, cor = '>0.6', stringsAsFactors = FALSE
        )
      )
    } else {
      result <- result
    }
  }
}

# corr.result <- result %>% mutate(
#   cor = ifelse(
#     gene_P == gene_R & abs(corr) < 0.3, '<0.3', ifelse(
#       gene_P != gene_R & abs(corr) > 0.6, '>0.6', NA
#     ))
# ) %>% filter(!is.na(cor))