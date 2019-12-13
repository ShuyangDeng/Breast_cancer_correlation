# install.packages('tidyverse')
library(tidyverse)

## 读取数据
Protein_cormat <- read.csv('correlation_data/Protein_cormat.csv', stringsAsFactors = FALSE)
RNA_cormat <- read.csv('correlation_data/RNA_cormat.csv', stringsAsFactors = FALSE)

# ## RNA_cormat 删除 X. 行列数据
# RNA_cormat <- RNA_cormat %>% select(-X.) %>% filter(X != 'X.')

## 数据处理
# gene.r表示行基因, gene.c 表示列基因 最终数据匹配为(gene.r, gene.c)
# 列转行

Protein_cormat.tran <- Protein_cormat %>% 
  gather(key, value, -X) %>%
  select(gene.r = X, gene.c = key, Protein_value = value)
RNA_cormat.tran <- RNA_cormat %>% 
  gather(key, value, -X) %>% 
  select(gene.r = X, gene.c = key, RNA_value = value)

## 数据合并
plot.data <- Protein_cormat.tran %>% 
  inner_join(
    RNA_cormat.tran, 
    by = c("gene.r", "gene.c")
  )

plot.data <- plot.data %>% filter(!is.na(Protein_value)) %>% 
  filter(!is.na(RNA_value))

## 画图
ggplot(plot.data %>% head(1000), aes(x = Protein_value, y = RNA_value)) +
  geom_point() + 
  xlim(-1, 1) + 
  ylim(-1, 1) + 
  xlab("Protein value") + 
  ylab("RNA value") + 
  ggtitle("gene correlation")
