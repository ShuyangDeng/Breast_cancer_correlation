##homework2019.11.13

#function-3
library(stringr)
func.cor3 <- function(x,y){#func.cor
  na1 <- which(is.na(x)==TRUE)
  na2 <- which(is.na(y)==TRUE)
  nas <- union(na1,na2)
  if(length(nas)!=0){
    x <- x[-nas]
    y <- y[-nas]#?????????????????
  }
  
  
  nn <- cor(x,y)
  nn <- round(nn,5)#??????5?????????
  nn <- format(nn, nsmall = 5)#????????????????????????5?????????
  mn <- cor.test(x,y,alternative = "two.sided",method = "pearson",conf.level = 0.95)
  if(mn$p.value<0.001){
    w <-as.character(str_c(nn,"***"))
    return(w)
  }
  else if(mn$p.value<0.01 & mn$p.value>0.001){
    w <-as.character(str_c(nn,"**"))
    return(w)
  }else if(mn$p.value<0.05 & mn$p.value>0.01){
    w <-as.character(str_c(nn,"*"))
    return(w)
  }else if(mn$p.value>0.05){
    w <- as.character(nn)
    return(w)
  }
}

P <- read.csv(file = "/Users/liuce/Desktop/homework/19.11.13/Protein_Breast_2.csv",header = TRUE)
head(P)

P2 <- P[,-1]
rownames(P2) <- P[,1]

nn <- setdiff(colnames(P2),colnames(R2))
ns <- vector();n=1
for (i in 1:length(nn)) {
  nn.i <- nn[i]
  w.i <- which(colnames(P2)==nn.i)
  ns[n] <- w.i
  n=n+1
}
P2 <- P2[,-ns]#去掉多余样品

nk <- colnames(P2)
ns <- order(nk)
P2 <- P2[,ns]##调整样品顺序


P2_1000 <- P2[1:1000,1:77]##前1000行，可以修改

R <- read.csv(file = "/Users/liuce/Desktop/homework/19.11.13/RNA_Breast_2.csv",header = TRUE)
R2 <- R[,-1]
rownames(R2) <- R[,1]

nk <- colnames(R2)
ns <- order(nk)
R2 <- R2[,ns]##调整样品顺序




R2_1000 <- R2[1:1000,]##前1000行，可以修改

D_M <- matrix(rep(NA,3*nrow(R2_1000)*nrow(P2_1000)),ncol =3 )
colnames(D_M) <- c("gene1_RNA","gene_Protein","correlation")

n=1
for (i in 1:nrow(R2_1000)) {
  for (j in 1:nrow(P2_1000)) {
    D_M[n,1] <- rownames(R2_1000)[i]
    D_M[n,2] <- rownames(P2_1000)[j]
    D_M[n,3] <- func.cor3(as.numeric(R2_1000[i,]),as.numeric(P2_1000[j,]))
    n=n+1
  }
}

D_M

