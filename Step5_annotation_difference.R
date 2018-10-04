rm( list = ls() )
load( "./g2t.Rdata" )
load( "./gset.Rdata" )
library( "GEOquery" )
library("hgu133plus2.db")
library( "plyr" )
library ("VennDiagram")

## 用不同的方法注释探针的lincRNA
## 1.bioconductor包hgu133plus2.db
{
  ids1 = toTable(hgu133plus2SYMBOL)
  ids1$type = g2t[match(ids1[,2], g2t$gene),2]
  ids1 = na.omit(ids1)
  ids1 = ids1[ids1[,3]=='lincRNA',]
}

dim(ids1)
length(unique(ids1[,2]))
table(sort(table(ids1[,2])))


## 2.getGEO, soft
gpl570 <- getGEO('GPL570', destdir = ".")
## save(gpl570, file = 'gpl570.Rdata')
## 加载gpl570
load(file = 'gpl570.Rdata')

GPL = gpl570@dataTable@table
colnames(GPL)
GPL[1:2,1:12]

ids2 = GPL[,c(1,11)]

ids=ids2
{
  b<-strsplit(as.character(ids[,2]), " /// ")
  d<-mapply(cbind, ids[,1], b)
  df <- ldply (d, data.frame)
  ids = df[,2:3]
  ids$type = g2t[match(ids[,2], g2t$gene),2]
  ids = na.omit(ids)
}
ids2 = ids
dim(ids2)
table(sort(table(ids2$type)))


## 3.Supplementary Table 1
ids3 = read.csv('20226-290534-1-SP.csv', stringsAsFactors = F)
ids3 = ids3[,c(1,4)]
colnames(ids3) = c("id", "gene")

ids=ids3
{
  b<-strsplit(as.character(ids[,2]), " /// ")
  d<-mapply(cbind, ids[,1], b)
  df <- ldply (d, data.frame)
  ids = df[,2:3]
  ids$type = g2t[match(ids[,2], g2t$gene),2]
  ids = na.omit(ids)
}
ids3 = ids
dim(ids3)
table(sort(table(ids3$type)))


## 4.matrix的注释信息
ids4 <- gset@featureData@data[c("ID", "Gene Symbol")]

ids=ids4
{
  b<-strsplit(as.character(ids[,2]), " /// ")
  d<-mapply(cbind, ids[,1], b)
  df <- ldply (d, data.frame)
  ids = df[,2:3]
  ids$type = g2t[match(ids[,2], g2t$gene),2]
  ids = na.omit(ids)
}
ids4 = ids
dim(ids4)
table(sort(table(ids4$type)))

## 不同注释方法比较
colnames(ids1) = c('ID', 'hgu133plus2', 'typ1')
colnames(ids2) = c('ID', 'soft', 'typ2')
colnames(ids3) = c('ID', 'Supplementary', 'typ3')
colnames(ids4) = c('ID', 'matrix', 'typ4')
ID1A2 = merge(ids1, ids2, all=TRUE) 
ID2A3 = merge(ID1A2, ids3, all=TRUE) 
ID3A4 = merge(ID2A3, ids4, all=TRUE) 

A = na.omit(ID3A4[,2])
B = na.omit(ID3A4[,4])
C = na.omit(ID3A4[,6])
D = na.omit(ID3A4[,8])

## Venn Diagram
venn.plot <- venn.diagram(x= list(A = A, B = B, C = C, D = D), 
             filename = "DIFF.png", height = 450, width = 450,
             resolution =300, imagetype="png", col="transparent",
             fill=c("cornflowerblue","green","yellow","darkorchid1"),alpha = 0.50, cex=0.45, cat.cex=0.45)

