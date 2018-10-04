rm( list = ls() )
options( stringsAsFactors = F )

load( './gset.Rdata' )
library( "GEOquery" )
## 取表达矩阵和样本信息表
{
  gset = gset[[1]]
  exprSet = exprs( gset )
  pdata = pData( gset )
  chl = length( colnames( pdata ) )
  group_list = as.character( pdata[, chl] )
}
dim( exprSet )
exprSet[ 1:5, 1:5 ]
table( group_list )

## 取出研究样本
{
  n_expr = exprSet[ , grep( "healthy",         group_list )]
  g_expr = exprSet[ , grep( "normal|GBM",      group_list )]
  g_expr = g_expr[,-5]
  a_expr = exprSet[ , grep( "fetal",           group_list )]
  exprSet = cbind( n_expr, g_expr, a_expr )
}


## 样本分组
{
  group_list = c(rep( 'normal', ncol( n_expr ) ), 
                 rep( 'gbm',    ncol( g_expr ) +1 ) )
}

dim( exprSet )
exprSet[ 1:5, 1:5 ]
table( group_list )

save( exprSet, group_list, file = 'exprSet_by_group.Rdata')


## 筛选探针
GPL = gset@featureData@data
colnames( GPL )
view( GPL )

ids = GPL[ ,c( 1, 8 ) ]
ids = ids[ ids[ , 2 ] != '---' , ]

## 一个探针对应一个基因
library( "stringr" )
a= str_split( ids[,2], ' // ' , simplify = T )[ , 2 ]
tmp = cbind( ids[,1], a )
## 一个探针对应多个基因
a<-strsplit(as.character(ids[,2]), " /// ")
tmp <- mapply( cbind, ids[,1], a ) 
ID2gene <- as.data.frame( tmp )
colnames( ID2gene ) = c( "id", "gene" )
load( "./Relationship_protein_coding_gene.Rdata" )
ID2gene$type = gene2type[ match( ID2gene[ , 2 ], gene2type[ , 1 ] ), 2 ]
ID2gene = na.omit( ID2gene )

dim(ID2gene)
save(ID2gene, file = 'ID2gene.Rdata')


load('./ID2gene.Rdata')
load('./exprSet_by_group.Rdata')
## 去除没有注释的数据集
{
  exprSet = exprSet[ rownames(exprSet) %in% ID2gene[ , 1 ], ]
  ID2gene = ID2gene[ match(rownames(exprSet), ID2gene[ , 1 ] ), ]
}

dim( exprSet )
dim( ID2gene )
tail( sort( table( ID2gene[ , 2 ] ) ), n = 12L )

## 相同基因的表达数据取最大值
{
  MAX = by( exprSet, ID2gene[ , 2 ], 
  	        function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
  MAX = as.character(MAX)
  exprSet = exprSet[ rownames(exprSet) %in% MAX , ]
  rownames( exprSet ) = ID2gene[ match( rownames( exprSet ), ID2gene[ , 1 ] ), 2 ]
}
exprSet = log(exprSet)
dim(exprSet)
exprSet[1:5,1:5]

save(exprSet, group_list, file = 'final_exprSet.Rdata')


## plot
load('./final_exprSet.Rdata')
library( "ggfortify" )
## 聚类
{
  colnames( exprSet ) = paste( group_list, 1:ncol( exprSet ), sep = '_' )
  nodePar <- list( lab.cex = 0.3, pch = c( NA, 19 ), cex = 0.3, col = "red" )
  hc = hclust( dist( t( exprSet ) ) )
  png('hclust.png', res = 250, height = 1800)
  plot( as.dendrogram( hc ), nodePar = nodePar, horiz = TRUE )
  dev.off()
}
## PCA
data = as.data.frame( t( exprSet ) )
data$group = group_list
png( 'pca_plot.png', res=80 )
autoplot( prcomp( data[ , 1:( ncol( data ) - 1 ) ] ), data = data, colour = 'group', 
		  label =T, frame = T) + theme_bw()
dev.off()

