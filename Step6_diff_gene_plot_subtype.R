rm( list = ls() )

load( "./gset.Rdata" )
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
  n_expr = exprSet[ , grep( 'non-tumor',         group_list )]
  g_expr = exprSet[ , grep( 'glioblastoma',      group_list )]
  a_expr = exprSet[ , grep( 'astrocytoma',       group_list )]
  o_expr = exprSet[ , grep( 'oligodendroglioma', group_list )]
  exprSet=cbind(n_expr, o_expr, a_expr, g_expr)
}


## 样本分组
{
  group_list=c(rep('NC',ncol(n_expr)),
               rep('od',ncol(o_expr)),
               rep('a',ncol(a_expr)),
               rep('GBM',ncol(g_expr)))
}

load('./ID2gene.Rdata')
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

save(exprSet, file = 'all_exprSet.Rdata')


library( "ggstatsplot" )
load( './final_exprSet.Rdata' )
load( "./nrDEG.Rdata" )
nrDEG_Z = nrDEG[ order( nrDEG$logFC ), ]
nrDEG_F = nrDEG[ order( -nrDEG$logFC ), ]
special_gene = c( rownames( nrDEG_Z )[1:5], rownames( nrDEG_F )[1:5] )
## 挑出基因作图
for( gene in special_gene ){
  filename <- paste( gene, '.png', sep = '' )
  TMP = exprSet[ rownames( exprSet ) == gene, ]
  data = as.data.frame(TMP)
  data$group = group_list
  p <- ggbetweenstats(data = data, x = group,  y = TMP )
  ggsave( p, filename = filename)
}



