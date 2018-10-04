rm( list = ls() )

load( "./nrDEG.Rdata" )

## 注释
library( "clusterProfiler" )
library( "org.Hs.eg.db" )
df <- bitr( rownames( nrDEG ), fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
head( df )
{
  nrDEG$SYMBOL = rownames( nrDEG )
  nrDEG = merge( nrDEG, df, by='SYMBOL' )
}
head( nrDEG )

{
  gene_up = nrDEG[ nrDEG$change == 'UP', 'ENTREZID' ] 
  gene_down = nrDEG[ nrDEG$change == 'DOWN', 'ENTREZID' ]
  gene_diff = c( gene_up, gene_down )
  gene_all = as.character(nrDEG[ ,'ENTREZID'] )
}

{
  geneList = nrDEG$logFC
  names( geneList ) = nrDEG$ENTREZID
  geneList = sort( geneList, decreasing = T )
}

{
  ## KEGG pathway analysis
  kk.up <- enrichKEGG(   gene          =  gene_up    ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.8        ,
                         qvalueCutoff  =  0.8        )
  kk.down <- enrichKEGG( gene          =  gene_down  ,
                         organism      =  'hsa'      ,
                         universe      =  gene_all   ,
                         pvalueCutoff  =  0.05       )
}

head( kk.up )[ ,1:6 ]
head( kk.down )[ ,1:6 ]

library( "ggplot2" )
{
  kegg_down_dt <- as.data.frame( kk.down )
  kegg_up_dt <- as.data.frame( kk.up )
  down_kegg <- kegg_down_dt[ kegg_down_dt$pvalue < 0.05, ]
  down_kegg$group = -1
  up_kegg <- kegg_up_dt[ kegg_up_dt$pvalue < 0.05, ]
  up_kegg$group = 1
  
  dat = rbind( up_kegg, down_kegg )
  dat$pvalue = -log10( dat$pvalue )
  dat$pvalue = dat$pvalue * dat$group
  
  dat = dat[ order( dat$pvalue, decreasing = F ), ]
  
  g_kegg <- ggplot( dat, 
                    aes(x = reorder( Description, order( pvalue, decreasing=F ) ), y = pvalue, fill = group)) + 
                    geom_bar( stat = "identity" ) + 
                    scale_fill_gradient( low = "blue", high = "red", guide = FALSE ) + 
                    scale_x_discrete( name = "Pathway names" ) +
                    scale_y_continuous( name = "log10P-value" ) +
                    coord_flip() + theme_bw() + theme( plot.title = element_text( hjust = 0.5 ) ) +
                    ggtitle( "Pathway Enrichment" ) 
  print( g_kegg )
  ggsave( g_kegg, filename = 'kegg_up_down.png' )
}

### GO database analysis 
g_list = list( gene_up = gene_up, gene_down = gene_down, gene_diff = gene_diff)

go_enrich_results <- lapply( g_list, function( gene ) {
  lapply( c( 'BP', 'MF', 'CC' ) , function( ont ) {
    cat( paste( 'Now process', ont ) )
    ego <- enrichGO( gene          =  gene,
                     universe      =  gene_all,
                     OrgDb         =  org.Hs.eg.db,
                     ont           =  ont ,
                     pAdjustMethod =  "BH",
                     pvalueCutoff  =  0.99,
                     qvalueCutoff  =  0.99,
                     readable      =  TRUE)
    print( head( ego ) )
    return( ego )
  })
})
save( go_enrich_results, file = 'go_enrich_results.Rdata' )

n1 = c( 'gene_up', 'gene_down', 'gene_diff' )
n2 = c( 'BP', 'MF', 'CC' ) 
for ( i in 1:3 ){
  for ( j in 1:3 ){
    fn = paste0( 'dotplot_', n1[i], '_', n2[j], '.png' )
    cat( paste0( fn, '\n' ) )
    png( fn, res = 150, width = 1080 )
    print( dotplot( go_enrich_results[[i]][[j]] ) )
    dev.off()
  }
}

