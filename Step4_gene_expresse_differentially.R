rm( list = ls() )

## 差异分析
load( "./final_exprSet.Rdata" )
library( "limma" )
{
  design <- model.matrix( ~0 + factor( group_list ) )
  colnames( design ) = levels( factor( group_list ) )
  rownames( design ) = colnames( exprSet )
}
design

contrast.matrix <- makeContrasts( "gbm-normal", levels = design )
contrast.matrix

load( "./ID2gene.Rdata" )
{
  fit <- lmFit( exprSet, design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  nrDEG = topTable( fit2, coef = 1, n = Inf )
  nrDEG$geneType = g2t[ match( rownames( nrDEG ), g2t$gene ), 2] 
  write.table( nrDEG, file = "nrDEG.out")
}
head(nrDEG)

## heatmap
library( "pheatmap" )
{
  choose_gene = head( rownames( nrDEG ), 50 )
  choose_matrix = exprSet[ choose_gene, ]
  choose_matrix = t( scale( t( exprSet ) ) )

  choose_matrix[choose_matrix > 2] = 2
  choose_matrix[choose_matrix < -2] = -2

  annotation_col = data.frame( CellType = factor( group_list ) )
  rownames( annotation_col ) = colnames( exprSet )
  pheatmap( fontsize = 2, choose_matrix, annotation_col = annotation_col, show_rownames = F, 
            annotation_legend = F, filename = "heatmap.png")
}

## volcano plot
library( "ggplot2" )
logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
logFC_cutoff
logFC_cutoff = 1.2

{
  nrDEG$change = as.factor( ifelse( nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) > logFC_cutoff,
                                  ifelse( nrDEG$logFC > logFC_cutoff , 'UP', 'DOWN' ), 'NOT' ) )
  
  save( nrDEG, file = "nrDEG.Rdata" )

  this_tile <- paste0( 'Cutoff for logFC is ', round( logFC_cutoff, 3 ),
                       '\nThe number of up gene is ', nrow(nrDEG[ nrDEG$change =='UP', ] ),
                       '\nThe number of down gene is ', nrow(nrDEG[ nrDEG$change =='DOWN', ] ) )
  
  volcano = ggplot(data = nrDEG, aes( x = logFC, y = -log10(P.Value), color = change)) +
				   geom_point( alpha = 0.4, size = 1.75) +
				   theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
				   xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
				   ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
				   scale_colour_manual( values = c('blue','black','red') )
  print( volcano )
  ggsave( volcano, filename = 'volcano.png' )
  dev.off()
}
