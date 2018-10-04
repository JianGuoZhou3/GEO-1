rm( list = ls() )

## 建立基因和基因类型的对应关系

## 提取GTF文件中基因和基因类型的对应关系的shell脚本
## awk '{if(!NF || /^#/){next}}1' /public/reference/gtf/gencode/gencode.v25lift37.annotation.gtf | cut -f9 | sed 's/\"//g'| sed 's/;//g' | 
## awk '{ print $4"\t"$8 }' |awk '{if(/^E/){next}}1'|  awk '{ print $2"\t"$1 }' | sort -k 1 | uniq > gencode.v25lift37.annotation.gtf.gene2type

{
  gene2type = read.table( 'gencode.v25lift37.annotation.gtf.gene2type' )
  colnames( gene2type ) = c( "gene", "type" )
}

dim( gene2type )
sort( table( gene2type$type ) )
save( gene2type, file = 'Relationship_all_gene.Rdata' )
  
## 挑选基因类型为“protein_coding”的对应关系
gene2type = gene2type[ gene2type[,2] == 'protein_coding', ]
length( unique( gene2type$gene ) )
save( gene2type, file = 'Relationship_protein_coding_gene.Rdata' )

## 剔除基因类型为“protein_coding”的对应关系
gene2type = gene2type[ gene2type[,2] != 'protein_coding', ]
length( unique( gene2type$gene ) )
save( gene2type, file = 'Relationship_others_gene.Rdata' )

