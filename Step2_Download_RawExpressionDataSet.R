rm( list = ls() )

## 使用"GEOquery"包下载GEO数据
## 无特殊情况，仅修改“GSE_name”即可
library( "GEOquery" )
GSE_name = 'GSE90604'
options( 'download.file.method.GEOquery' = 'libcurl' )
{
  possibleError <- tryCatch(
    gset <- getGEO( GSE_name, getGPL = T ),
    error = { function(e) e },
  )
  if( inherits( possibleError, "error" ) ){

    ## 仅适用一个平台注释的数据
    stub = gsub( "\\d{1,3}$", "nnn", GSE_name, perl = TRUE )
    GSEfilename <- paste( GSE_name, 'series_matrix.txt.gz', sep = '_' )
    url <- paste( 'https://ftp.ncbi.nlm.nih.gov/geo/series/', 
                    stub, '/', GSE_name, '/matrix/', GSEfilename, sep = '' )
      
    ## 多个平台手动复制下载地址
    ## url <- ""
    ## GSEfilename = strsplit( url, split = "/" )[[1]][9]
      
    download.file( url, GSEfilename, method = "libcurl" )
    GSEfilePATH <- paste( './', GSEfilename, sep = '' )
    gset <- getGEO( filename = GSEfilePATH, getGPL = F )
  }
}
save( gset, file = 'gset.Rdata' )

## trying URL 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90604/matrix/GSE90604-GPL17692_series_matrix.txt.gz'
