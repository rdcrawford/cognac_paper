#  -----------------------------------------------------------------------------
# 
# 2020/06/01
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( genomeToolBox )

# ---- Variable initializations ------------------------------------------------

# Path to the patric geomes
patricDir = "../../patric_genomes/"

# Get the paths to the complete Kleb genomes
completePaths = system(
  "ls ../../complete_kleb_genomes/*/*.gff", intern = TRUE
  )

# Paths to the patric genomes for the enterobacter project
enteroBacterPaths = system(
  "ls ../../2020_04_13_CRE_EIP_core_gene_analysis/data/patric_genomes/*.gff",
  intern = TRUE
  )

# Paths to the patric genomes
patricPaths = system( paste0( "ls ", patricDir, "*/*.gff" ), intern = TRUE )

# Get the paths to the downloaded genomes
dwnldPaths = system( "ls ../data/genomes/*/*.gff", intern = TRUE )

# Total the number of geomes already downloaded
totals = sum(
  length( dwnldPaths ),
  length( patricPaths ),
  length( enteroBacterPaths ),
  length( completePaths )
  )

cat( 
  "Complete Paths: ",     length( completePaths ),     '\n',
  "Enterobacter Paths: ", length( enteroBacterPaths ), '\n',
  "Patric Paths: ",       length( patricPaths ),       '\n',
  "Downloaded Paths: ",   length( dwnldPaths ),        '\n',
  "Total: ",              totals,                      '\n',
  sep = ''
  )

genomeIds = sapply( 
  c( patricPaths, dwnldPaths, completePaths, enteroBacterPaths ), 
  ExtractGenomeNameFromPath, USE.NAMES = FALSE 
  )
patricDataPaths = system(
  "ls ../data/2020_05_27_*_PATRIC_genome.tsv", intern = TRUE
  )

# ---- Function definitions ----------------------------------------------------

source( "R/ReadPatricData.R" )

# ------------------------------------------------------------------------------

for ( tsvPath in patricDataPaths )
{
  patricData = ReadPatricData( tsvPath )
  
  for ( x in  patricData$Genome.ID )
  {
    if ( !x %in% genomeIds )
    {
      DownloadPatricGenome( x, "../data/genomes" )
    }
  }
}

# ------------------------------------------------------------------------------

toMove    = c( enteroBacterPaths, completePaths )
toMoveIds = sapply( toMove, ExtractGenomeNameFromPath, USE.NAMES = FALSE )

for ( i in 2:length( toMove ) )
{
  genomeDir = paste0( patricDir, toMoveIds[i], '/' )
  system( paste( "mkdir", genomeDir ) )
  system( paste("cp", toMove[i], genomeDir ) )
  system( paste("cp", gsub( "gff", "fasta", toMove[i] ), genomeDir ) )
}
  
# ------------------------------------------------------------------------------

patricData = ReadPatricData( patricDataPaths[1] )
for ( i in 2:length( patricDataPaths ) ) 
{
  patricData = rbind.data.frame( 
    patricData, ReadPatricData( patricDataPaths[i] ), stringsAsFactors = FALSE
    )
}
dim( patricData )  

# ------------------------------------------------------------------------------

metaData = read.table( 
  "../data/2020_06_02_isolate_meta_data.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
  )
dim( metaData )
isLtach = metaData$Source == "LTACH"
sum( isLtach )
tail( metaData[ isLtach, ], 90 )

ltachData = data.frame( matrix(
  NA, nrow = sum( isLtach ), ncol = ncol( patricData ) 
  ))
dim( ltachData )
colnames( ltachData ) = colnames( patricData ) 
ltachData[ , 1 ] = metaData[ isLtach, 1 ]
ltachData[ , 2 ] = paste( metaData[ isLtach, 1 ], metaData[ isLtach, 2 ] )
dim( ltachData )

isolateMetaData = rbind.data.frame( 
  patricData, ltachData, stringsAsFactors = FALSE
  )
write.table( 
  isolateMetaData, 
  file = "../data/Supplementary_table_1.tsv",
  row.names = FALSE,
  col.names = TRUE,
  sep = '\t'
  )
