#  -----------------------------------------------------------------------------
#
# 2020/05/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#  
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( genomeToolBox )

# ---- Variable initializations ------------------------------------------------

patricDataPaths = system( "ls ../data/2020_*_PATRIC_genome.tsv", intern = TRUE )

columnNames = c( "GenomeId", "Species", "Genus", "Source", "Fasta", "Gff" )

# ---- Constant declarations --------------------------------------------------

GENOME_ID   = 1
GENOME_NAME = 2
SPECIES     = 2
GENUS       = 3
SOURCE      = 4
FASTA       = 5
GFF         = 6

# ---- Function definitions ----------------------------------------------------

# Function to read in the patric tsvs
source( "R/ReadPatricData.R" )

# This function returns a length 2 character vector: 1) species name and 2)
# the genus name
GetSpecies = function( genomeName )
{
  species = strsplit( genomeName, ' ')[[1]]
  if ( length( species) > 1 ) species = species[ c( 1, 2 ) ]
  species = gsub("[[:punct:]]", "", species) 
  return( c( paste( species, collapse = ' '), species[1] ) )
}

# ------------------------------------------------------------------------------

# Merge the individual data frames
patricData = ReadPatricData( patricDataPaths[1] )
for ( i in 2:length(patricDataPaths) )
  patricData = rbind( patricData, ReadPatricData( patricDataPaths[i] ) )
head( patricData )

# Remove any genomes which are duplicated
isDuplicated = duplicated( patricData[ , GENOME_ID ] )
sum(isDuplicated)
patricData = patricData[ !isDuplicated, ]
nrow( patricData )

# Get the paths to the genomic data
patricDataPath = "../../patric_genomes/"
patricFas = sapply( patricData[ , GENOME_ID ], 
  function(x) paste0( patricDataPath, x, '/', x, ".fasta" )
  )
patricGffs = sapply( patricData[ , GENOME_ID ], 
  function(x) paste0( patricDataPath, x, '/', x, ".gff" )
  )

# Create the meta data for the patric genomes
patricMetaData = cbind.data.frame(
  patricData[ , GENOME_ID ], 
  t( sapply( patricData[ , GENOME_NAME ], GetSpecies, USE.NAMES = FALSE ) ),
  rep( "PATRIC", nrow( patricData ) ),
  patricFas,
  patricGffs,
  stringsAsFactors = FALSE
  )
colnames( patricMetaData ) = columnNames 


# Remove any phage genus
isPhage        = grepl( "phage",  patricMetaData[ , SPECIES ] )
patricMetaData = patricMetaData[ !isPhage, ]
patricData     = patricData[ !isPhage, ]

# ------------------------------------------------------------------------------

# Read in the meta-data on the LTACH genomes
ltachMetaData = read.table(
  "../../hgt_project/data/isolateMetaData.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
  )

ltachPath = "../../hgt_project/data/ltach_genomes/"
ltachFas  = sapply( ltachMetaData[ , GENOME_ID ],
  function(x) paste0( ltachPath, x, ".fasta" )
  )
ltachGffs  = sapply( ltachMetaData[ , GENOME_ID ],
  function(x) paste0( ltachPath, x, ".gff" )
  )
ltachMetaData = cbind.data.frame(
  ltachMetaData[ , 1:2 ],
  sapply( ltachMetaData[ , SPECIES ], function(x) strsplit(x, ' ')[[1]][1] ),
  rep( "LTACH", nrow( ltachMetaData ) ),
  ltachFas,
  ltachGffs,
  stringsAsFactors = FALSE
  )
colnames( ltachMetaData ) = columnNames
head( ltachMetaData )

# ------------------------------------------------------------------------------

isRunGenus     = patricMetaData[ , GENUS ] %in% ltachMetaData[ , GENUS ]
patricMetaData = patricMetaData[ isRunGenus, ]
patricData     = patricData[ isRunGenus, ]

isolateMetaData = rbind.data.frame(
  patricMetaData,
  ltachMetaData,
  stringsAsFactors = FALSE
  )

specisCounts = table( isolateMetaData[ , SPECIES ] )
specisCounts = specisCounts[ order(specisCounts, decreasing = TRUE ) ]
specisCounts
genusCounts = table( isolateMetaData[ , GENUS ] )
genusCounts = genusCounts[ order(genusCounts, decreasing = TRUE ) ]
genusCounts
write.table( 
  isolateMetaData,
  file = "../data/2020_06_02_isolate_meta_data.tsv",
  sep = '\t'
  )
write.table(
  patricData,
  file = "../data/2020_06_02_PATRIC_genome_meta_data_merged.tsv"
  )

# ------------------------------------------------------------------------------
