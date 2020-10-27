# ------------------------------------------------------------------------------
# Run cognac for genus
# 2020/06/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This script runs cognac for the specific genus input in the command line 
# argument. Cognac is run removing any genomes missing greater than 1% of 
# genes from the analysis and alloing core genes to be present at a frequency
# of 5 per thousand genomes, and requiring that 1000 genes at minimum are 
# included in the alignment. The alignmet is then reverse translated and
# used to create a neighbor joining tree. The enviroment containing the 
# alignment data is then saved in an "Rdata" file.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )

# ---- Constant declarations --------------------------------------------------

# Column names in the meta-data 
GENOME_ID   = 1
SPECIES     = 2
GENUS       = 3
SOURCE      = 4
FASTA       = 5
GFF         = 6

# ---- Variable initializations ------------------------------------------------

# Get the genus for this run -- the only command line argument
runGenus = commandArgs()[6]

# Read in the meta-data
isolateMetaData = read.table( 
  "../data/2020_06_02_isolate_meta_data.tsv",
  sep = '\t',
  header = TRUE,
  stringsAsFactors = FALSE  
  )

# Subset to only include the genomes for this run
isGenus = isolateMetaData$Genus == runGenus
isolateMetaData = isolateMetaData[ isGenus, ]

# Create the output directory
analysisDir = "../analysis/2020_06_02_genus_cognac_aligments/"
outDir = paste0( analysisDir, runGenus, '/' )
if ( !file.exists(outDir) ) system( paste( "mkdir", outDir ) )

# Create the run ID to prepend to outpt files
runId = paste0( runGenus, "_cognac_algn" )

# ---- Run cognac --------------------------------------------------------------

algnEnv = cognac(
  fastaFiles    = isolateMetaData[ , FASTA ],
  featureFiles  = isolateMetaData[ , GFF ],
  genomeIds     = isolateMetaData[ , GENOME_ID ],
  runId         = runId,
  minGeneNum    = 1000,
  maxMissGenes  = 0.01,
  copyNumTresh  = 0.005,
  njTree        = TRUE,
  revTranslate  = TRUE,
  keepTempFiles = TRUE,
  outDir        = outDir
  )

# ---- Save the data -----------------------------------------------------------

save( file = paste0( outDir, runId, ".Rdata" ), list = ls() )

# ------------------------------------------------------------------------------

