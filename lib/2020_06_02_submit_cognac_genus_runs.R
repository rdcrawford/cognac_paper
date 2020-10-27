# ------------------------------------------------------------------------------
# Submit genus cognac jobs
# 2020/06/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This script writes and submits the slurm scripts for the genus specific 
# cognac runs 
# ------------------------------------------------------------------------------

# ---- Constant declarations --------------------------------------------------

GENUS = 3

# ---- Variable initializations ------------------------------------------------

# Read in the meta-data
isolateMetaData = read.table( 
  "../data/2020_06_02_isolate_meta_data.tsv",
  sep = '\t',
  header = TRUE,
  stringsAsFactors = FALSE  
  )

slurmTemplate = readLines( "slurm_template.sh" )

# ------------------------------------------------------------------------------

CreateJob = function( genus, slurmTemplate )
{
  # Find the indicies of relevant parameters
  jIdx  = grep( "#SBATCH --job-name=", slurmTemplate )
  cIdx  = grep( "#SBATCH --cpus-per-task=", slurmTemplate )
  mIdx  = grep( "#SBATCH --mem-per-cpu=", slurmTemplate )
  pIdx  = grep( "#SBATCH --partition=standard", slurmTemplate )
  runId = paste0( "2020_06_02_run_cognac_for_genus", "_", genus )
  
  # Set the run id
  slurmTemplate[ jIdx ] = paste0( "#SBATCH --job-name=", runId )
  
  # Get the number of genomes to be included
  numGenomes = sum( isolateMetaData[ , GENUS ] == genus )
  
  # IF the number of genomes is greater than 10,000, submit a job on a 
  # large memor node
  if ( numGenomes  > 10000 )
  {
   
   slurmTemplate[ cIdx ] = "#SBATCH --cpus-per-task=12"
   slurmTemplate[ mIdx ] = "#SBATCH --mem-per-cpu=32g"
   slurmTemplate[ pIdx ] = "#SBATCH --partition=largemem"
   
  # Submit a job on 12 cores with 8 gigs if this is a sm
  } else if ( numGenomes  > 1000 ) {
  
   slurmTemplate[ cIdx ] = "#SBATCH --cpus-per-task=12"
   slurmTemplate[ mIdx ] = "#SBATCH --mem-per-cpu=7900m"
    
  # Submit a job on 4 cores with 4 gigs if this is a s
  } else if ( numGenomes  > 100 ) {
  
   slurmTemplate[ cIdx ] = "#SBATCH --cpus-per-task=12"
   slurmTemplate[ mIdx ] = "#SBATCH --mem-per-cpu=1000m"
    
  # Submit a job on one core with 8 gigs if this is a small job
  } else {
    
   slurmTemplate[ cIdx ] = "#SBATCH --cpus-per-task=12"
   slurmTemplate[ mIdx ] = "#SBATCH --mem-per-cpu=100m"
  }
  
  slurmTemplate[ length(slurmTemplate) + 1 ] = 
    paste( 
      "\nml cd-hit\nml mafft\n\nRscript 2020_06_02_run_cognac_for_genus.R", 
      genus
      )
  
  sbatPath = paste0( runId, ".sh" )
  sink( sbatPath )
  for ( x in slurmTemplate ) cat( x, '\n' )
  sink()
  
  runId = system( paste( "sbatch", sbatPath ), intern = TRUE )
  return( runId )
}

# ------------------------------------------------------------------------------

runIds = sapply( unique( isolateMetaData[ , GENUS ] ),
  function(x) CreateJob( x, slurmTemplate )
  )

# ---- Save the data -----------------------------------------------------------

save( file = "../data/2020_06_02_submit_cognac_genus_runs.Rdata", list = ls() )

# ------------------------------------------------------------------------------