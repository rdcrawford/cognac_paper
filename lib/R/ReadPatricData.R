# ------------------------------------------------------------------------------
# ReadPatricData
# 2020/06/01
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Read in PATRIC data as a data.frame
# ------------------------------------------------------------------------------

ReadPatricData = function( tsvPath )
{
  read.table(
    tsvPath,
    sep              = '\t',
    header           = TRUE,
    colClasses       = c( "Genome.ID" = "character"),
    stringsAsFactors = FALSE
    )
}

# ------------------------------------------------------------------------------