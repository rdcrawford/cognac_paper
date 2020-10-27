#  -----------------------------------------------------------------------------
# 
# 2020/05/20
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#  
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )

# ---- Variable initializations ------------------------------- -----------------

aaAlgnPaths = paste0(
  "../analysis/2020_06_02_genus_cognac_aligments/*/",
  "*_cognac_algn_concatenated_gene_aa_alignment.fasta"
  )
aaAlgnPaths = system( paste( "ls", aaAlgnPaths), intern = TRUE )

ntAlgnPaths = paste0(
  "../analysis/2020_06_02_genus_cognac_aligments/*/",
  "*_cognac_algn_concatenated_gene_nt_alignment.fasta"
  )
ntAlgnPaths = system( paste( "ls", ntAlgnPaths), intern = TRUE )

# ------------------------------------------------------------------------------

ntAlgnDists = sapply( ntAlgnPaths, function(algn)
{
  list( cognac::CreateAlgnDistMat( algn, "raw", FALSE ) )
})
names( ntAlgnDists ) = 
  sapply( ntAlgnPaths, function(x) strsplit( x, '/' )[[ 1 ]][ 4 ] )

# ---- Save the data -----------------------------------------------------------

save( 
  file = paste0( "../data/2020_06_09_nt_algn_dists.Rdata"),
  list = "ntAlgnDists"
  )

# ------------------------------------------------------------------------------

# aaAlgnDists = sapply( aaAlgnPaths, function(algn)
# {
#   list( cognac::CreateAlgnDistMat( algn, "raw", FALSE ) )
# })
# names( aaAlgnDists ) = 
#   sapply( aaAlgnPaths, function(x) strsplit(x, '/')[[1]][4] )
# 
# # ---- Save the data -----------------------------------------------------------
# 
# save( 
#   file = paste0( "../data/2020_06_09_aa_algn_dists.Rdata"),
#   list = "aaAlgnDists"
#   )

# ------------------------------------------------------------------------------
