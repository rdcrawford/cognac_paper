# ------------------------------------------------------------------------------
# Sort Matrix
# 2020/04/13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function takes a matrix and another matrix to use as a reference. The 
# rows and columns of the first matrix are sorted by the row names and 
# column names of the second matrix. The sorted matrix is returned.
# ------------------------------------------------------------------------------

SortMatrix = function( mat, refMat )
{
  # Make sure there are actually row names and column names to sort by
  if ( 
    is.null(row.names(mat)) || is.null(row.names(refMat)) ||
    is.null(colnames(mat)) || is.null(colnames(refMat))
    ) 
    stop( 
      "Both the reference matrix and matrix to be sorted must have row and",
      " column names" 
      )
  
  # Find the order of the rows and columns 
  refRowOrder = sapply( row.names( refMat ),
    function(x) which( row.names( mat ) == x )
    )
  refColOrder = sapply( colnames( refMat ),
    function(x) which( colnames( mat ) == x )
    )
  mat = mat[ refRowOrder, refColOrder ]
  if ( !identical( dimnames(mat), dimnames(refMat) ) )
    stop( "You did not sort this properly" )
  return( mat )
}

# ------------------------------------------------------------------------------