#  -----------------------------------------------------------------------------
# 
# 2020/05/20
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#  
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( genomeToolBox )

# ---- Function definitions ----------------------------------------------------

ChangeHexTone = function( hexSymbols, percIntensity )
{
  rgbVals = col2rgb( hexSymbols )
  rgbVals = round( rgbVals * percIntensity, 0 )
  
  if ( TRUE %in% ( rgbVals > 255 ) )
  {
    for ( i in seq( 3 ) ) 
    {
      isTooHigh = rgbVals[ i , ] > 255
      if ( TRUE %in% isTooHigh )
      {
        for ( j in which( isTooHigh ) ) rgbVals[ i, j ] = 255
      }
    }
  }
  
  hexVals = sapply( seq( ncol(rgbVals) ), function( j)
    rgb( 
      rgbVals[ 1, j ], rgbVals[ 2, j ], rgbVals[ 3, j ], maxColorValue = 255
      )
    )
  return( hexVals )
}

lightCols = plotCols
darkCols = plotCols

PlotAxis = function( )
{

}

GetDistMat = function( rData )
{
  load( rData )
  return( algnEnv$distMat )
}

plotCols = c( plotCols, sapply( plotCols, function(x) ChangeHexTone( x, 0.33 )))

# ---- Variable initializations -------------------------------------------------

load("../data/2020_06_09_aa_algn_dists.Rdata")

load( "../data/2020_06_09_nt_algn_dists.Rdata")


load("../data/rdata/plotCols.Rdata")



rDataFiles = system(
  "ls ../analysis/2020_06_02_genus_cognac_aligments/*/*_cognac_algn.Rdata",
  intern = TRUE
  )
names( rDataFiles ) = sapply( rDataFiles, 
  function(x) strsplit( x, '/' )[[1]][4], USE.NAMES = FALSE
  )
isOut = names( rDataFiles ) %in% c( "Escherichia", "Klebsiella" )

rDataFiles = rDataFiles[ !isOut ]
ntAlgnDists = sapply(rDataFiles, function(x) list( GetDistMat( x ) ) )

aaAlgnDists = aaAlgnDists[ 
  !names(aaAlgnDists) %in% c( "Escherichia", "Klebsiella" )
  ]


# ------------------------------------------------------------------------------

# Assign default values to missing arguments
nMats = length( aaAlgnDists )


# Create a list where each elment contains a vector of the pairwise
# distances between genomes
aaDistValList = sapply( seq( nMats ), function(i)
{
  distVals = ConvertDistMatToVec( aaAlgnDists[[i]] )
  return( list( distVals[ order(distVals) ] ) )
})

ntDistValList = sapply( seq( nMats ), function(i)
{
  distVals = ConvertDistMatToVec( ntAlgnDists[[i]] )
  return( list( distVals[ order(distVals) ] ) )
})


binSize = 1000
xLabBin = 5000

xLimMax = max(sapply( seq(nMats),
  function(i) aaDistValList[[ i ]][ length( aaDistValList[[ i ]] ) ]
  ))
xLimMax = max( ceiling( xLimMax / xLabBin ) * xLabBin )

# Create a vector with the size of the bins
binVals = seq( 0, xLimMax, binSize )

# Make a list with the bin values
binCountList =  sapply( seq(nMats), function(i)
{
  distVals  = aaDistValList[[i]]
  binCounts = sapply( 2:length(binVals),
    function(j) 
      sum( distVals[ distVals >= binVals[ j - 1 ] ] < binVals[ j ] )
    )
  return( list( binCounts ) )
})


# Get he x-axis labels 
xLabels   = seq( 0, ceiling( xLimMax / xLabBin ) * xLabBin, by = xLabBin )

png( 
  filename = "../figures/2020_06_10_pairwise_aa_algn_hist.png",
  height = 1000,
  width = 750
  )
# Set up the plotting parameters
par( mfrow = c( nMats / 2, 2 ) )
for ( i in seq( nMats ) )
{
  yRange   = GetAxisRange( max( binCountList[[i]] ) )
  xAxisPos = yRange[2] * -0.05
  
  # Create a vector of colors where white is 0
  histCols = sapply( binCountList[[i]], function(x)
  {
    if ( x ) return( plotCols[i] )
    return( "white" )
  })
  
  # Make the barplot
  barplot(
    binCountList[[i]],
    main   = names(aaAlgnDists)[i],
    col    = histCols,
    border = histCols,
    ylab   = "",
    ylim   = yRange,
    las    = 2,
    )
  
  axisStart = par()$ylbias #* par()$plt[1]
  axisEnd   = ( par()$usr[2] * ( par()$plt[2] + 0.025 ) ) + axisStart
  atVals    = seq( axisStart, axisEnd, length.out = length(xLabels) )
  axis(
    side   = 1,
    at     = atVals,
    labels = as.character( xLabels ),
    las    = 3,
    pos    = xAxisPos
    )
}
dev.off()

# ------------------------------------------------------------------------------



# ---- Save the data -----------------------------------------------------------

save( file = paste0( "../data/2020_06_09_.Rdata"), list = ls() )

# ------------------------------------------------------------------------------
