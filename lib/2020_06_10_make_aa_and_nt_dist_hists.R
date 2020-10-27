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

GetBinCounts = function( distVals, binVals )
{
  binCounts = sapply( 2:length(binVals),
    function(j) 
      sum( distVals[ distVals >= binVals[ j - 1 ] ] < binVals[ j ] )
    )
  return( binCounts )
}

GetAlgnLen = function( logFile )
{
  log  = readLines( logFile )
  idx  = grep( "alignment length is", log )
  strs = strsplit( log[idx], ' ' )[[1]]
  return( as.numeric( strs[length(strs) ] ) )
}

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

logFiles = system( "ls 2020_06_02_run_cognac_for_genus_*.log", intern = TRUE )
logFiles = logFiles[ -c(8) ]


MakeTransParent = function( color, percVal )
{
  rgbVals = col2rgb( color )
  
  ## Make new color using input color as base and alpha set by transparency
  tCol = rgb(
    rgbVals[1],
    rgbVals[2],
    rgbVals[3],
    max = 255,
    alpha = ( 100 - percVal ) * 255 / 100,
    )
  return( tCol )
}

lightCols = c( 
  ChangeHexTone( plotCols[1], 2),
  plotCols[2],
  ChangeHexTone( plotCols[3], 5),
  plotCols[4],
  ChangeHexTone( plotCols[5], 10),
  ChangeHexTone( plotCols[6], 2),
  ChangeHexTone( plotCols[7], 5),
  ChangeHexTone( plotCols[8], 2)
  )
darkCols  = c( 
  plotCols[1], 
  ChangeHexTone( plotCols[2], 0.15), 
  plotCols[3], 
  ChangeHexTone( plotCols[4], 0.65),
  plotCols[5], 
  plotCols[6], 
  plotCols[7], 
  plotCols[8]
  )
    

# lightCols = sapply( lightCols, function(x) MakeTransParent( x, 45 ) )
darkCols  = sapply( darkCols, function(x) MakeTransParent( x, 45 ) )
# ------------------------------------------------------------------------------

PlotDists = function( aaAlgnDists, ntAlgnDists, filename )
{
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
  
  #   binSize = 5000
  # xLabBin = 50000
  binSize = 0.001
  xLabBin = 0.01
  
  xLimMax = max(sapply( seq(nMats),
    function(i) ntDistValList[[ i ]][ length( ntDistValList[[ i ]] ) ]
    ))
  xLimMax =  ceiling(  xLimMax / xLabBin ) * xLabBin 
  
  # Create a vector with the size of the bins
  binVals = seq( 0, xLimMax, binSize )
  
  # Make a list with the bin values
  binCountList =  sapply( seq(nMats), function(i)
  {
    binCounts = rbind(
      GetBinCounts( aaDistValList[[i]], binVals ),
      GetBinCounts( ntDistValList[[i]], binVals )
      )
    return( list( binCounts ) )
  })
  
  
  # Get he x-axis labels 
  xLabels   = seq( 0, ceiling( xLimMax / xLabBin ) * xLabBin, by = xLabBin )
  filename = "../figures/2020_06_10_pairwise_aa_and_nt_algn_hist_normalized.png"
  png( 
    filename = filename,
    height = 750,
    width = 750
    )
  # Set up the plotting parameters
  par( mfrow = c( nMats / 2, 2 ) )
  for ( i in seq( nMats ) )
  {
    yRange   = GetAxisRange( max( binCountList[[i]] ) )
    xAxisPos = yRange[2] * -0.05
    
    aaCol = lightCols[i]
    ntCol = darkCols[i]
    
    invisible( sapply( 1:nrow( binCountList[[i]] ), function(j)
    {
      if ( j == 1 )
      {
        histCol = aaCol
      } else {
        histCol = ntCol
      }
      
      # Create a vector of colors where white is 0
      histCols = sapply( binCountList[[i]][ j, ], function(x)
      {
        if ( x ) return( histCol )
        return( "white" )
      })
      
      barplot(
        binCountList[[i]][ j, ], 
        main   = names( aaAlgnDists )[i], 
        cex.main = 1.5,
        border = NA,
        col    = histCols, 
        ylim   = yRange, 
        las    = 2,
        add    = ifelse( j == 1, FALSE, TRUE )
        )
    }))
    
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
    
    legend(
      "topright",
      legend = c( "Amino Acid", "Nucleotide" ),
      col = c( aaCol, ntCol ),
      pch = 15,
      cex = 1.5,
      bty = "n"
      )
  }
  dev.off()
}

# ------------------------------------------------------------------------------

PlotDists( 
  aaAlgnDists, 
  ntAlgnDists, 
  "../figures/2020_06_10_pairwise_aa_and_nt_algn_hist.png" 
  )

# ------------------------------------------------------------------------------

algnLens = sapply( logFiles, GetAlgnLen )

for ( i in 1:length(aaAlgnDists) ) 
  aaAlgnDists[[i]] = aaAlgnDists[[i]] / algnLens[i]

for ( i in 1:length(ntAlgnDists) ) 
  ntAlgnDists[[i]] = ntAlgnDists[[i]] / ( algnLens[i] * 3 )

PlotDists( aaAlgnDists, ntAlgnDists )

PlotDists( 
  aaAlgnDists, 
  ntAlgnDists, 
  "../figures/2020_06_10_pairwise_aa_and_nt_algn_hist_normalized.png" 
  )
sapply( 1:length(aaAlgnDists), function(i) max(aaAlgnDists[[i]]) )

sapply( ntAlgnDists, summary )

# ---- Save the data -----------------------------------------------------------

save( file = paste0( "../data/2020_06_09_.Rdata"), list = ls() )

# ------------------------------------------------------------------------------
