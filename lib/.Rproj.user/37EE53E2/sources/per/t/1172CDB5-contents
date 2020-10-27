
load( "../data/rdata/plotCols.Rdata" )
load( "../data/2020_06_17_get_resource_usage.Rdata" )

logFiles = system("ls 2020_06_02_run_cognac_for_genus_*.log", intern = TRUE )
rId = "2020_06_02_run_cognac_for_genus_"
genusNames = sapply( logFiles, 
  function(x) strsplit( gsub(rId, '', x), '-' )[[1]][1], USE.NAMES = FALSE
  )

# ------------------------------------------------------------------------------

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

# This function 
GetElementAtIdx = function( inLine, idx )
{
  line = strsplit( inLine, ' ' )[[1]]
  line = line[ line != '' ]
  line = line[ 2:length(line) ]
  return( as.numeric( line[ idx ] ) )
}

GetValFromLog = function( str, idx, log )
{
  lineIdxs =  grep( str, log )
  lineVals = sapply( lineIdxs, function(i) GetElementAtIdx( log[ i ], idx ) )
  return( list(lineVals) )
}

GetRunStats = function( log )
{
  lineStrs = c( 
    "genomes were input", 
    "Running on:", 
    "Parsed the data for", 
    "to be included in the alignment", 
    "The total alignment length",
    "Finished in"
    )
  
  lineIdxs = c( 1, 3, 5, 1, 6, 3 )
  
  runStatList = sapply( 1:length(lineStrs),
    function(i) GetValFromLog( lineStrs[i], lineIdxs[i], log )
    )
  
  lineNames = c(
    "Input genomes", 
    "Number of cores", 
    "Total number of genes",
    "Number of core genes", 
    "Allignment length",
    "Step times"
    )
  
  names( runStatList ) = lineNames

  return( list( runStatList ) )
}

MakeBarPlot = function( runStats, lIdx, yLab )
{
  plotVals = sapply( 1:length(runStats), 
    function(i) runStats[[i]][[ lIdx ]] 
    )
  names(plotVals) = names( runStats )
  par( mar = c(7.1, 8.5, 4.1, 1.1) )
  barplot(
    plotVals,
    col    = plotCols,
    border = plotCols,
    main   = names( runStats[[1]] )[lIdx],
    las    = 2,
    ylim   = genomeToolBox::GetAxisRange( max(plotVals), 5 ),
    ylab = yLab
    )
}

# ------------------------------------------------------------------------------

plotCols[7] = ChangeHexTone( plotCols[7], 3)
runStats = sapply( logFiles, function(x) GetRunStats( readLines(x) ) )
sapply( runStats, length )
names(runStats) = genusNames

nGenomes = sapply( 1:length(runStats), 
  function(i) runStats[[i]][[ 1 ]] 
  )
# runStats = runStats[ order( nGenomes, decreasing = TRUE ) ]

timeMat = sapply( 1:length(runStats), 
  function(i) runStats[[i]][[6]]
  )


statMatrix = sapply( 1:length(runStats), function(i) 
{
  sapply( seq(5), function(j) runStats[[i]][j] )
})

totalTime = timeMat[ 10,  ]
statMatrix = rbind( statMatrix, totalTime, resources[ 2, ] )

# Set the names for the matrix
stepLabs = c(
  "Number of genomes", 
  "Number of cores", 
  "Total number of genes",
  "Number of core genes", 
  "Allignment length (aa residues)",
  "Run time (minutes)",
  "Memory Usage (GB)"
  )

row.names( statMatrix ) = stepLabs
colnames( statMatrix ) = genusNames

statMatrix = t( statMatrix )
statMatrix = statMatrix[ , -c(2) ]

write.table( statMatrix, file = "../data/2020_06_24_run_stats.tsv", sep = '\t')

# ------------------------------------------------------------------------------

CEX_LAB  = 2
CEX_AXIS = 2
CEX_PCH  = 3
CEX_LGD  = 2

par( mar = c( 7.1, 7.5, 4.1, 4) )

plot(
  x        = statMatrix[ , 1 ],
  y        = totalTime,
  col      = plotCols,
  pch      = 19,
  ylab     = "Total run time (minutes)\n",
  xlab     = "Number of genomes",
  cex      = CEX_PCH,
  las      = 1,
  cex.lab  = CEX_LAB,
  cex.axis = CEX_AXIS
  )
legend(
  "topleft",
  col = plotCols,
  legend = genusNames,
  pch = 15,
  cex = CEX_LGD
  )


timeMat[ 6, ] = sapply( 1:ncol( timeMat), function(j) sum( timeMat[ 6:9, j] ) )
timeMat = timeMat[ 1:6, ]

colnames(timeMat) = genusNames
set.seed( 144 )
bpCols = plotCols[ sample.int( length(plotCols) ) ]

steps = c(
  "Parsing input files",
  "Finding orthologs with cd-hit",
  "Filtering for single copy genes",
  "Selecting  genes to align",
  "Aligning and concatenating genes",
  "Creating output files"
  )
par( mar = c( 12.2, 7.5, 4.1, 0.5) )
barplot(
  timeMat[ , order(totalTime, decreasing = TRUE) ],
  col      = bpCols,
  border   = bpCols,
  las      = 2,
  ylim     = genomeToolBox::GetAxisRange( max( colSums(timeMat) ), 5 ),
  ylab     = "Run-time (minutes)\n",
  cex      = 2,
  cex.lab  = CEX_LAB,
  cex.axis = CEX_AXIS
  )

legend( 
  "topright",
  legend = steps,
  col = bpCols,
  pch = 15,
  cex = CEX_LGD
  )

# ------------------------------------------------------------------------------

timeFracMat = sapply( 1:ncol(timeMat),
  function(j) timeMat[ , j ] / sum( timeMat[ , j ] )
  )
timeFracMat = timeFracMat * 100
colnames(timeFracMat) = genusNames
par( mfrow = c(1,1))
par( mar = c( 12.2, 7.5, 4.1, 0.5) )
barplot(
  timeFracMat[ , order(totalTime, decreasing = TRUE) ],
  col    = bpCols,
  border = bpCols,
  las    = 2,
  ylim   = c( 0 , 100),
  ylab   = "Percent of Run-time",
  cex      = 2,
  cex.lab  = CEX_LAB,
  cex.axis = CEX_AXIS
  )

fnames( seqCountList )

x = which.max( seqCountList[[4]] )
seqCountList[[4]][x]

# ------------------------------------------------------------------------------

load("../data/2020_06_10_count_unique_alleles.Rdata")

# Normalize the number of unique alleles to the number of genomes 
seqCountList = sapply( 1:length(seqCountList), 
  function(i) seqCountList[[i]] / unlist( statMatrix[ i, 1 ] ) 
  )

summaryVals = sapply( seqCountList, summary )
colnames(summaryVals) = genusNames
summaryVals * 100
maxVals = sapply( seqCountList, max )

x = rbind( minVals, maxVals ) * 100
colnames( x ) = genusNames 

par(mfrow = c( 4, 2 ) )
for ( i in 1:length( seqCountList ) )
{
  uniqueFrac = seqCountList[[i]]
  barplot(
    uniqueFrac[ order(uniqueFrac, decreasing = TRUE) ],
    col    = plotCols[i],
    border = plotCols[i],
    main   = genusNames[i],
    ylim   = c( 0, 1 ),
    las    = 2,
    cex.main = 2
    )
}

# ------------------------------------------------------------------------------

medSeqCount = sapply( 1:length(seqCountList), 
  function(i) median( seqCountList[[i]] ) 
  )

nGenomes = unlist( statMatrix[ , 1 ] )

par(mfrow = c( 1, 1 ) )
plot(
  x        = medSeqCount,
  y        = nGenomes,
  col      = plotCols,
  pch      = 19,
  ylab     = "Number of genomes",
  xlab     = "Median fraction of unique allelels",
  xlim     = c( 0 , 0.25),
  cex      = 2,
  cex.lab  = 1.5,
  cex.axis = 1.5
  )
legend(
  "topright",
  col = plotCols,
  legend = genusNames,
  pch = 15,
  cex = 1.5
  )

normsapply( 1:length(seqCountList), function(i) medSeqCount[i] / nGenomes[i] )
# ------------------------------------------------------------------------------

