library(rgrass7)
library(rgdal)
library(raster)
library(cluster)
library(randomForest)


loc <- initGRASS("/home/epilib/Envs/env1/grass-7.1.svn", home=tempdir())
execGRASS('g.proj',epsg=as.integer(32633), location='utm_wgs84_33N')
execGRASS('g.mapset', mapset='PERMANENT', location='utm_wgs84_33N')

# run only once ...
## create a raster for use in grass
#execGRASS('r.in.gdal', input="bathy20_1.tif",
#          flags=c('overwrite','o','e'), output='elev')     # o= override the prjection parameters, and e=extend the location

## set the region to the extent of the file
#execGRASS('g.region', raster='elev', res='20', flags=c('a','p'))    # align cells with region extent

# NOTE:
# bathy20_1 is a large raster. 
# Better work on a sub-dataset to test the code first, 
# once happy with the results on a sub-dataset run the analysis on the full dataset
# this 2 lines will generate a sub-datset "minielev" with small extent and same resolution :

execGRASS('g.region', raster='elev', res='20', w='580000',e='625000',s='7770000',n='7840000', flags=c('a','p')) 
execGRASS('r.mapcalc', expression = 'minielev=elev', flags='overwrite') 

source('morphoclara.R')
s <- morphoclara(elevation='minielev', remove=TRUE)

# we can't operate on the entire set of cells,
# so we use one of the central concepts from statistics: sampling
# sample 10000 random points
s.r <- as.data.frame(sampleRegular(s, 10000))

# clara() function: need to remove NA from training set
s.r <- na.omit(s.r)

# from :
# list with silhouette width information for the best sample, see partition.object.
#
# https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/partition.object.html
# https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/clara.object.html

ncl <- numeric(8)
for (i in 2:9) ncl[i] <- clara(s.r,stand=TRUE,k=i)$ silinfo $ avg.width
plot(2:10,ncl,type="b")

s.clara <- clara(s.r, stand=TRUE, k=6)
s.r$cluster <- factor(s.clara$clustering)

rf <- randomForest(cluster ~ er + xslope + profc + crosc + minic + maxic
                   + longc, data=s.r, importance=TRUE, ntree=201)

# make predictions from rf model, along all cells of input stack
p <- predict(s, rf, type='response', progress='text')

# variable importance: all contribute about the same... good
#importance(rf)
#varImpPlot(rf)

# customized plot (needs an extra library)
#par(mar=c(0,0,0,0))
#plot(p, maxpixels=50000, axes=FALSE, legend=FALSE, col=brewer.pal('Set1', n=5))

writeRaster(p, "TerClass_sub.asc", overwrite=TRUE)

# variable importance: all contribute about the same... good
#importance(rf)
#varImpPlot(rf)

# customized plot (needs an extra library)
#par(mar=c(0,0,0,0))
#plot(p, maxpixels=50000, axes=FALSE, legend=FALSE, col=brewer.pal('Set1', n=5))

writeRaster(p, "TerClass.asc")
