library(spgrass6)
library(rgdal)
library(raster)
library(cluster)
library(randomForest)

loc <- initGRASS("/home/epilib/Envs/env1/grass-7.1.svn", home=tempdir())
execGRASS('g.proj',epsg=as.integer(32633), location='utm_wgs84_33N')
execGRASS('g.mapset', mapset='PERMANENT', location='utm_wgs84_33N') 

## create a raster for use in grass
execGRASS('r.in.gdal', input="subbathy.tif",
          flags=c('overwrite','o','e'), output='elev')                                  # o= override the prjection parameters, and e=extend the location

# should we aggregate the bathy at this point?

## set the region to the extent of the file
execGRASS('g.region', raster='elev', res='20', flags=c('a','p'))                             # align cells with region extent

## calculate the smooth raster
execGRASS('r.neighbors', input='elev', output='avg', size=as.integer(9),
          method='average', flags='overwrite')

execGRASS('r.neighbors', input='elev', output='min', size=as.integer(9),
          method='minimum', flags='overwrite')

execGRASS('r.neighbors', input='elev', output='max', size=as.integer(9),
          method='maximum', flags='overwrite')

#Elevation-relief ratio  
execGRASS('r.mapcalc', expression = 'er=1.0*(avg-min)/(max-min)', flags='overwrite') 

gc()       

#slope, and normalized slope, using the raster package
execGRASS('r.out.gdal', type='Float64', input="elev", output='elev.tif', flags='overwrite')
elev <- raster('elev.tif')
slope<-terrain(elev,opt=c('slope'),unit='degrees')
xslope<-slope/74.82501              #divide by max slope

#Profile curvature
execGRASS('r.param.scale', input='elev', output='profc', size=as.integer(11),
          slope_tolerance=0.1, curvature_tolerance=0.0001, method='profc', exponent=0.0, zscale=1.0, flags='overwrite')

#Cross-sectional curvature 
execGRASS('r.param.scale', input='elev', output='crosc', size=as.integer(11),
          slope_tolerance=0.1, curvature_tolerance=0.0001, method='crosc', exponent=0.0, zscale=1.0, flags='overwrite')

#Minimum curvature
execGRASS('r.param.scale', input='elev', output='minic', size=as.integer(11),
          slope_tolerance=0.1, curvature_tolerance=0.0001, method='minic', exponent=0.0, zscale=1.0, flags='overwrite')

#Maximum curvature
execGRASS('r.param.scale', input='elev', output='maxic', size=as.integer(11),
          slope_tolerance=0.1, curvature_tolerance=0.0001, method='maxic', exponent=0.0, zscale=1.0, flags='overwrite')

#Longitudinal curvature  
execGRASS('r.param.scale', input='elev', output='longc', size=as.integer(11),
          slope_tolerance=0.1, curvature_tolerance=0.0001, method='longc', exponent=0.0, zscale=1.0, flags='overwrite')                                                                                

#export, import
execGRASS('r.out.gdal', type='Float64', input="er", output='er.tif', flags='overwrite')
er <- raster('er.tif')

execGRASS('r.out.gdal', type='Float64', input="profc", output='profc.tif', flags='overwrite')
profc <- raster('profc.tif')

execGRASS('r.out.gdal', type='Float64', input="crosc", output='crosc.tif', flags='overwrite')
crosc <- raster('crosc.tif')

execGRASS('r.out.gdal', type='Float64', input="minic", output='minic.tif', flags='overwrite')
minic <- raster('minic.tif')

execGRASS('r.out.gdal', type='Float64', input="maxic", output='maxic.tif', flags='overwrite')
maxic <- raster('maxic.tif')

execGRASS('r.out.gdal', type='Float64', input="longc", output='longc.tif', flags='overwrite')
longc <- raster('longc.tif')

# stack
s <- stack(er,xslope,profc,crosc,minic,maxic,longc)
names(s) <- c('er','xslope','profc','crosc','minic','maxic','longc')

# we can't operate on the entire set of cells,
# so we use one of the central concepts from statistics: sampling
# sample 10000 random points
s.r <- as.data.frame(sampleRegular(s, 10000))

# clara() function: need to remove NA from training set
s.r <- na.omit(s.r)

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

writeRaster(p, "TerClass.asc")
