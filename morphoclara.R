library(rgrass7)
library(rgdal)
library(raster)

    # Input:
    #    elevation: grid used as input
    #    resolution: (FALSE, Int) manualy set the resolution, FALSE use native dataset resolution
    #    win_neighbors:  int, default:9 - Size in pixel for the neighbors moving window operator.
    #                    (used to compute: average, minimun, maximum elevation)
    #    win_param_scale: int, default:15 - Size in pixel for the param_scale moving window operator.
    #                    (used to compute: profile, cross, minimum, longitudinal curvatures)
    #    curvature_tolerance: float, default: 0.0001 (*)
    #    slope_tolerance: float, default: 0.1 (*)
    #    exponent: float, default: 0.0 (*)
    #    zscale: float, default: 1.0 (*)
    #    Remove: boolean, If TRUE remove intermediate products. default FALSE
    #    
    #    (*) : more info at http://grass.osgeo.org/grass71/manuals/r.param.scale.html
    #
    # Output:
    #    brick raster. Layers: 
    #    'er','xslope','profc','crosc','minic','maxic','longc'

morphoclara <- function(elevation, resolution=FALSE, win_neighbors=9, win_param_scale=15, curvature_tolerance=0.0001, slope_tolerance=0.1, exponent=0.0, zscale=1.0, remove=FALSE ){
if (resolution) {
    execGRASS('g.region',  raster=elevation, res=toString(resolution), flags='p') 
    } else if (!resolution) {
    execGRASS('g.region', raster=elevation, flags='p')
} 
print('average elevation')
execGRASS('r.neighbors', input=elevation, output='avg', size=as.integer(win_neighbors),
          method='average', flags='overwrite')
print('minimum elevation')
execGRASS('r.neighbors', input=elevation, output='min', size=as.integer(win_neighbors),
          method='minimum', flags='overwrite')
print('maximum elevation')
execGRASS('r.neighbors', input=elevation, output='max', size=as.integer(win_neighbors),
          method='maximum', flags='overwrite')
print('normalized elevation')
execGRASS('r.mapcalc', expression = 'er=1.0*(avg-min)/(max-min)', flags='overwrite') 
print('slope')
execGRASS('r.slope.aspect', elevation=elevation, slope='slope', flags='overwrite')
minmax <- read.table(textConnection(execGRASS('r.info', map='slope', flags=c('r'), intern=TRUE)), 
                     header=FALSE, sep="=")
maxslope = minmax[2,2]
expression <- paste('xslope=slope/',maxslope)
execGRASS('r.mapcalc', expression = expression, flags='overwrite') 
print('Profile curvature')
execGRASS('r.param.scale', input=elevation, output='profc', size=as.integer(win_param_scale),
          slope_tolerance=as.numeric(slope_tolerance), curvature_tolerance=as.numeric(curvature_tolerance), 
          method='profc', exponent=as.numeric(exponent), zscale=as.numeric(zscale), flags='overwrite')
print('Cross-sectional curvature')
execGRASS('r.param.scale', input=elevation, output='crosc', size=as.integer(win_param_scale),
          slope_tolerance=as.numeric(slope_tolerance), curvature_tolerance=as.numeric(curvature_tolerance), 
          method='crosc', exponent=as.numeric(exponent), zscale=as.numeric(zscale), flags='overwrite')
print('Minimum curvature')
execGRASS('r.param.scale', input=elevation, output='minic', size=as.integer(win_param_scale),
          slope_tolerance=as.numeric(slope_tolerance), curvature_tolerance=as.numeric(curvature_tolerance), 
          method='minic', exponent=as.numeric(exponent), zscale=as.numeric(zscale), flags='overwrite')
print('Maximum curvature')
execGRASS('r.param.scale', input=elevation, output='maxic', size=as.integer(win_param_scale),
          slope_tolerance=as.numeric(slope_tolerance), curvature_tolerance=as.numeric(curvature_tolerance), 
          method='maxic', exponent=as.numeric(exponent), zscale=as.numeric(zscale), flags='overwrite')
print('Longitudinal curvature')
execGRASS('r.param.scale', input=elevation, output='longc', size=as.integer(win_param_scale),
          slope_tolerance=as.numeric(slope_tolerance), curvature_tolerance=as.numeric(curvature_tolerance), 
          method='longc', exponent=as.numeric(exponent), zscale=as.numeric(zscale), flags='overwrite')                                                                                       
rasters <- readRAST(c('er','xslope','profc','crosc','minic','maxic','longc'), cat=NULL, 
                    ignore.stderr = get.ignore.stderrOption(), NODATA=NULL,
                    mapset=NULL, useGDAL=get.useGDALOption(), close_OK=TRUE, drivername="GTiff",
                    driverFileExt=NULL, return_SGDF=TRUE)
print('Raster stack')
s <- stack(rasters)
names(s) <- c('er','xslope','profc','crosc','minic','maxic','longc')
if (remove) {execGRASS('g.remove', type='raster', name=c('er','xslope','profc','crosc','minic','maxic','longc'), flags='f')
}
return(s)
}
