######################################################################################################################
# Fits a Conditional Inference Forest and calculates probability of occurrence of a sessile species                  #
# (for which P/A data is supplied, and predictors are available as raster datasets)                                  #
######################################################################################################################  

# inputs/requirements:
#
# predictor data as geotiffs or ascii files (stored in folder named "PredictorData")
#
# response data (one per species) as .txt or .csv with observations as rows, and the following columns:
#                                                           "RefNo": ID for the sampling reference,
#                                                           "PointID": ID for each individual observation,
#                                                           "Count" and/or "PA": count or Presence/Absence data
#                                                           "x_coord", "y_coord", Easting and Northing
# (stored in folder named "ResponseData")
#
# optinally, a shapefile for the region to be modelled (stored in folder named "Shapes")


# Libraries and paths -----------------------------------------------------

require(raster)
require(party)
require(pROC)

#setwd()=getwd()

# change as needed!!!
dataPath <- "G:/R_projects/DWCDistributionModel/PredictorData"
resPath <- "G:/R_projects/DWCDistributionModel/ResponseData"
shapePath <- "G:/R_projects/DWCDistributionModel/shapes"

outPath <- "G:/R_projects/DWCDistributionModel/Outputs"
dir.create(file.path(outPath))

# change as needed!!!
projection <- "+proj=utm +north +zone=33 +datum=WGS84 +ellps=WGS84 
  +towgs84=0,0,0"

# comment or uncomment as needed
#modelregion <- readOGR(dsn = shapePath, layer = "mask") 

# Read in predictor data --------------------------------------------------

## bathymetry and terrain

# change as needed!!!
pred1<-raster(file.path(dataPath,'bathy20_1.tif'))
projection(pred1)<-projection
#pred1<-mask(pred1,modelregion)

# change as needed!!!
pred2<-raster(file.path(dataPath,'TerClass.asc'),datatype=INT1U)
projection(pred2)<-projection
#pred2<-mask(pred2,modelregion)

## ocean climate

# change as needed!!!
pred3<-raster(file.path(dataPath,'Smax_ann.asc'))
projection(pred3)<-projection

# change as needed!!!
pred4<-raster(file.path(dataPath,'Smean_ann.asc'))
projection(pred4)<-projection

# change as needed!!!
pred5<-raster(file.path(dataPath,'Smin_ann.asc'))
projection(pred5)<-projection

# change as needed!!!
pred6<-raster(file.path(dataPath,'SPDmax_ann.asc'))
projection(pred6)<-projection

# change as needed!!!
pred7<-raster(file.path(dataPath,'SPDmean_ann.asc'))
projection(pred7)<-projection

# change as needed!!!
pred8<-raster(file.path(dataPath,'Tmax_ann.asc'))
projection(pred8)<-projection

# change as needed!!!
pred9<-raster(file.path(dataPath,'Tmean_ann.asc'))
projection(pred9)<-projection

# change as needed!!!
pred10<-raster(file.path(dataPath,'Tmin_ann.asc'))
projection(pred10)<-projection

# change as needed!!!
n.pred <- 10 # number of predictors

pred.names <- c(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8,pred9,pred10)

## crop all layers to same extent and resample to same resolution (credits: Ben K,
## http://stackoverflow.com/questions/20733555/how-to-create-a-raster-brick-with-rasters-of-different-extents)

extent_list <-list()
length(extent_list)<-n.pred

for(i in 1:n.pred){
  extent_list[[i]]<-extent(pred.names[[i]])
  }

# make a matrix out of it, each column represents a raster, rows the values
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=n.pred)
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")

# create an extent with the extrem values of your extent
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]),
                    min(matrix_extent[2,]), max(matrix_extent[4,]))

# the range of your extent in degrees
ranges<-apply(as.matrix(best_extent), 1, diff)
# the resolution of your raster (pick one) or add a desired resolution
reso<-res(pred1)
# dividing the range by your desired resolution gives you the number of rows and columns
nrow_ncol<-ranges/reso

s<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=projection)

# resample predictor layers
pred1<-resample(pred1, s)
pred2<-resample(pred2, s, method="ngb")

pred3<-resample(pred3, s)
pred4<-resample(pred4, s)
pred5<-resample(pred5, s)
pred6<-resample(pred6, s)
pred7<-resample(pred7, s)
pred8<-resample(pred8, s)
pred9<-resample(pred9, s)
pred10<-resample(pred10, s)

# or (but not tested yet):
#library(oce) #needed?
#for (name in pred.names) {
#  load(paste(name, '.rda', sep='')
#       d <- get(name)
#       eval(parse(text=paste('rm(', name, ')')))
#       d <- resample(d,s)
#       assign(name, d)
#}

# brick the predictors
predbrick <-  brick(pred1,
                    pred2,
                    pred3,
                    pred4,
                    pred5,
                    pred6,
                    pred7,
                    pred8,
                    pred9,
                    pred10)

# Read in response data ---------------------------------------------------

# change as needed!!!
respdata<-read.table(file.path(resPath,"L_pertusa.csv"),header=T, sep=",")
#respdata<-read.table(file.path(resPath,"P_resedaef.csv"),header=T, sep=",")


## define Training and Evaluation data, where E data should not break up reference stations

respdata$split <- "T"
p <- 25 # proportion of data that should be used for Evaluation, as percent (minimum acceptable)

# the number of E data depends on p and on the average number of observations per station
propE <- as.integer((dim(respdata)[1]*p)/100)
use <- sample(unique(respdata$RefNo),propE/mean(table(respdata$RefNo)))
respdata$split[respdata$RefNo %in% use] <- "E"

percentE <- (table(respdata$split)[1]*100)/table(respdata$split)[2]
percentE # percent of the datapoints being used for evaluation

percentE_P <- (dim(respdata[respdata$PA=="1" & respdata$split=="E",])[1]/dim(respdata[respdata$PA=="1",])[1])*100
percentE_P # percent of presence datapoints being used for evaluation


## define Training and Evaluation data, completely randomly

respdata$split <- "T"

# now the number of E datapoints only depends on p
propE <- as.integer((dim(respdata)[1]*p)/100)
use <- sample(unique(respdata$PointID),propE)
respdata$split[respdata$PointID %in% use] <- "E"

percentE <- (table(respdata$split)[1]*100)/table(respdata$split)[2]
percentE # percent of the datapoints being used for evaluation

percentE_P <- (dim(respdata[respdata$PA=="1" & respdata$split=="E",])[1]/dim(respdata[respdata$PA=="1",])[1])*100
percentE_P # percent of presence datapoints being used for evaluation


# Extract predictor values ------------------------------------------------

coordinates(respdata)<-~x_coord+y_coord
proj4string(respdata) <- CRS(projection)

v <- cbind(respdata@data, extract(predbrick, respdata))
#v$TerClass<-as.factor(v$TerClass)


# Fitting the model -------------------------------------------------------

CIF.pa <- cforest(PA~Smax_ann + Smean_ann, 
                  data=v[which(v$split=="T"),],
                  control = cforest_unbiased(ntree=1000,mtry = 2))

sort(varimp(CIF.pa))

#check whether you get the same results
#with a different random seed before interpreting
#the variable importance ranking!
#If the ranking of even the top scoring predictor
#variables depends on the choice of the random
#seed, increase the number of trees (argument
#ntree in randomForest and cforest_control).


# Model performance -------------------------------------------------------

t<-v[which(v$split=="E"),6:dim(v)[2]]                                                                                
pred.prob<-predict(CIF.pa,t,OOB=TRUE)                                        # Predicted probability for E points


o<-v[which(v$split=="E"),4]                                                  # Observed presences for E points                           

crosscheck<-cbind(pred.prob,o) 
colnames(crosscheck)<-c('predicted','observed')                                                                                
accuracy<-auc(observed~predicted,crosscheck)


z<-respdata@data[which(respdata@data$split=="T"),which(colnames(respdata@data)=="PA")]
prevalence<-length(which(z>0))/length(z)

## is the accuracy greater than the no information rate (which is the prevalence)? 
roc.model<-roc(observed~predicted,crosscheck)
crosscheck.shuff<-cbind(crosscheck[,1],sample(crosscheck[,2]))
colnames(crosscheck.shuff)<-c('predicted','observed')
roc.null<-roc(observed~predicted,crosscheck.shuff)
roc.test(roc.model, roc.null) 


# Other info to report ------------------------------------------------------------------

CIF.pa
table(respdata$split)
sort(varimp(CIF.pa))
accuracy
summary(z)
prevalence
