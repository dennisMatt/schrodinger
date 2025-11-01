#install.packages("terra")
library(terra)
library(betareg)
library(distributions3)
library(dplyr)
library(raster)
library(gdistance)
library(sf)
library(parallel)
library(doParallel)
library(spatstat)
library(spatstat.geom)

setwd("/Users/mattdennis/downloads/btp")

rs<-sample.int(100,1)

set.seed(rs)
woodWeight=0.8


##################################################################
#########Upland Grass


grassWeight<-0.5





#######################################
#############Water



waterWeight<-0.3



##########################urban




urbanWeight<-0.2


#############################Nesting

#plot(dWoodNesting$prob,dWoodNesting$density)

woodWeightNesting<-0.9


grassWeightNesting<-0.05


#######################################
#############Water


waterWeightNesting<-0


#####################################
#########################Urban


urbanWeightNesting<-0.1

#####################Exp filter functions


filter_create <- function(r = NULL,
                          radius = NULL,
                          type = c("exp_decay", "bartlett", "circle", "threshold_decay",
                                   "gaussian_decay", "Gauss", "rectangle")[1],
                          zoi_limit = 0.05,
                          half_life = NULL,
                          zoi_hl_ratio = NULL,
                          sigma = NULL,
                          min_intensity = 0.01,
                          max_dist = 5000,
                          normalize = FALSE,
                          divisor = 1,
                          round_vals = NULL,
                          save_txt = FALSE,
                          save_format = c("GRASS_rmfilter", "raw")[1],
                          save_folder = NULL,
                          save_file = NULL,
                          parallel = TRUE,
                          ...) {

  # check the input data class of r
  if(class(r) %in% c("RasterLayer", "RasterBrick", "RasterStack", "SpatRaster")) {
    res <- terra::res(r)[1]
  } else {
    if(is.numeric(r) & r > 0) {
      res <- r
    } else
      stop("'r' must be either an input raster map or a numeric value corresponding to the resolution of a raster.")

  }

  # apply function
  if(type == "exp_decay") {
    parms <- set_filt_exp_decay(radius = radius,
                                zoi_limit = zoi_limit,
                                res = res,
                                half_life = half_life,
                                zoi_hl_ratio = zoi_hl_ratio,
                                min_intensity = min_intensity,
                                max_dist = max_dist)
  }

  if(type %in% c("step", "threshold", "circle", "threshold_decay", "step_decay")) {
    parms <- set_filt_step(radius = radius, res = res)
  }

  if(type %in% c("bartlett", "batlett_decay", "tent_decay", "linear_decay")) {
    parms <- set_filt_bartlett(radius = radius, res = res)
  }

  if(type %in% c("rectangle", "box")) {
    parms <- set_filt_rectangle(radius = radius, res = res)
  }

  if(type %in% c("Gauss", "gauss", "gaussian", "gaussian_decay")) {
    parms <- set_filt_gassian_decay(radius = radius,
                                    zoi_limit = zoi_limit,
                                    res = res,
                                    sigma = sigma,
                                    min_intensity = min_intensity,
                                    max_dist = max_dist)
  }

  # get parameters
  radius <- parms$radius
  radius_pix <- parms$radius_pix
  size_pix <- parms$size_pix

  # create distance matrix
  # distance in pixels to the central cell of the matrix
  dist_mat <- sqrt((matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = F) - (radius_pix + 1))^2+
                     (matrix(c(1:size_pix), nrow = size_pix, ncol = size_pix, byrow = T) - (radius_pix + 1))^2)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # apply function
  if(type == "exp_decay") {
    dist_mat <- exp(-parms$lambda * dist_mat)
  }

  if(type %in% c("step", "threshold", "circle", "threshold_decay", "step_decay")) {
    dist_mat <- 1 * (dist_mat*res <= radius)
  }

  if(type %in% c("bartlett", "batlett_decay", "tent_decay", "linear_decay")) {
    dist_mat <- pmax((1 + parms$lambda * dist_mat), 0)
  }

  if(type %in% c("rectangle", "box")) {
    dist_mat[] <- 1
  }

  if(type %in% c("Gauss", "gauss", "gaussian", "gaussian_decay")) {
    dist_mat <- exp(-parms$lambda * dist_mat**2)
  }
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # normalize
  if(normalize)
    # dist_mat <- dist_mat/sum(dist_mat[1+radius_pix,])
    dist_mat <- dist_mat/sum(dist_mat)

  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  # round decimals
  if(!is.null(round_vals))
    if(round_vals >= 0) dist_mat <- round(dist_mat, round_vals)
  # image(dist_mat)
  # plot(terra::rast(dist_mat))

  if(save_txt) {
    # save matrix outside R for use within GRASS GIS
    oneimpact::filter_save(filt = dist_mat, radius = radius, type = type,
                           save_format = save_format, save_folder = save_folder,
                           save_file = save_file, parallel = parallel,
                           divisor = divisor, separator = " ")

  }

  dist_mat
}

set_filt_exp_decay <- function(radius = NULL,
                               zoi_limit = 0.05,
                               half_life = NULL,
                               res = 100,
                               zoi_hl_ratio = NULL,
                               min_intensity = 0.01,
                               max_dist = 200){

  # define lambda depending on the input parameter
  if(!is.null(radius)) {

    # define radius in terms on number of pixels
    radius <- radius/res

    if(is.null(zoi_hl_ratio)) {
      lambda <- log(1/zoi_limit) / radius
    } else {
      half_life <- radius/zoi_hl_ratio
      lambda <- log(2)/half_life
    }

  } else {

    if(!is.null(half_life)) {
      # define radius or half life, depending on which is given as input
      half_life <- half_life/res
      lambda <- log(2)/half_life
    } else {
      stop("Either both 'radius' and 'zoi_limit' must be specified, or both 'half_life' and 'zoi_hl_ratio'.")
    }
  }

  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*radius)))
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}


set_filt_step <- function(radius, res){

  # define radius and size (diameter)
  radius_pix <- ceiling(radius/res)
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_rectangle <- function(radius, res){

  # define radius and size (diameter)
  radius_pix <- floor(radius/res)
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = NULL))
}

set_filt_bartlett <- function(radius, res){

  # define radius and size (diameter)
  radius_pix <- ceiling(radius/res)
  size_pix <- 2*radius_pix + 1
  # define beta (beta = -b/a or beta = -1/radius)
  lambda <- -1/(radius/res)

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}

set_filt_gassian_decay <- function(radius = NULL,
                                   zoi_limit = 0.05,
                                   res = 100,
                                   sigma = NULL,
                                   min_intensity = 0.01,
                                   max_dist = 50000){


  if(!is.null(radius)) {
    # define radius in terms on number of pixels
    radius <- radius/res
    lambda = log(1/zoi_limit) / (radius**2)
  } else {
    if(!is.null(sigma)) {
      # define sigma in terms on number of pixels
      sigma <- sigma/res
      lambda = 1/(2*sigma**2)
    } else {
      stop("Either 'radius' or 'sigma' must be specified.")
    }

  }

  # tmp <- exp(-lambda * c(0:round(half_life*6))/half_life)
  # define radius and size (diameter)
  tmp <- exp(-lambda * c(0:round(2*radius))**2)
  radius_pix <- min(which(tmp < min_intensity)[1], round(max_dist/res))
  size_pix <- 2*radius_pix + 1

  return(list(radius = radius, radius_pix = radius_pix, size_pix = size_pix, lambda = lambda))
}



#########################layers

simRasters<-rast("simRasters.tif")

plot(simRasters)

valList<-list()

for(i in 1:ncell(simRasters$layer.1)){
  print(i)
  woodVal<-simRasters$layer.1[i]
  grassVal<-simRasters$layer.2[i]
  waterVal<-simRasters$layer.3[i]
  urbanVal<-simRasters$layer.4[i]
  allVals<-c(woodVal,grassVal,waterVal,urbanVal)
  maxVal<-which.max(unlist(allVals))
  if(maxVal==1){
    valList[[i]]<-1
  }
  if(maxVal==2){
    valList[[i]]<-3
  }
  if(maxVal==3){
    valList[[i]]<-5
  }
  if(maxVal==4){
    valList[[i]]<-8
  }
}

boolean<-simRasters$layer.1

boolean[]<-unlist(valList)

par(xpd=FALSE)

dev.off()

par(c("mar", "mai"))
par(mar=c(0,2,2,2))
dev.off()

plot(boolean,legend=F)

plot(boolean,
         col=c("forestgreen",
               "green",
               "blue",
               "grey"),legend=F)






####################weighted cost layer

costList<-list()

for(i in 1:ncell(simRasters$layer.1)){
  print(i)
  woodCost<-simRasters$layer.1[i]
  grassCost<-simRasters$layer.2[i]*3
  waterCost<-simRasters$layer.3[i]*5
  urbanCost<-simRasters$layer.4[i]*8
  costVals<-sum(woodCost,grassCost,waterCost,urbanCost)
  costList[[i]]<-costVals
}




costLayer<-simRasters$layer.1


costLayer[]<-unlist(costList)

colours()

cols <- brewer.pal(7, "BuPu")
#colsSD <- brewer.pal(7, "Greens")

pal <- colorRampPalette(cols)
#palSD <- colorRampPalette(colsSD)

plot(costLayer,box=F,axes=F,col = pal(30),
     plg=list(legend=c("Cost"), x="right"))
text(x=390, y=140, "Cost", srt=0, cex=1.5, xpd=NA, pos=4)
title("B",
      adj  = 0,
      line = -6)


plot(costLayer)



library(rasterVis)

bool<-as.factor(boolean)


tar<-levels(bool)[[1]]
tar$landcover<-c("woodland (cost=1)", "grassland (cost=3)","wetland (cost=5)","urban (cost=8)")
tar
levels(bool)<-tar

plot(bool)

levels(bool) <- data.frame(id=c(1,3,5,8), element=c("woodland (cost=1)",
                                            "grassland (cost=3)",
                                            "wetland (cost=5)",
                                            "urban (cost=8)"))
plot(bool,axes=F,box=F,
     col=c("forestgreen","green",
      "blue","grey"), plg=list(x=390,130,cex=1.2),mar=c(7,7,7,9))


text(x=12, y=135, "B", srt=0, cex=1.2, xpd=NA, pos=4)

##########################Boolean LCP

type1CostEdge<-raster(bool)
plot(type1CostEdge)

r.All<-simRasters$layer.1
r.All[r.All<0.5]<-NA
plot(r.All)

woodland<-bool
woodland[woodland!=1]<-NA
plot(woodland)
patches.i<-patches(woodland,direction=8,zeroAsNA=TRUE)
patches.i<-as.polygons(patches.i,dissolve=TRUE)
plot(patches.i,add=T)


      Interior<- extract(r.All,patches.i,method="simple",fun="sum")

      Interior<-Interior[2]*(res(r.All)[2]^2)
      patches.i$Area<-expanse(patches.i)
      head(patches.i)

  
      sites<-SpatialPoints(crds(centroids(patches.i)))
      
      costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))

      land_cost <- transition(1/type1CostEdge, transitionFunction=mean, 8)


      n.patch<-nrow(costMat)


      lineList=list()
      for (i in 1:n.patch) {
        c <- gdistance::costDistance(land_cost, sites[i], sites)
        costMat[i,] <- c
        l <- shortestPath(land_cost, sites[i,], sites,output="SpatialLines")
        lineList[[i]]<-l
        #print(i)
      }
      
      distMat<-as.matrix(costMat, nrow=nrow(patches.i), ncol=nrow(patches.i))
      distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
      distMat<-distMat*10
      
      
      # set alpha which determines colonization probability of the species
      alpha= -log(0.05) / 1000

      # init empty matrix for adjacency matrix
      A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))

      # negative exponential of colonization kernel x distance get probability
      A.prob <- exp(-alpha * distMat)

      # set diag to zero
      diag(A.prob) <- 1

      # final matrix for connectivity graph
      A.prob <- as.matrix(A.prob)
      #A.prob
      # final matrix for connectivity graph

      # final matrix for connectivity graph
      graph.Aprob <- igraph::graph_from_adjacency_matrix(A.prob, mode="undirected", weighted=T)

      ## calculate RHI


      # calculate all shortest paths between nodes
      pstar.mat <- distances(graph.Aprob, weights= -log(igraph::E(graph.Aprob)$weight))

      # back-transform to probabilities of connectedness
      pstar.mat <- exp(-pstar.mat)

      # get study area in m2
      AL <- ncell(r.All)*res(r.All)[1]^2

      # get area vector
      area <- patches.i$Area
      # sum areas from the vector
      areaSum <- sum(area)
    
      # get product of all patch areas ij and multiply by probabilities above
      PCmat <- outer(area,area) * pstar.mat
      
      # sum above
      pcMatSum <- sum(PCmat)

      # divide by total area of the study squared to get the PC metric
      EHI <- (sqrt(pcMatSum) / as.numeric(AL))*100
      

      plot(bool,axes=F,box=F,
           col=c("forestgreen","green",
                 "blue","grey"), plg=list(x=390,130,cex=1.2),mar=c(7,7,7,9))



      plot(l,add=T,col="red",lwd=3)






#####################


cellValsWoodSim<-list()

fixCellsWood<-list()

for(i in 1:ncell(wood)){
  cell.i<-wood[i]
  cell.i
  if(cell.i==0){
    cell.i=0.001}
  fixCellsWood[[i]]<-cell.i

}

wood[]<-unlist(fixCellsWood)


for(i in 1:ncell(wood)){
  cell.i<-wood[i]
  cell.q<-((sqrt(cell.i))-((cell.i^2)))*0.1
  #cell.q<-cell.q-0.1
  #if(cell.q>cell.i){
  # cell.q=cell.i/2}
  #if(cell.q>0.19){
  #cell.q=0.19}
  if(cell.q>(cell.i/2)){
    cell.q=cell.i/2
  }
  if(cell.q<=0){
    cell.q=0.00001}


  cellValsWoodSim[[i]]<-cell.q
}


valWoodSim<-unlist(cellValsWoodSim)

woodVar<-wood
woodVar[]<-valWoodSim

meansBetaWood<-values(wood)
varBetaWood<-values(woodVar)

library(betareg)
library(distributions3)

X_betaWood <- Beta01(
  mu =meansBetaWood[1:507],
  phi = meansBetaWood[1:507]*(1-meansBetaWood[1:507])/(varBetaWood[1:507])-1,
  p0 = 0,
  p1 = 0
)

X_betaWood

#apply function and create data frame



###########################GRASS SIMULATION


grass=simRasters$layer.2

cellValsGrassSim<-list()


for(i in 1:ncell(grass)){
  cell.i<-grass[i]
  cell.q<-((sqrt(cell.i))-((cell.i^2)))*0.1
  #cell.q<-cell.q-0.1
  #if(cell.q>cell.i){
  # cell.q=cell.i/2}
  #if(cell.q>0.19){
  #cell.q=0.19}
  if(cell.q>(cell.i/2)){
    cell.q=cell.i/2
  }
  if(cell.q<=0){
    cell.q=0.00001}


  cellValsGrassSim[[i]]<-cell.q
}


valGrassSim<-unlist(cellValsGrassSim)

grassVar<-wood
grassVar[]<-valGrassSim

meansBetaGrass<-values(grass)
varBetaGrass<-values(grassVar)


X_betaGrass <- Beta01(
  mu =meansBetaGrass[1:507],
  phi = meansBetaGrass[1:507]*(1-meansBetaGrass[1:507])/(varBetaGrass[1:507])-1,
  p0 = 0,
  p1 = 0
)



###########################WATER SIMULATION


water=simRasters$layer.3


cellValsWaterSim<-list()


for(i in 1:ncell(water)){
  cell.i<-water[i]
  cell.q<-((sqrt(cell.i))-((cell.i^2)))*0.1

  if(cell.q>(cell.i/2)){
    cell.q=cell.i/2
  }
  if(cell.q<=0){
    cell.q=0.00001}


  cellValsWaterSim[[i]]<-cell.q
}


valWaterSim<-unlist(cellValsWaterSim)


waterVar<-wood
waterVar[]<-valWaterSim

meansBetaWater<-values(water)
varBetaWater<-values(waterVar)

X_betaWater <- Beta01(
  mu =meansBetaWater[1:507],
  phi = meansBetaWater[1:507]*(1-meansBetaWater[1:507])/(varBetaWater[1:507])-1,
  p0 = 0,
  p1 = 0
)
######################################

## simulate random variables

p<-seq(from=0.001,to=0.999,0.001)

#X_wood

funBeta<-function(x){


  nrep=x
  qWood<-random(X_betaWood,1)

  wood.i<-wood

  wood.i[]<-qWood

  #plot(wood.i)
  #grass

  qGrass<-random(X_betaGrass,1)

  grass.i<-grass

  grass.i[]<-qGrass

  #plot(grass.i)


  #########################water
  qWater<-random(X_betaWater,1)

  water.i<-water

  water.i[]<-qWater
  #plot(water.i)

  #########################urban

  qUrban<-random(X_betaUrban,1)

  urban.i<-urban

  urban.i[]<-qUrban



  #####################Negative Urban


  urb2<-extend(urban.i,c(30,30),fill=0)
 
  filterUrban<-filter_create(r=urb2,type="exp_decay",radius=80, normalize = TRUE)

  x<-ceiling(ncol(filterUrban)/2)
  y<-ceiling(nrow(filterUrban)/2)

  filterUrban[x,y]<-0

  focalD5Sum<-focal(urb2,w=filterUrban,fun="sum")

  #focalD5Sum[focalD5Sum<0.001]<-NA
  focalD5Sum<-crop(focalD5Sum,ext(wood.i))

  #plot(focalD5Sum)
  #plot(urban.i)
  ############################

  woodCostWeighted<-wood.i
  grassCostWeighted<-grass.i*5
  urbanCostWeight<-urban.i*10
  waterCostWeight<-water.i*20


  type1Cost<-woodCostWeighted+grassCostWeighted+urbanCostWeight+waterCostWeight
  #plot(type1Cost)
  costFocal<-(focalD5Sum*5)
  type1CostEdge<-type1Cost+costFocal

  #plot(type1CostEdge)

  ####################Mask to potential patches

  maskPatch<-(wood.i*woodWeightNesting)+(grass.i*grassWeightNesting)+(water.i*waterWeightNesting)+(urban.i*urbanWeightNesting)


  #plot(maskPatch)
  rastEdge<-maskPatch
  rastEdge[rastEdge<0.5]=0
  rastEdge[rastEdge!=0]=1

  maskPatchOne<-maskPatch

  maskPatchOne[maskPatchOne>=0.5]=1
  maskPatchOne[maskPatchOne!=1]=0

  #plot(maskPatchOne)

  maskPatches<-maskPatch

  maskPatches[maskPatches<0.5]<-0

  maskPatches<-patches(maskPatches,directions=8,zeroAsNA=T)
  #plot(maskPatches)
  maskPatches<-as.polygons(maskPatches,dissolve=T)
  #nrow(maskPatches)
  maskRast<-rasterize(maskPatches,wood.i,field=0,background=1)


  costClip<-type1CostEdge*maskRast

  costComp<-costDist(costClip,target=0)
  #costComp[costComp==0]<-1
  allPatches<-patches(maskPatchOne,directions=8,zeroAsNA=T)

  allPatches.i<-as.polygons(allPatches,dissolve=TRUE)
  #allPatches.i$A<-expanse(allPatches.i)
  #allPatches.i<-allPatches.i[allPatches.i$A>=5000]

  nPatch<-nrow(allPatches.i)

  allPatches.i$ID<-1:nrow(allPatches.i)
  allPatches.r<-rasterize(allPatches.i,wood.i,field="ID",background=0)
  #plot(allPatches)
  #####

  alpha= -log(0.05) / 100

  #plot(costComp)

  # negative exponential of colonization kernel x distance get probability

  A.prob <- exp(-alpha * costComp)
  #plot(A.prob)
  A.probWood<-A.prob
  A.probWood[A.probWood==1]<-0

  urbanProb<-(A.prob*urban.i)

  grassProb<-A.prob*grass.i

  waterProb<-A.prob*water.i
  woodProb<-A.probWood*wood.i
  #probWood<-A.prob*wood.i
  #plot(woodProb)

  focWin<-filter_create(A.prob,100,type="circle",normalize = TRUE)

  x<-ceiling(ncol(focWin)/2)
  y<-ceiling(nrow(focWin)/2)

  focWin[x,y]<-0

  coverProbUrban<-focal(extend(urbanProb,c(30,30),fill=0),focWin,fun=sum)

  coverProbUrban<-crop(coverProbUrban,wood.i)

  coverProbGrass<-focal(extend(grassProb,c(30,30),fill=0),focWin,fun=sum)
  coverProbGrass<-crop(coverProbGrass,wood.i)

  coverProbWood<-focal(extend(woodProb,c(30,30),fill=0),focWin,fun=sum)

  coverProbWood<-crop(coverProbWood,wood.i)

  patchVals<-values(allPatches.r)

  patchFocList<-list()
  for(i in unique(patchVals[!is.na(patchVals)])){
    if(i!=0){
    r.i<-allPatches.r
    #r.i[is.na(r.i)]<-0
    r.i[r.i!=i]<-0
    r.i[r.i!=0]<-1
    #plot(r.i)
    wood.ri<-r.i*woodProb
    #plot(wood.ri)
    foc.i<-focal(wood.ri,focWin,fun=sum,na.rm=T)
    foc.i<-foc.i*r.i
    #r.i
    #plot(foc.i)

    #print(i)
    patchFocList[[i]]<-foc.i
  }
  }

  #if(length(patchFocList>0)){
  patchFocals<-rast(patchFocList)
  patchSum<-sum(patchFocals)
  #plot(patchSum)
  #plot(coverProbWood)
  wood.j<-(coverProbWood-patchSum)
  wood.j[wood.j<0]=0
  #plot(wood.j)
  allPatches1<-allPatches.r
  allPatches1[allPatches1!=0]<-1
  wood.j<-wood.j*allPatches1
  wood.j<-(wood.j*woodWeight)
  #plot(wood.j)
  #plot(coverProbGrass)
  #extWood<-extract(wood.j,allPatches.i)
  #sum(extWood$focal_sum)
  wood.k<-wood.j
  forage<-coverProbGrass*grassWeight+coverProbUrban*urbanWeight+wood.j*woodWeight

  forage[forage<0]<-0

  #plot(forage)
  #plot(focalD5Sum)
  #library(RColorBrewer)
  #colours()


  forage<-forage*allPatches1
  maskPatch<-maskPatch*allPatches1
  #plot(maskPatch)
  r.All<-maskPatch+forage-focalD5Sum
  #plot(r.All)
  r.All[r.All<0]<-0
  r.All[r.All>1]<-1

  r.All[r.All<0.5]<-0

  # get study area in m2
  AL <- ncell(r.All)*res(r.All)[1]^2

  #plot(r.All,col=palSD(20),box=T,axes=T)
  #plot(maskPatch)
  #plot(patches.i,add=T)
  #########################


  type1CostEdge<-raster(type1CostEdge)


  patches.i<-patches(r.All,direction=8,zeroAsNA=TRUE)
  patches.i<-as.polygons(patches.i,dissolve=TRUE)

  #patches.i$A<-expanse(patches.i)

  #patches.i<-patches.i[patches.i$A>=5000,]

  #print(paste("N sites =",nrow(patches.i)))
  #patches.i

  if(nrow(patches.i)==0){
  EHI=0}else{
  if(nrow(patches.i)==1){
   EHI=sqrt(expanse(patches.i))/AL
  }else{


  Interior<- extract(r.All,patches.i,method="simple",fun="sum")

  Interior<-Interior[2]*(res(r.All)[2]^2)
  patches.i$Area<-Interior


  #print(paste("total area =",sum(patches.i$Area,na.rm=T)))
  sites<-SpatialPoints(crds(centroids(patches.i)))

  costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))

  land_cost <- transition(1/type1CostEdge, transitionFunction=mean, 8)


  n.patch<-nrow(costMat)
  #print(n.patch)
  lineList<-list()
  for (i in 1:n.patch) {
    c <- gdistance::costDistance(land_cost, sites[i], sites)
    costMat[i,] <- c
    l <- shortestPath(land_cost, sites[i,], sites, output="SpatialLines")
    lineList[[i]]<-l
    #print(i)
  }


  #writeRaster(rasHolder,paste("lcp",i,".tif",sep=""))
  distMat<-as.matrix(costMat, nrow=nrow(patches.i), ncol=nrow(patches.i))
  distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
  distMat<-distMat*10
  meanDist<-mean(distMat)
  # set alpha which determines colonization probability of the species
  alpha= -log(0.05) / 1000

  # init empty matrix for adjacency matrix
  A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))

  # negative exponential of colonization kernel x distance get probability
  A.prob <- exp(-alpha * distMat)

  # set diag to zero
  diag(A.prob) <- 1

  # final matrix for connectivity graph
  A.prob <- as.matrix(A.prob)

  # final matrix for connectivity graph
  graph.Aprob <- graph_from_adjacency_matrix(A.prob, mode="undirected", weighted=T)

  ## calculate RHI

  # calculate all shortest paths between nodes
  pstar.mat <- distances(graph.Aprob, weights= -log(E(graph.Aprob)$weight))

  # back-transform to probabilities of connectedness
  pstar.mat <- exp(-pstar.mat)



  # get area vector
  area <- patches.i$Area
  # sum areas from the vector
  areaSum <- sum(area)
  #areaSum/AL*100
  #outer(area,area)
  # get product of all patch areas ij and multiply by probabilities above
  PCmat <- outer(area,area) * pstar.mat
  #hist(PCmat)
  # sum above
  pcMatSum <- sum(PCmat)

  # divide by total area of the study squared to get the PC metric
  EHI <- (sqrt(pcMatSum) / as.numeric(AL))*100

  }
  }
  #print(paste("rep = ",x))
  #print(paste("EHI=",EHI))

  df<-data.frame(sum(expanse(patches.i)/10000),nPatch,cbind(nrow(patches.i),sum(patches.i$Area,na.rm=T)/10000),EHI,meanDist)
  colnames(df)<-c("BoolArea","allPatches","nPatches","Area","EHI","meanDist")
  #df
  lcps<-vect(lineList)
  #plot(lcps)
  print(nrep)
  resList<-list(df,lineList)
  return(resList)

}


iters=1:10000

resSim<-lapply(iters,FUN=funBeta)

sv<-lapply(resSim,'[[',2)

sv<-vect(sv)

veSP<-as(sv,"Spatial")

sfVect<-st_as_sf(veSP)

sfVect<-sfVect[!st_is_empty(sfVect),]

spatLine<-spatstat.geom::as.psp(sf::st_geometry(sfVect))

plot(spatLine)

p<-density.psp(spatLine,sigma=1,method="FFT")

plot(p)


############################



dfs<-lapply(resSim,'[[',1)

dfBind<-do.call("rbind",dfs)

hist(dfBind$EHI)








