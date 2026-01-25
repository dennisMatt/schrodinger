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

# Linear interpolation function
#   rescale value c from range a-b to range y-z
###
lerp <- function(c, a, b, y, z) {
  (c - a) * (z - y) / (b - a) + y
}


logistic <- function(x, L=1.0, b=0.0, k, x0) {
  L / (1.0 + exp(-k*(x-x0))) + b
}

###
# This makes a call to the logistic function and then scales the result to enforce a value of 1 for the
#   stated maximum value ([max]), 0 for 0 and 0.5 where x = x0 (maximum val/2). This is achieved using a
#   linear  interpolation function to scale either from 0-0.5 (where x<0.5) or 0.5-[max] (otherwise). This
#   is necessary as scaling from 0-max can displace the case so that y!=0.5 where x=x0 if the adjustment
#   required at x=0 and x=[max] are not equal.
#
# Params:
#   `x` is the value for which you want to return the y
#   `fullEdgeEffectArea` is the x value that should be coerced to 1
#
###
coercedLogistic <- function(x, fullEdgeEffectArea) {
  
  # work out order of magnitude and convert to scale factor
  # TODO: this is an arbitrary parameter
  k <- 1 / (10^(floor(log(fullEdgeEffectArea, 10))-1))
  
  # x0 is the midpoint of the scale
  x0 <- fullEdgeEffectArea*0.5
  
  # get the 'raw' logistic value for the current x
  vx = logistic(x, k=k, x0=x0)
  
  # stretch it either up towards 1 or down towards 0 to coerce the scale to (x=0)=0, (x-x0)=0.5, (x=fullEdgeEffectArea)=1
  if (vx < 0.5) lerp(vx, logistic(0, k=k, x0=x0), 0.5, 0, 0.5) else lerp(vx, 0.5, logistic(fullEdgeEffectArea, k=k, x0=x0), 0.5, 1)
}



#set alpha and beta parameters of foraging values

alphaWood<-4.926017992

betaWood<-1.829663826

rs<-sample.int(100,1)

set.seed(rs)

dWood<-dbeta(seq(0.001,1,0.001), shape1 = alphaWood, shape2 = betaWood)

dWood<-data.frame(cbind(dWood,seq(0.001,1,0.001)))

colnames(dWood)<-c("density","prob")

#get final weight

woodWeight<-dWood$prob[which.max(dWood$density)]



##################################################################
######### Grassland 

alphaGrass<-0.019093878


betaGrass<-1.890293878

dGrass<-dbeta(seq(0.001,1,0.001), shape1 = alphaGrass, shape2 = betaGrass)

dGrass<-data.frame(cbind(dGrass,seq(0.001,1,0.001)))

colnames(dGrass)<-c("density","prob")



grassWeight<-dGrass$prob[which.max(dGrass$density)]



#################################################
#shrub


alphashrub<-17.74833333


betashrub<-5.401666667

dshrub<-dbeta(seq(0.001,1,0.001), shape1 = alphashrub, shape2 = betashrub)

dshrub<-data.frame(cbind(dshrub,seq(0.001,1,0.001)))

colnames(dshrub)<-c("density","prob")

shrubWeight<-dshrub$prob[which.max(dshrub$density)]


#######################################
#############Water



alphaWater<-0.019093878




betaWater<-1.890293878




dWater<-dbeta(seq(0.001,1,0.001), shape1 = alphaWater, shape2 = betaWater)

dWater<-data.frame(cbind(dWater,seq(0.001,1,0.001)))

colnames(dWater)<-c("density","prob")


waterWeight<-dWater$prob[which.max(dWater$density)]



##########################urban

alphaUrban<-4.981049563

betaUrban<-18.2638484



dUrban<-dbeta(seq(0.001,1,0.001), shape1 = alphaUrban, shape2 = betaUrban)

dUrban<-data.frame(cbind(dUrban,seq(0.001,1,0.001)))

colnames(dUrban)<-c("density","prob")


urbanWeight<-dUrban$prob[which.max(dUrban$density)]


#############################Nesting weights


woodWeightNesting<-0.78125

###############################


##################################################################
######### Grass


grassWeightNesting<-0

#################################################
#shrub


shrubWeightNesting<-0.43333



#######################################
#############Water



waterWeightNesting<-0




#####################################
#########################Urban




urbanWeightNesting<-0.21111






#########################load in layers containing membership and variance from Random Forest classifier

woodMean<-rast("woodMeanClipped.tif")
woodVar<-rast("woodVarClipped.tif")
grassMean<-rast("grassMeanClipped.tif")
grassVar<-rast("grassVarClipped.tif")
shrubMean<-rast("shrubMeanClipped.tif")
shrubVar<-rast("shrubVarClipped.tif")
waterMean<-rast("waterMeanClipped.tif")
waterVar<-rast("waterVarClipped.tif")
urbanMean<-rast("urbanMeanClipped.tif")
urbanVar<-rast("urbanVarClipped.tif")
classif<-rast("rfClasses.tif")

##########Generate Beta Distributions

#get memberships and variance of woodland cover
woodMeans<-values(woodMean)
woodVars<-values(woodVar)

#build data frame
woodData<-data.frame(cbind(woodMeans,woodVars))

# build beta distribution
X_wood <- Beta01(
  mu =woodData$compositeGB_1,
  phi = woodData$compositeGB_1*(1-woodData$compositeGB_1)/(woodData$compositeGB_1.1)-1,
  p0 = 0,
  p1 = 0
)

#glassland cover beta
grassMeans<-values(grassMean)
grassVars<-values(grassVar)

grassMeans[grassMeans==0]<-0.001
grassVars[grassVars==0]<-0.0001


grassData<-data.frame(cbind(grassMeans,grassVars))


X_grass <- Beta01(
  mu =grassData$compositeGB_1,
  phi = grassData$compositeGB_1*(1-grassData$compositeGB_1)/(grassData$compositeGB_1.1)-1,
  p0 = 0,
  p1 = 0
)


#shrub cover beta

shrubMeans<-values(shrubMean)
shrubVars<-values(shrubVar)

shrubMeans[shrubMeans==0]<-0.001
shrubVars[shrubVars==0]<-0.0001

shrubData<-data.frame(cbind(shrubMeans,shrubVars))



X_shrub <- Beta01(
  mu =shrubData$compositeGB_1,
  phi = shrubData$compositeGB_1*(1-shrubData$compositeGB_1)/(shrubData$compositeGB_1.1)-1,
  p0 = 0,
  p1 = 0
)


#water cover beta

waterMeans<-values(waterMean)
waterVars<-values(waterVar)

waterMeans[waterMeans==0]<-0.001
waterVars[waterVars==0]<-0.0001


waterData<-data.frame(cbind(waterMeans,waterVars))

waterData[waterData$compositeGB_1==0,]<-0.0001


X_water <- Beta01(
  mu =waterData$compositeGB_1,
  phi = waterData$compositeGB_1*(1-waterData$compositeGB_1)/(waterData$compositeGB_1.1)-1,
  p0 = 0,
  p1 = 0
)

#urban beta

urbanMeans<-values(urbanMean)
urbanVars<-values(urbanVar)

urbanMeans[urbanMeans==0]<-0.001
urbanVars[urbanVars==0]<-0.0001


urbanData<-data.frame(cbind(urbanMeans,urbanVars))


X_urban <- Beta01(
  mu =urbanData$compositeGB_1,
  phi = urbanData$compositeGB_1*(1-urbanData$compositeGB_1)/(urbanData$compositeGB_1.1)-1,
  p0 = 0,
  p1 = 0
)



## simulate random variables

p<-seq(from=0.001,to=0.999,0.001)

#function to simulate landscape composition and calculate patch density, area and connectivity 

funBeta<-function(x){
  #generate random variable and select from ppf of woodland membership values
  qWood<-quantile(X_wood,sample(p,1), drop = FALSE, elementwise = FALSE)
  
  
  wood.i<-woodMean
  
  wood.i[]<-qWood
  
 
  
  #generate random variable and select from ppf of grassland membership values
  qGrass<-quantile(X_grass,sample(p,1), drop = TRUE, elementwise = FALSE)
  
  grass.i<-woodMean
  
  grass.i[]<-qGrass
  
  #generate random variable and select from ppf of shrub membership values
  qShrub<-quantile(X_shrub,sample(p,1), drop = TRUE, elementwise = FALSE)
  
  shrub.i<-woodMean
  
  shrub.i[]<-qShrub
  

  
  
  #generate random variable and select from ppf of water membership values
  qWater<-quantile(X_water, sample(p,1), drop = TRUE, elementwise = FALSE)
  
  water.i<-woodMean
  
  water.i[]<-qWater
 
  #
  #generate random variable and select from ppf of urban membership values
  qUrban<-quantile(X_urban, sample(p,1), drop = TRUE, elementwise = FALSE)
  
  urban.i<-woodMean
  
  urban.i[]<-qUrban
  

  
  
  
  
  #####################Calculate proportions per cell
  
  
  sumVals<-(wood.i+shrub.i+grass.i+
              urban.i+water.i)
  wood.i<-wood.i/(sumVals)
  
  shrub.i<-shrub.i/sumVals
  
  grass.i<-grass.i/sumVals
  
  water.i<-water.i/sumVals
  
  urban.i<-urban.i/sumVals
  
  #calculate negative urban neighbourhood effects
  #dge effect radius
  radiusE<-76
  
  #neighbourhood weights matrix
  nPixE<-round(radiusE/res(wood.i)[1])
  nPixE<-(nPixE*2)+1
  
  
  matBaseE<-matrix(1:nPixE^2,nrow=nPixE,ncol=nPixE)
  
  x<-ceiling(ncol(matBaseE)/2)
  y<-ceiling(nrow(matBaseE)/2)
  #get focal cell
  focalCellE<-matBaseE[x,y]
  indFocalE<-which(matBaseE==focalCellE,arr.ind = T)
  
  distListE<-list()
  #compute distances
  for(i in 1:nPixE^2){
    ind.i=which(matBaseE==i,arr.ind=T)
    diffX<-abs(ind.i[1,1]-indFocalE[1,1])*res(wood.i)[1]
    diffY<-abs(ind.i[1,2]-indFocalE[1,2])*res(wood.i)[1]
    
    dist.i<-sqrt(diffX^2+diffY^2)
    distListE[[i]]<-dist.i
    
  }
  #init new matrix
  distMatE<-matBaseE
  #add distance values
  distMatE[]<-unlist(distListE)
  #set cells outside search radius to zero
  distMatE[distMatE>radiusE]=NA
  #calculate edge effect decsay rate based on Dennis et al. 2024
  logN <- coercedLogistic(res(wood.i)[1]^2, fullEdgeEffectArea=1000000)
  #get constant for kernel shape based on max 76 metre edge effect (for urban)
  alphaE=log(logN)/76
  
  #compute cell edge effect values
  matKernelE<-exp(alphaE * distMatE)
  #set focal cell to zero
  matKernelE[x,y]<-0
  
  #sum edge effects from all surrounding cells
  focalD5Sum<-focal(extend(urban.i,c(60,60),fill=0),matKernelE,fun=sum)
  focalD5Sum<-crop(focalD5Sum,ext(wood.i))
  
  ############repeat for grassland cover
  
  radiusG<-30
  nPixG<-round(radiusG/res(wood.i)[1])
  nPixG<-(nPixG*2)+1
  
  
  matBaseG<-matrix(1:nPixG^2,nrow=nPixG,ncol=nPixG)
  
  xG<-ceiling(ncol(matBaseG)/2)
  yG<-ceiling(nrow(matBaseG)/2)
  
  focalCellG<-matBaseG[xG,yG]
  indFocalG<-which(matBaseG==focalCellG,arr.ind = T)
  
  distListG<-list()
  
  for(i in 1:nPixG^2){
    ind.i=which(matBaseG==i,arr.ind=T)
    diffX<-abs(ind.i[1,1]-indFocalG[1,1])*res(wood.i)[1]
    diffY<-abs(ind.i[1,2]-indFocalG[1,2])*res(wood.i)[1]
    
    dist.i<-sqrt(diffX^2+diffY^2)
    distListG[[i]]<-dist.i
    
  }
  
  distMatG<-matBaseG
  distMatG[]<-unlist(distListG)
  #distMatG
  #distMat[x,y]<-1
  distMatG[distMatG>radiusG]=NA
  #plot(rast(distMatG))
  #res(wood.i)[1]^2
  
  alphaG=log(logN)/30
  
  matKernelG<-exp(alphaG * distMatG)
 
  matKernelG[xG,yG]<-0
  sum(matKernelG,na.rm=T)
 
  focalD5Grass<-focal(extend(grass.i,c(30,30),fill=0),matKernelG,fun=sum)
  focalD5Grass<-crop(focalD5Grass,ext(wood.i))
  
  ###################
  
  
  
  ############################Create cost layers
  #get mask of motorway (high cost)
  m62<-vect("M62.shp")
  
  m62Cost<-rasterize(m62,field=40,wood.i,background=0)
  m62Hab<-rasterize(m62,field=0,wood.i,background=1)
  
  #get contribite of motorway to functional cost
  focalM62Cost<-focal(extend(m62Cost,c(30,30),fill=0),matKernelE,fun="sum")
 
  focalM62Cost<-crop(focalM62Cost,wood.i)
  
  #assign cost weights to other cover types
  woodCostWeighted<-wood.i
  grassCostWeighted<-grass.i*10
  urbanCostWeight<-urban.i*5
  waterCostWeight<-water.i*10
  shrubCostWeight<-shrub.i*2
  
  
  
  type1Cost<-woodCostWeighted+grassCostWeighted+urbanCostWeight+waterCostWeight+shrubCostWeight
 
  costFocal<-(focalD5Sum*5)
  #add in neighbourhood effects on cost and the motorway cover
  type1CostEdge<-type1Cost+costFocal+m62Cost+focalM62Cost+focalD5Grass
  
  #remove motorway from urban foraging potential
  urbanFin<-urban.i*m62Hab
 
  
  ####################Mask to potential patches
  
  maskPatch<-(wood.i*woodWeightNesting)+(grass.i*grassWeightNesting)+(shrub.i*shrubWeightNesting)+(water.i*waterWeightNesting)+(urbanFin*urbanWeightNesting)
  
  
  rastEdge<-maskPatch
  rastEdge[rastEdge<0.5]=0
  rastEdge[rastEdge!=0]=1
  
  #init binary raster 0-1
  maskPatchOne<-maskPatch
  
  #init mask of habitat patches
  maskPatchOne[maskPatchOne>=0.5]=1
  maskPatchOne[maskPatchOne!=1]=0
 
  #combine adjacent habitat cells
  allPatches<-patches(maskPatchOne,directions=8,zeroAsNA=T)
  
  #polygonize patches for cost distance analysis
  allPatches.i<-as.polygons(allPatches,dissolve=TRUE)
  allPatches.i$A<-expanse(allPatches.i)
  
  #subset patches >= 0.5 ha
  allPatches.i<-allPatches.i[allPatches.i$A>=5000]
  
  #get number of patches to loop over for cost distance computation
  nPatch<-nrow(allPatches.i)
  
  allPatches.i$ID<-1:nrow(allPatches.i)
  allPatches.r<-rasterize(allPatches.i,wood.i,field="ID",background=0)
  
  
  ############################################
  #############Foraging kernels
  ##############################################
  
  #assume foraging distance of 250 metres
  radiusF<-250
  
  #get constant to set kernel shape
  alpha= -log(0.05) / radiusF
  
  #object to iterate of cell IDs
  patchVals<-values(allPatches.r)
  
  #foraging values per cell
  fWeights<-(grass.i*grassWeight+water.i*waterWeight+
               urban.i*urban.i+wood.i*woodWeight)
  
  #init list to store total foraging scores per focal cell 
  forageList<-list()
  
  for(i in unique(patchVals[!is.na(patchVals)])){
    if(i!=0){
      r.i<-allPatches.r
      
    
      
      r.i[r.i!=i]<-0
      r.i[r.i!=0]<-1
      r.iCost<-r.i-1
      r.iCost[r.iCost==-1]<-1
   
      costClip.i<-(type1CostEdge)*r.iCost
      
      costComp.i<-costDist(costClip.i,target=0)
      
      # this sets shape of neg. exp. curve
      distForage = -log(0.05) / radiusF
      # this creates the raster kernel
      forageRast <- exp(-distForage * costComp.i)
     
      # add zeros back in so we can sum rasters later
     
      forageRast[forageRast==1]<-0
     
      forage.i<-(forageRast*fWeights)

      #Normalise bu N cells within foraging distance
      win=matrix(1,ncol=round(radiusF/res(forage.i)[1])+1,nrow=round(radiusF/res(forage.i)[1])+1)
      win=win/sum(win)
    
      focalForageSum<-focal(extend(forage.i,10,fill=0),w=win,fun="sum")
      focalForageSum<-crop(focalForageSum,wood.i)
     
      
      forageList[[i]]<-focalForageSum
      
    }
  }
  
  #remove empty list entries
  forageList<-forageList[lapply(forageList,length)>0]
  #sum foraging scores
  forageSum<-sum(rast(forageList))
  
  #keep only habitat cells
  forage<-forageSum*maskPatchOne
  
  maskPatch<-maskPatch*maskPatchOne
  
  #get final habitat values by remove negative and adding positive neighbourhood effects
  r.All<-maskPatch-focalD5Sum-focalD5Grass+forage
  r.All[r.All<0]<-0
  r.All[r.All>1]<-1
  
  r.All[r.All<0.5]<-0

  # get study area in m2
  AL <- ncell(r.All)*res(r.All)[1]^2
  
  #create raster object for use with gdistance package
  type1CostEdge<-raster(type1CostEdge)
  
  #get patches of adjacent cells
  patches.i<-patches(r.All,direction=8,zeroAsNA=TRUE)
  patches.i<-as.polygons(patches.i,dissolve=TRUE)
  
  #subset patches >=0.5 ha
  patches.i$A<-expanse(patches.i)
  
  patches.i<-patches.i[patches.i$A>=5000,]
  
  #get habitat ammount based on habitat values per cell
  Interior<- extract(r.All,patches.i,method="simple",fun="sum")
  
  Interior<-Interior[2]*(res(r.All)[2]^2)
  
  patches.i$Area<-Interior
  
  
  #get patch centroids for least cost path delineation
  sites<-SpatialPoints(crds(centroids(patches.i)))
  
  #create cost distance matrix
  costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))
  #init transition layer
  land_cost <- transition(1/type1CostEdge, transitionFunction=mean, 8)
  
  #number of patches to compute
  n.patch<-nrow(costMat)
  
  #init list for least cost paths 
  lineList<-list()
  for (i in 1:n.patch) {
    c <- gdistance::costDistance(land_cost, sites[i], sites)
    costMat[i,] <- c
    l <- shortestPath(land_cost, sites[i,], sites,output="SpatialLines")
    lineList[[i]]<-l
    
  }
  
  
  distMat<-as.matrix(costMat, nrow=nrow(patches.i), ncol=nrow(patches.i))
  # make sure all elements are numeric here
  distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) 
  
  #get mean distances
  meanDist<-mean(distMat)
  # set alpha to determine colonization probability of the species (max dispersal capacity as 10km after Dennis et al. 2024)
  alpha= -log(0.05) / 10000
  
  # init empty matrix for adjacency matrix
  A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))
  
  # negative exponential of colonization kernel x distance get probability
  A.prob <- exp(-alpha * distMat)
  
  # set diag to one
  diag(A.prob) <- 1
  
  # final matrix for connectivity graph
  A.prob <- as.matrix(A.prob)

  
  
  # adjacency graph
  graph.Aprob <- igraph::graph_from_adjacency_matrix(A.prob, mode="undirected", weighted=T)
  
  ## calculate Connectivity metric (FHI)

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
 
  pcMatSum <- sum(PCmat)
  
  # divide by total area of the study squared to get the PC metric
  EHI <- (sqrt(pcMatSum) / as.numeric(AL))*100
  
  #print(paste("rep = ",x))
  print(paste("EHI=",EHI))
  df<-data.frame(sum(expanse(patches.i))/10000,nPatch,cbind(nrow(patches.i),sum(patches.i$Area,na.rm=T)/10000),EHI,meanDist)
  colnames(df)<-c("BoolArea","allPatches","nPatches","Area","EHI","meanDist")
  print(x)
  
  resList<-list(df,lineList)
  return(resList)
}


iters=1:10000

resColne<-lapply(iters,FUN=funBeta)



sv<-lapply(resColne,'[[',2)

sv<-vect(sv)

writeVector(sv,paste("all",sample.int(10000,1),"_lcps.shp",sep=""))


dfs<-lapply(resColne,'[[',1)

dfBind<-do.call("rbind",dfs)

write.csv(dfBind,paste("all",sample.int(10000,1),"_all.csv",sep=""))





