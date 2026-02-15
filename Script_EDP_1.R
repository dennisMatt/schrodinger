#install.packages("terra")
library(terra)
library(betareg)
library(distributions3)
library(dplyr)
library(raster)
library(gdistance)
library(sf)
library(spatstat)
library(spatstat.geom)
library(RColorBrewer)



################### Set foraging values


##########Woodland

woodWeight=0.8


##########################Grass


grassWeight<-0.5



#########################Water



waterWeight<-0.3



####################urban


urbanWeight<-0.2


#############################Nesting values (from Gardner et al.)



woodWeightNesting<-0.9


grassWeightNesting<-0.05


waterWeightNesting<-0


urbanWeightNesting<-0.1

#simulated habitat layers

simRasters<-rast("simRasters.tif")

names(simRasters)=c("woodland","grassland","water","urban")

plot(simRasters)


#init list for Boolean class values and assign values corresponding to movement cost
valList<-list()

for(i in 1:ncell(simRasters$woodland)){
  print(i)
  woodVal<-simRasters$woodland[i]
  grassVal<-simRasters$grassland[i]
  waterVal<-simRasters$water[i]
  urbanVal<-simRasters$urban[i]
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

boolean<-simRasters$woodland

boolean[]<-unlist(valList)

#Plot as cost values

plot(boolean,
         col=c("forestgreen",
               "green",
               "blue",
               "grey"),legend=T)



####################weighted cost layer

costList<-list()

for(i in 1:ncell(simRasters$woodland)){
  print(i)
  woodCost<-simRasters$woodland[i]
  grassCost<-simRasters$grassland[i]*3
  waterCost<-simRasters$water[i]*5
  urbanCost<-simRasters$urban[i]*8
  costVals<-sum(woodCost,grassCost,waterCost,urbanCost)
  costList[[i]]<-costVals
}

#Pass to new raster layer

costLayer<-simRasters$woodland


costLayer[]<-unlist(costList)


cols <- brewer.pal(7, "BuPu")


pal <- colorRampPalette(cols)


plot(costLayer,box=F,axes=F,col = pal(20),
     plg=list(legend=c("Cost"), x="right"))
text(x=390, y=140, "Cost", srt=0, cex=1.5, xpd=NA, pos=4)


############Categorical raster
library(rasterVis)

bool<-as.factor(boolean)


tar<-levels(bool)[[1]]
tar$landcover<-c("woodland (cost=1)", "grassland (cost=3)","wetland (cost=5)","urban (cost=8)")

levels(bool)<-tar


levels(bool) <- data.frame(id=c(1,3,5,8), element=c("woodland (cost=1)",
                                            "grassland (cost=3)",
                                            "wetland (cost=5)",
                                            "urban (cost=8)"))
plot(bool,axes=F,box=F,
     col=c("forestgreen","green",
      "blue","grey"), plg=list(x=390,130,cex=1.2),mar=c(7,7,7,9))



##########################Boolean PC metric and Least cost path


cost<-raster(bool)


wood<-simRasters$woodland

wood[wood<0.5]<-0


#Habitat patches
patches.i<-patches(wood,direction=8,zeroAsNA=TRUE)
patches.i<-as.polygons(patches.i,dissolve=TRUE)

      #calculate area
      patches.i$Area<-expanse(patches.i)
      
      #site centroids as nodes
      sites<-SpatialPoints(crds(centroids(patches.i)))
      
      #init matrix for cost distance values
      costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))

      #transition matrix
      land_cost <- transition(cost, transitionFunction = function(x) 1/mean(x), 8)
     
      #get number of patches
      n.patch<-nrow(costMat)

      #get least cost patch
     
      lcp <- shortestPath(land_cost,sites[1],sites,output="SpatialLines")
        
     
      distMat <- apply(costMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
      
      
      # set alpha which determines colonization probability of hypothetical species
      alpha= -log(0.05) / 100

      # init empty matrix for adjacency matrix
      A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))

      # negative exponential of colonization kernel x distance to get probability
      A.prob <- exp(-alpha * distMat)

      # set diag to one (patch is connected to itself)
      diag(A.prob) <- 1

      # final matrix for connectivity graph
      A.prob <- as.matrix(A.prob)
     
      graph.Aprob <- igraph::graph_from_adjacency_matrix(A.prob, mode="undirected", weighted=T)

      ## calculate modified PC metric

      # calculate all shortest paths between nodes
      pstar.mat <- distances(graph.Aprob, weights= -log(igraph::E(graph.Aprob)$weight))

      # back-transform to probabilities of connectedness
      pstar.mat <- exp(-pstar.mat)

      # get study area in m2
      AL <- ncell(wood)*res(wood)[1]^2

      # get area vector
      area <- patches.i$Area
      # sum areas from the vector
      areaSum <- sum(area)
    
      # get product of all patch areas ij and multiply by probabilities above
      PCmat <- outer(area,area) * pstar.mat
      
      # sum above
      pcMatSum <- sum(PCmat)

      # divide by total area of the study squared to get the EHI metric
      PC <- (sqrt(pcMatSum) / as.numeric(AL))*100
      

      plot(bool,axes=F,box=F,
           col=c("forestgreen","green",
                 "blue","grey"), plg=list(x=390,130,cex=1.2),mar=c(7,7,7,9))



      plot(lcp,add=T,col="red",lwd=3)



#####################Create variance layers

cellValsWoodSim<-list()

fixCellsWood<-list()
#Initial pass to set zeroes to some small number
for(i in 1:ncell(wood)){
  cell.i<-wood[i]
  cell.i
  if(cell.i==0){
    cell.i=0.001}
  fixCellsWood[[i]]<-cell.i

}

wood[]<-unlist(fixCellsWood)

#generate variance values
for(i in 1:ncell(wood)){
  cell.i<-wood[i]
  cell.q<-((sqrt(cell.i))-((cell.i^2)))*0.1
 
  if(cell.q>(cell.i/2)){
    cell.q=cell.i/2
  }
  if(cell.q<=0){
    cell.q=0.00001}


  cellValsWoodSim[[i]]<-cell.q
}

#add to new raster
valWoodSim<-unlist(cellValsWoodSim)

woodVar<-wood
woodVar[]<-valWoodSim

meansBetaWood<-values(wood)
varBetaWood<-values(woodVar)

#create beta distribution object from means and variance
X_betaWood <- Beta01(
  mu =meansBetaWood[1:507],
  phi = meansBetaWood[1:507]*(1-meansBetaWood[1:507])/(varBetaWood[1:507])-1,
  p0 = 0,
  p1 = 0
)


###########################GRASS SIMULATION


grass=simRasters$grassland

cellValsGrassSim<-list()

for(i in 1:ncell(grass)){
  cell.i<-grass[i]
  cell.q<-((sqrt(cell.i))-((cell.i^2)))*0.1

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


water=simRasters$water


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


###########################URBAN SIMULATION


urban=simRasters$urban


cellValsUrbanSim<-list()


for(i in 1:ncell(urban)){
  cell.i<-urban[i]
  cell.q<-((sqrt(cell.i))-((cell.i^2)))*0.1
  
  if(cell.q>(cell.i/2)){
    cell.q=cell.i/2
  }
  if(cell.q<=0){
    cell.q=0.00001}
  
  
  cellValsUrbanSim[[i]]<-cell.q
}


valUrbanSim<-unlist(cellValsUrbanSim)


urbanVar<-wood
urbanVar[]<-valUrbanSim

meansBetaUrban<-values(urban)
varBetaUrban<-values(urbanVar)



X_betaUrban <- Beta01(
  mu =meansBetaUrban[1:507],
  phi = meansBetaUrban[1:507]*(1-meansBetaUrban[1:507])/(varBetaUrban[1:507])-1,
  p0 = 0,
  p1 = 0
)


## function to simulate random variables, generate functional habitat patches,
#and compute connectivity metric

funBetaHab<-function(x){
  print(x)
  #random variable woodland distribution
  qWood<-random(X_betaWood,1)

  wood.i<-wood

  wood.i[]<-qWood

  
  #grass

  qGrass<-random(X_betaGrass,1)

  grass.i<-grass

  grass.i[]<-qGrass

  

  ##water
  qWater<-random(X_betaWater,1)

  water.i<-water

  water.i[]<-qWater
  

  ##urban

  qUrban<-random(X_betaUrban,1)

  urban.i<-urban

  urban.i[]<-qUrban

  ##Negative Urban effects

  #urban copy
  urb2<-extend(urban.i,c(30,30),fill=0)
 
  #calculate negative urban neighbourhood effects
  #edge effect radius
  radiusUrban<-80
  
  #neighbourhood weights matrix
  nPixE<-round(radiusUrban/res(urban.i)[1])
  nPixE<-(nPixE*2)+1
  
  matBaseUrban<-matrix(1:nPixE^2,nrow=nPixE,ncol=nPixE)
  
  #get focal cell
  x<-ceiling(ncol(matBaseUrban)/2)
  y<-ceiling(nrow(matBaseUrban)/2)
  
  focalCellUrban<-matBaseUrban[x,y]
  indFocalUrban<-which(matBaseUrban==focalCellUrban,arr.ind = T)
  
  #compute distances
  #init distance list
  distListUrban<-list()
  
  for(i in 1:nPixE^2){
    ind.i=which(matBaseUrban==i,arr.ind=T)
    diffX<-abs(ind.i[1,1]-indFocalUrban[1,1])*res(urban.i)[1]
    diffY<-abs(ind.i[1,2]-indFocalUrban[1,2])*res(urban.i)[1]
    
    dist.i<-sqrt(diffX^2+diffY^2)
    distListUrban[[i]]<-dist.i
    
  }
  #init new matrix
  distMatUrban<-matBaseUrban
  #add distance values
  distMatUrban[]<-unlist(distListUrban)
  #set cells outside search radius to zero
  distMatUrban[distMatUrban>radiusUrban]=NA
  #calculate edge effect decsay rate based on Dennis et al. 2024
  logN <- coercedLogistic(res(urban.i)[1]^2, fullEdgeEffectArea=1000000)
  #get constant for kernel shape based on max 80 metre edge effect (for urban)
  alphaUrban=log(logN)/80
  
  #compute cell edge effect values
  matKernelUrban<-exp(alphaUrban * distMatUrban)
  
  #set focal cell to zero
  matKernelUrban[x,y]<-0
  matKernelUrban[!is.na(matKernelUrban)] =1/length(matKernelUrban[!is.na(matKernelUrban)])
  
  
  
  #sum edge effects from all surrounding cells
  focalEdgeUrban<-focal(extend(urban.i,c(60,60),fill=0),matKernelUrban,fun=sum)
  focalEdgeUrban<-crop(focalEdgeUrban,ext(urban.i))
  
  
  ############################
  #build cost layer
  woodCostWeighted<-wood.i
  grassCostWeighted<-grass.i*3
  urbanCostWeighted<-urban.i*8
  waterCostWeighted<-water.i*5


  type1Cost<-woodCostWeighted+grassCostWeighted+urbanCostWeighted+waterCostWeighted
  
  #add additional cost from urban neighbourhood effect
  costFocal<-(focalEdgeUrban*8)
  
  #final cost
  type1CostEdge<-type1Cost+costFocal
  
  ####################Mask to potential patches (mutivariate habitat)
  #habitat values
  maskPatch<-(wood.i*woodWeightNesting)+(grass.i*grassWeightNesting)+(water.i*waterWeightNesting)+(urban.i*urbanWeightNesting)

  #binary habitat mask 
  maskPatchOne<-maskPatch

  maskPatchOne[maskPatchOne>=0.5]=1
  maskPatchOne[maskPatchOne!=1]=0

  #zero habitat layer as movement origins/destinations
  maskPatchZero<-maskPatch

  
  maskPatchZero[maskPatchZero<0.5]<-1
  maskPatchZero[maskPatchZero!=1]=0

  #cost layer
  costClip<-type1CostEdge*maskPatchZero
  
  #coast distance layer
  costComp<-costDist(costClip,target=0)
  

  # negative exponential of foraging kernel x distance to get probability

  radiusForage=100
  
  alphaForage= -log(0.05) / radiusForage
  
  #Kernel
  A.prob <- exp(-alphaForage * costComp)
  
  #kernel for woodland (minus patch to avoid double counting resources)
  A.probWood<-A.prob
  A.probWood[A.probWood==1]<-0

  #Availability probability for each resource type
  urbanProb<-A.prob*urban.i

  grassProb<-A.prob*grass.i

  waterProb<-A.prob*water.i
  
  woodProb<-A.probWood*wood.i
  
  
  
  #neighbourhood weights matrix to sum all available resources for each cell
  nPixF<-round(radiusForage/res(grass.i)[1])
  nPixF<-(nPixF*2)+1
  
  #buiild weights matrix
  matBaseForage<-matrix(1:nPixF^2,nrow=nPixF,ncol=nPixF)
  
  #get focal cell
  x<-ceiling(ncol(matBaseForage)/2)
  y<-ceiling(nrow(matBaseForage)/2)
  
  focalCellForage<-matBaseForage[x,y]
  indFocalForage<-which(matBaseForage==focalCellForage,arr.ind = T)
  #compute distances
  distListForage<-list()
  
  for(i in 1:nPixF^2){
    ind.i=which(matBaseForage==i,arr.ind=T)
    diffX<-abs(ind.i[1,1]-indFocalForage[1,1])*res(grass.i)[1]
    diffY<-abs(ind.i[1,2]-indFocalForage[1,2])*res(grass.i)[1]
    
    dist.i<-sqrt(diffX^2+diffY^2)
    distListForage[[i]]<-dist.i
    
  }
  #init new matrix
  focWin<-matBaseForage
  #add distance values
  focWin[]<-unlist(distListForage)
  #set cells outside search radius to NA
  focWin[focWin>radiusForage]=NA
  
  
  
  focWin[!is.na(focWin)]=1/length(focWin[!is.na(focWin)])
  
  #set focal cell to zero
  focWin[x,y]<-0
  
  
  #sum neighbourhood effects from all surrounding cells
  
  #extend to avoid losing raster edge
  coverProbUrban<-focal(extend(urbanProb,c(30,30),fill=0),focWin,fun=sum)

  coverProbUrban<-crop(coverProbUrban,wood.i)

  coverProbGrass<-focal(extend(grassProb,c(30,30),fill=0),focWin,fun=sum)
  
  coverProbGrass<-crop(coverProbGrass,wood.i)
  
  coverProbWood<-focal(extend(woodProb,c(30,30),fill=0),focWin,fun=sum)

  coverProbWood<-crop(coverProbWood,wood.i)
  
  #final foraging values
  forage<-coverProbGrass*grassWeight+coverProbUrban*urbanWeight+coverProbWood*woodWeight
  
  
  #make sure min is zero
  forage[forage<0]<-0

  #clip to original habitat cells
  
  maskPatch<-maskPatch*maskPatchOne
  
  forage<-forage*maskPatchOne
  
  #final habitat layer accounting for positive and negative neighbourhood effects
  r.All<-maskPatch+forage-focalEdgeUrban
  
  
  #Create patches
  r.All[r.All<0]<-0
  r.All[r.All>1]<-1

  r.All[r.All<0.5]<-0

  # get study area in m2
  AL <- ncell(r.All)*res(r.All)[1]^2


  #########################

  #gdistance needs raster objects (not spatraster)
  type1CostEdge<-raster(type1CostEdge)

  #get patches
  patches.i<-patches(r.All,direction=8,zeroAsNA=TRUE)
  patches.i<-as.polygons(patches.i,dissolve=TRUE)
 
  #catch special cases of connectivity where number of patches is zero or one
  if(nrow(patches.i)==0){
  RFH=0}else{
  if(nrow(patches.i)==1){
   RFH=sqrt(expanse(patches.i))/AL
  }else{

  #get continuous habitat amount
  contHab<- extract(r.All,patches.i,method="simple",fun="sum")

  contHab<-contHab[2]*(res(r.All)[2]^2)
  patches.i$Area<-contHab
  
  #get site centroids as network nodes
  sites<-SpatialPoints(crds(centroids(patches.i)))

  #init cost distance matrix
  costMat <- matrix(0, nrow=nrow(patches.i), ncol=nrow(patches.i))

  #transition matrix
  land_cost <- transition(type1CostEdge, transitionFunction = function(x) 1/mean(x), 8)

  #get number of patches
 
  
    costMat <- gdistance::costDistance(land_cost, sites, sites)
   
    l <- shortestPath(land_cost, sites, sites, output="SpatialLines")
    
    
  
  


  
  # make sure all elements in distance matrix are numeric 
  distMat <- apply(costMat, MARGIN=1, FUN=as.numeric) 
  

  # set alpha which determines colonization probability of the hypothetical species
  alpha= -log(0.05) / 100

  # init empty matrix for adjacency matrix
  A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))

  # negative exponential of colonization kernel x distance get probability
  A.prob <- exp(-alpha * distMat)

  # set diag to zero
  diag(A.prob) <- 1

  # final matrix for connectivity graph
  A.prob <- as.matrix(A.prob)

  # final matrix for connectivity graph
  graph.Aprob <- graph_from_adjacency_matrix(A.prob, mode="max", weighted=T)

  ## calculate FHI

  # calculate all shortest paths between nodes
  pstar.mat <- distances(graph.Aprob, weights= -log(E(graph.Aprob)$weight))

  # back-transform to probabilities of connectedness
  pstar.mat <- exp(-pstar.mat)



  # get area vector
  area <- patches.i$Area
  # sum areas from the vector
  areaSum <- sum(area)
 
  # get product of all patch areas ij and multiply by probabilities above
  PCmat <- outer(area,area) * pstar.mat
  
  # sum above
  pcMatSum <- sum(PCmat)

  # divide by total area of the study squared to get the PC metric
  RFH <- (sqrt(pcMatSum) / as.numeric(AL))*100

  }
  }
  
  #collect metrics in a data frame
  df<-data.frame(sum(expanse(patches.i)/10000),cbind(nrow(patches.i),sum(patches.i$Area,na.rm=T)/10000),RFH)
  colnames(df)<-c("BoolArea","nPatches","Area","RFH")
  
  #collate shortest paths as spatVector lines
  lcps<-vect(l)
  
  #plot(type1CostEdge)
  #plot(lcps,add=T)
  resList<-list(df,lcps)
  return(resList)

}


#iterate over N random variables
iters=1:10000

resSim<-lapply(iters,FUN=funBetaHab)

#access all shortest paths
sv<-lapply(resSim,'[[',2)

#convert to vector
sv<-vect(sv)

#make spatial then simple feature
veSP<-as(sv,"Spatial")

#sf
sfVect<-st_as_sf(veSP)

#remove any empty entries
sfVect<-sfVect[!st_is_empty(sfVect),]

#sf line object
spatLine<-spatstat.geom::as.psp(sf::st_geometry(sfVect))

#line density 
p<-density.psp(spatLine,sigma=2,method="FFT")

plot(p) 

############################
 
#access data frame and plot metric values
dfs<-lapply(resSim,'[[',1)

dfBind<-do.call("rbind",dfs)

hist(dfBind$RFH,breaks=20)






