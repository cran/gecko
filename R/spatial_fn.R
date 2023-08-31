#####GECKO - Geographical Ecology and Conservation Knowledge Online
#####Version 0.1.2 (2022-05-28)
#####By Vasco Branco, Pedro Cardoso, Luís Correia
#####Maintainer: vasco.branco@helsinki.fi
#####Changed from v0.1.1:
#####Replaced package raster functions with new ones from package terra

############################################################################
##################################PACKAGES##################################
############################################################################
#' @importFrom terra rast vect mask trim terrain as.points as.data.frame rasterize distance extract plot spatSample predict layerCor values sbar xmin xmax classify crds mean where.min xyFromCell
#' @importFrom grDevices dev.copy dev.off chull pdf
#' @importFrom sp SpatialPolygons Polygons Polygon
#' @importFrom graphics par text lines points
#' @importFrom stats as.dist dist prcomp
#' @importFrom utils data read.csv
#' @importFrom geosphere areaPolygon
#' @importFrom methods is
NULL
#> NULL 

############################################################################
##################################DATASETS##################################
############################################################################
#' Example data packaged with *gecko*
#' @description Load data included in the package. This includes *gecko.records*,
#' a matrix of longitude and latitude (two columns) occurrence records for
#' Hogna maderiana (Walckenaer, 1837); *gecko.range*, a SpatRaster object, as
#' defined by package terra, of the geographic range of Hogna maderiana
#' (Walckenaer, 1837); *gecko.layers*, a SpatRaster object with layers 
#' representing the average annual temperature, total annual precipitation,
#' altitude and landcover for Madeira Island
#' (Fick & Hijmans 2017, Tuanmu & Jetz 2014); and *worldborders* is a small vector
#' of world country borders. 
#' @param data Name of data in quotes. E.g.: `"gecko.records"`
#' If `NULL`, the example files will be listed.
#' @export
#' @examples
#' gecko.examples()
#' gecko.examples("gecko.range")
#' @source This function is inspired by `palmerpanguins::path_to_file()`
#' which in turn is based on `readxl::readxl_example()`.
gecko.examples <- function(data = NULL) {
  if (is.null(data)) {
    print(
      c(
       "gecko.records", "gecko.range", "gecko.layers", "worldborders" 
      )
    )
    return(NULL)
  } else {
    if(data == "gecko.records"){
      path = system.file(paste0("extdata/gecko.records.csv"), package = "gecko")
      out = utils::read.csv(path)
    } else if (data == "gecko.range") {
      path = system.file(paste0("extdata/gecko.range.tif"), package = "gecko")
      out = terra::rast(x = path)
    } else if (data == "gecko.layers") {
      path = system.file(paste0("extdata/gecko.layers.", c(1:4), ".tif"), package = "gecko")
      out = terra::rast(x = path)
    } else if (data == "worldborders") {
      path = system.file(paste0("extdata/worldborders"), package = "gecko")
      out = terra::vect(x = path)
    }
  }
  return(out)
}


##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################

#' Uniformize raster layers.
#' @description Crop raster layers to minimum size possible and uniformize NA values across layers.
#' @param layers Raster* object as defined by package raster.
#' @details Excludes all marginal rows and columns with only NA values and change values to NA if they are NA in any of the layers.
#' @return A Raster* object, same class as layers.
#' @examples data = gecko.examples("gecko.layers")
#' terra::plot(clean(data))
#' @export
clean <- function(layers){
  
  ##apply mask to have NAs everywhere where any layer has NAs
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  layers <- terra::mask(layers, maskLayer)
  
  ##crop by excluding external rows and columns with NAs only
  layers <- terra::trim(layers)
  
  return(layers)
}


#' Create eastness layer.
#' @description Create a layer depicting eastness based on an elevation layer.
#' @param dem RasterLayer object of elevation (a digital elevation model - DEM) as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples data = gecko.examples("gecko.layers")
#' terra::plot(create.east(data[[3]]))
#' @export
create.east <- function(dem){
  asp <- terra::terrain(dem, v = "aspect")
  return(sin(asp))
}


#' Create latitude layer.
#' @description Create a layer depicting latitude based on any other.
#' @param layers Raster* object as defined by package raster.
#' @details Using latitude (and longitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data = gecko.examples("gecko.layers")
#' terra::plot(create.lat(data[[1]]))
#' @export
create.lat <- function(layers){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
  x <- terra::as.points(layers)[,1:2]
  lat <- terra::rasterize(terra::crds(x), y = layers, values = terra::crds(x)[,1] )
  lat <- terra::mask(lat, layers)
  names(lat) <- "latitude"
  return(lat)
}


#' Create longitude layer.
#' @description Create a layer depicting longitude based on any other.
#' @param layers Raster* object as defined by package raster.
#' @details Using longitude (and latitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples data = gecko.examples("gecko.layers")
#' terra::plot(create.long(data))
#' @export
create.long <- function(layers){
  if(dim(layers)[3] > 1) {
    layers <- layers[[1]]
  }
  x <- terra::as.points(layers)[,1:2]
  long <- terra::rasterize(terra::crds(x), y = layers, values = terra::crds(x)[,2] )
  long <- terra::mask(long, layers)
  names(long) <- "longitude"
  return(long)
}


#' Create northness layer.
#' @description Create a layer depicting northness based on an elevation layer.
#' @param dem RasterLayer object of elevation (a digital elevation model - DEM) as defined by package raster.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return A RasterLayer object.
#' @examples data = gecko.examples("gecko.layers")
#' terra::plot(create.north(data[[3]]))
#' @export
create.north <- function(dem){
  asp <- terra::terrain(dem, v = "aspect")
  return(cos(asp))
}


#' Create distance layer.
#' @description Creates a layer depicting distances to records using the minimum, average, distance to the minimum convex polygon or distance taking into account a cost surface.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster to serve as model to create distance layer.
#' @param type text string indicating whether the output should be the "minimum", "average" or "mcp" distance to all records. "mcp" means the distance to the minimum convex polygon encompassing all records.
#' @details Using distance to records in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return A RasterLayer object.
#' @examples userpar <- par(no.readonly = TRUE) 
#' layers = gecko.examples("gecko.layers")
#' alt = layers[[3]]
#' records = gecko.examples("gecko.records")
#' par(mfrow=c(3,2))
#' terra::plot(alt)
#' points(records)
#' terra::plot(distance(records, alt))
#' terra::plot(distance(records, alt, type = "average"))
#' par(userpar)
#' @export
distance <- function(longlat, layers, type = "minimum"){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
    
  #layers[] <- 0
  
  layers_template = terra::classify(!is.na(layers), c(TRUE, 0))
  
  if(type == "average"){
    for (d in 1:nrow(longlat)) {
      # layers <- layers + raster::distanceFromPoints(layers, longlat[d,]) #  he distance unit is in meters if the coordinate reference system (crs) of the Raster* object is (+proj=longlat) or assumed to be if the crs is NA. In all other cases it is in the units defined by the crs (which typically is meters).
      layers <- c(
        layers,
        terra::distance(
          layers_template,
          terra::vect(longlat[d, ],
            geom = colnames(longlat),
            crs = terra::crs(layers)
          )
        )
      )
    }
    layers = layers[[2:dim(layers)[3]]]
    # No mean?
    layers = terra::mean(layers)
    
    # layers <- layers/nrow(longlat)
    names(layers) <- "average distance"
  } else {
    #layers <- raster::mask(raster::distanceFromPoints(layers, longlat), layers)
    layers <- terra::mask(
      terra::distance(
        layers,
        terra::vect(longlat,
          geom = colnames(longlat),
          crs = terra::crs(layers)
        )
      ),
      layers
    )
    names(layers) <- "minimum distance"
  }
  return(layers)
}


#' Spatial thinning of occurrence records.
#' @description Thinning of records with minimum distances either absolute or relative to the species range.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param distance Distance either in relative terms (proportion of maximum distance between any two records) or in raster units.
#' @param relative If TRUE, represents the proportion of maximum distance between any two records. If FALSE, is in raster units.
#' @param runs Number of runs
#' @details Clumped distribution records due to ease of accessibility of sites, emphasis of sampling on certain areas in the past, etc. may bias species distribution models.
#' The algorithm used here eliminates records closer than a given distance to any other record. The choice of records to eliminate is random, so a number of runs are made and the one keeping more of the original records is chosen.
#' @return A matrix of species occurrence records separated by at least the given distance.
#' @examples userpar <- par(no.readonly = TRUE)
#' records <- matrix(sample(100), ncol = 2)
#' par(mfrow=c(1,2))
#' graphics::plot(records)
#' records <- thin(records, 0.1)
#' graphics::plot(records)
#' par(userpar)
#' @export
thin <- function(longlat, distance = 0.01, relative = TRUE, runs = 100){
  longlat = longlat[!duplicated(longlat),]                #first, remove duplicate rows
  nSites = nrow(longlat)
  if(nSites < 4)
    return(longlat)
  
  ##if relative, calculate maxDist between any two points
  if(relative){
    if(nSites < 40){ #if limited number of sites use all data
      maxDist = 0
      for(x in 1:(nSites-1)){
        for(y in (x+1):nSites){
          maxDist = max(maxDist,((longlat[x,1]-longlat[y,1])^2+(longlat[x,2]-longlat[y,2])^2)^.5)
        }
      }
    } else { #if many sites use hypothenusa of square encompassing all of them
      horiDist = max(longlat[,1]) - min(longlat[,1])
      vertDist = max(longlat[,2]) - min(longlat[,2])
      maxDist = (horiDist^2 + vertDist^2)^0.5
    }
    distance = maxDist*distance
  }
  
  listSites = matrix(longlat[1,], ncol=2, byrow = TRUE)
  for (r in 1:runs){
    longlat = longlat[sample(nSites),]       ##shuffle rows (sites)
    rndSites = longlat[1,]                   ##start with first random site
    for(newSite in 2:nSites){
      for(oldSite in 1:(newSite-1)){
        addSite = TRUE
        dist = ((longlat[newSite,1]-longlat[oldSite,1])^2+(longlat[newSite,2]-longlat[oldSite,2])^2)^.5
        if(dist < distance){
          addSite = FALSE
          break
        }
      }
      if(addSite)
        rndSites = rbind(rndSites, longlat[newSite,])
    }
    if(nrow(rndSites) > nrow(listSites))
      listSites = rndSites
  }
  return(as.matrix(listSites))
}


#' Move records to closest non-NA cell.
#' @description Identifies and moves presence records to cells with environmental values.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster.
#' @param buffer Maximum distance in map units that a record will move. If 0 all NA records will be changed.
#' @details Often records are in coastal or other areas for which no environmental data is available. This function moves such records to the closest cells with data so that no information is lost during modelling.
#' @return A matrix with new coordinate values.
#' @examples rast <- terra::rast(matrix(c(rep(NA,100), rep(1,100), rep(NA,100)), ncol = 15))
#' pts <- cbind(runif(100, 0, 0.55), runif(100, 0, 1))
#' terra::plot(rast)
#' points(pts)
#' pts <- move(pts, rast)
#' terra::plot(rast)
#' points(pts)
#' @export
move <- function(longlat, layers, buffer = 0){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
  
  if(is(longlat, "matrix")){
    longlat = as.data.frame(longlat)
  }
  
  if(terra::crs(layers) == ""){
    terra::crs(layers) = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  }
  
  # layers <- layers[[1]]
  values <- terra::extract(layers, longlat)   #get values of each record
  suppressWarnings(
    for(i in which(is.na(values))){    #if a value is NA, move it
      # Distance does not work when the crs is "". needs an exception
      distRaster <- terra::distance(
        layers,
        terra::vect(as.data.frame(longlat)[i, ], # remove if enforced at start
          geom = colnames(longlat),
          crs = terra::crs(layers)
        )
      )
      distRaster <- terra::mask(distRaster, layers)
      # vmin = distRaster@ptr$range_min
      vmin <- terra::where.min(distRaster)
      
      if(buffer <= 0 || buffer > vmin){
        # vmin = terra::as.points(distRaster, function(x) x == vmin)
        longlat[i,] = terra::xyFromCell(distRaster, vmin[2]) # vmin[1,1:2]
      }
      
      
    }
  )
  return(longlat)
}


#' Visual detection of outliers.
#' @description Draws plots of sites in geographical (longlat) and environmental (2-axis PCA) space.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers Raster* object as defined by package raster. It can be any set of environmental layers thought to allow the identification of environmental outliers.
#' @details Erroneous data sources or errors in transcriptions may introduce outliers that can be easily detected by looking at simple graphs of geographical or environmental space.
#' @return A data.frame with coordinate values and distance to centroid in pca is returned. Two plots are drawn for visual inspection. The environmental plot includes row numbers for easy identification of possible outliers.
#' @examples records = gecko.examples("gecko.records")
#' layers = gecko.examples("gecko.layers")
#' outliers(records, layers[[1:3]])
#' @export
outliers <- function(longlat, layers){
  userpar <- par(no.readonly = TRUE) 
  on.exit(par(userpar))
  if(dim(layers)[3] == 33){
    pca <- reduce(layers[[1:19]], n = 2) #if layers come from read
  } else {
    pca <- reduce(layers, n = 2)
  }
    
  
  ##extract pca values from longlat
  pca <- as.data.frame(terra::extract(pca, longlat))
  goodRows <-  which(!is.na(pca[,1]))
  pca <- pca[goodRows,]
  longlat <- longlat[goodRows,]
  par(mfrow = c(1,2))
  map.draw(longlat, layers[[1]], spName = "Geographical")
  # plot(pca, main = "Environmental", type = "n")
  plot(pca, main = "Environmental")
  centroid = colMeans(pca)
  text(centroid[1], centroid[2], label = "X")
  for(i in 1:nrow(pca)){
    text(pca[i,1], pca[i,2], label = row.names(longlat)[i])
  }
  
  ##build new matrix ordered by distance to centroid
  dist2centroid = apply(pca, 1, function(x) dist(rbind(x, centroid)))
  out = as.data.frame(cbind(longlat, dist2centroid))
  out = out[order(-dist2centroid),]
  return(out)
}


#' Reduce dimensionality of raster layers.
#' @description Reduce the number of layers by either performing a PCA on them or by eliminating highly correlated ones.
#' @param layers Raster* object as defined by package raster.
#' @param method Either Principal Components Analysis ("pca", default) or Pearson's correlation ("cor").
#' @param n Number of layers to reduce to.
#' @param thres Value for pairwise Pearson's correlation above which one of the layers (randomly selected) is eliminated.
#' @details Using a large number of explanatory variables in models with few records may lead to overfitting. This function allows to avoid it as much as possible.
#' If both n and thres are given, n has priority. If method is not recognized and layers come from read function, only landcover is reduced by using only the dominating landuse of each cell.
#' @return A RasterStack object.
#' @export
reduce <- function(layers, method = "pca", n = NULL, thres = NULL){
  ##method = "pca, cor", if unrecognized method only reduce landcover but not climate

  if(dim(layers)[3] == 33){          ##check if layers are obtained with read
    out <- c(layers[[33]])
    layers = layers[[1:19]]
  }
  # HANDLE COR
  if(method == "cor"){                       ##if correlation
    if(is.null(n)){
      if(is.null(thres))
        thres = 0.7
      for(i in 1:dim(layers)[3]){                  ##delete layers until none are correlated above threshold
        # cor = as.matrix(as.dist(layerStats(layers, 'pearson', na.rm = TRUE)[[1]]))
        cor = as.matrix(
          as.dist(terra::layerCor(layers, 'pearson', na.rm = TRUE)[[1]])
          )
        
        if(max(cor) < thres)
          break
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    } else {
      while (dim(layers)[3] > n){                   ##delete layers until reaching n layers
        cor = abs(as.matrix(as.dist(
          terra::layerCor(layers, 'pearson', na.rm = TRUE)[[1]])))
        
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    }
  } else if(method == "pca"){                                  ##if pca
    if(is.null(n))
      n = 3
    if(sum(!is.na(terra::values(layers[[1]], mat = FALSE))) > 2000)
      sr <- terra::spatSample(layers, 1000)
    else
      sr <- terra::spatSample(layers, as.integer(sum(!is.na(terra::values(layers[[1]], mat = FALSE)))/2), na.rm = TRUE) # added na.rm
    pca <- prcomp(sr)
    layers <- terra::predict(layers, pca, index = 1:n)
    for(i in 1:n){
      names(layers[[i]]) <- paste("pca",i)
    }
  }
  
  if(dim(layers)[3] == 33){
    out <- c(layers, out)
  } else {
    out <- layers
  }
  return(out)
}