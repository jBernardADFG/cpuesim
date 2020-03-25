# # Read in and prepare sampling area data
# library(rgdal)
# library(raster)
# library(rgeos)
# library(sp)
# setwd("C:/Users/jwbernard/Documents")
# lakelines <- readOGR("~/ArcGIS/Shapefiles/Lake Louise.shp")
# holes <- readOGR("~/ArcGIS/Shapefiles/Lake Louise Holes.shp")
# # Data structure manipulation cause I suck at GIS
# ll_coords <- lakelines@polygons[[1]]@Polygons[[1]]@coords
# h1_coords <- holes@polygons[[1]]@Polygons[[1]]@coords
# h2_coords <- holes@polygons[[2]]@Polygons[[1]]@coords
# ll_p <- Polygon(ll_coords)
# h1_p <- Polygon(h1_coords, hole=T)
# h2_p <- Polygon(h2_coords, hole=T)
# ps <- Polygons(list(ll_p, h1_p, h2_p), ID="sampArea")
# sp <- SpatialPolygons(list(ps))
# crs(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# plot(sp)
# 
# h1_coords <- locator()
# h2_coords <- locator()
# h3_coords <- locator()
# h1_p <- Polygon(h1_coords, hole=T)
# h2_p <- Polygon(h2_coords, hole=T)
# h3_p <- Polygon(h3_coords, hole=T)
# ps <- Polygons(list(ll_p, h1_p, h2_p, h3_p), ID="sampArea")
# sp <- SpatialPolygons(list(ps))
# crs(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# plot(sp)
# save(sp, file="C:/Jordy/louiseOP2020/Data.R")
# 
# help(save)
# 
# # To get coordinate locations for the 2005 data
# png <-raster("S:/Jordy/louiseOP2020/Data/map.png")
# e <- extent(95.43794,1562.79625,278.7775,2256.1850)
# png <- crop(png,e)
# extent(png) <- c(-146.6356, -146.4595, 62.27376, 62.37895)
# crs(png) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# # Read in the 2005 data
# library(readxl)
# dat <- as.data.frame(read_excel("S:/Jordy/louiseOP2020/Data/CPUE.xls", sheet = "CPUEL", col_names = F))
# colnames(dat) <- c("tran", "set", "ct")
# dat[,1][dat[,1]==2] <- 20
# dat <- dat[order(dat$tran, dat$set),]
# table(dat$tran)
# plot(png)
# # Use locator to determine coordinate locations
# load("S:/Jordy/LouiseOP2020/Data/cpue_coord.RData")
# plot(ct_dt[,1:2], col=ct_dt[,3])
# extent(png)
# ref_pt <- c(-146.6356 , 62.27376)
# new_x <- new_y <- rep(NA,nrow(dat))
# for(i in 1:nrow(dat)){
#   new_x[i] <- get.geo.dist(ref_pt[1],ref_pt[2], ct_dt[i,1], ref_pt[2])
#   new_y[i] <- get.geo.dist(ref_pt[1],ref_pt[2], ref_pt[1], ct_dt[i,2])
# }
# ct_dt <- data.frame(new_x, new_y, ct_dt$count)
# names(ct_dt) <- c("x", "y", "count")
# coordinates(ct_dt) <- c("x", "y")
# fit.variog(ct_dt)
# get.geo.dist(long1, lat1, long2, lat2, units = "km")
# get.geo.dist(-146.6356, 62.27376, -146.4595, 62.38155)
# fit.variog(ct_dt)
# coordinates(ct_dt) <- c("x", "y")
# 
# 
# 
# 
# 
# 
# # Get the empirical variogram
# tran <- unique(dat[,1])
# df <- data.frame()
# for(i_1 in 1:nrow(dat)){
#   for(i_2 in 1:nrow(dat)){
#     if (dat[i_1,1]==dat[i_2,1]){ # sets are in the same transect
#       d <- abs(dat[i_1,2]-dat[i_2,2])
#       d <- d*125
#       v <- (dat[i_1,3]-dat[i_2,3])^2
#       df <- rbind(df, c(d,v))
#     }
#   }
# }
# df <- df[order(df[,1]),]
# df <- df[!df[,1]==0,]
# rdf <- data.frame()
# for(d in unique(df[,1])){
#   ddf <- df[df[,1]==d,]
#   sv <- sum(ddf[,2])/(nrow(ddf)) # Empirical semivariance
#   n <- nrow(ddf)
#   rdf <- rbind(rdf, c(d, sv, n))
# }
# rdf <- data.frame(rdf)
# rdf <- rdf[rdf[,1]<3000,]
# with(rdf, symbols(x=rdf[,1], y=rdf[,2], circles=rdf[,3], inches=1/5,
#                   ann=F, bg="steelblue2", fg=NULL))
# 
# 
# 
# 
# 
# 
# 
# 
# library(sp)
# data(meuse)
# # no trend:
# coordinates(meuse) = ~x+y
# v<-variogram(log(zinc)~1, meuse)
# 
# 
# 
# 
# 
# 
# fit.variogram(vgm, vgm(1, "Sph", 300, 1))
# fit.variogram(vgm, vgm("Exp"))
# 
# 
# #
# # computing variograms:
# #
# # binned variogram
# vario.b <- variog(s100, max.dist=1)
# # variogram cloud
# vario.c <- variog(s100, max.dist=1, op="cloud")
# #binned variogram and stores the cloud
# vario.bc <- variog(s100, max.dist=1, bin.cloud=TRUE)
# # smoothed variogram
# vario.s <- variog(s100, max.dist=1, op="sm", band=0.2)
# #
# #
# # plotting the variograms:
# par(mfrow=c(2,2))
# plot(vario.b, main="binned variogram")
# plot(vario.c, main="variogram cloud")
# plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")
# plot(vario.s, main="smoothed variogram")
# 
# # computing a directional variogram
# vario.0 <- variog(s100, max.dist=1, dir=0, tol=pi/8)
# plot(vario.b, type="l", lty=2)
# lines(vario.0)
# legend("topleft", legend=c("omnidirectional", expression(0 * degree)), lty=c(2,1))
# 
# 
# 
# locator()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# dfx = data.frame(ev1=1:10, ev2=sample(10:99, 10), ev3=10:1)
# 
# with(dfx, symbols(x=ev1, y=ev2, circles=ev3, inches=1/3,
#                   ann=F, bg="steelblue2", fg=NULL))
# 
# 
# 
# 
# 
# 
# 
# ev_3
# 
# 
# dat_coords_df <- data.frame(dat_coords$x, dat_coords$y)
# count <- c(2,2,2,3,3,2,2,3,3,3,2,2,2,2,1)
# data <- data.frame(dat_coords,count)
# coordinates(data) <- c("x","y")
# vario_fit<-fit.variog(data)
# 
# 
# 
# 
# help("Emp.variog")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# data(slp)
# day <- slp$date.obs
# id <- slp$id.stat
# coord1 <- slp$lon.stat
# coord2 <- slp$lat.stat
# obs <- slp$obs
# forecast <- slp$forecast
# 
# ## Computing variogram
# ## No specified cutpoints, no specified maximum distance
# ## Default number of bins
# variogram <- Emp.variog(day=day,obs=obs,forecast=forecast,id=id,
#                         coord1=coord1,coord2=coord2,cut.points=NULL,max.dist=NULL,nbins=NULL)
# ## Plotting variogram
# plot(variogram$bin.midpoints,variogram$empir.variog,xlab="Distance",
#      ylab="Semi-variance",main="Empirical variogram")
# 
# 
# 
# Variog.fit(emp.variog, variog.model="exponential", max.dist.fit=NULL,
#            init.val=NULL, fix.nugget=FALSE)
# 
# 
# 
# variogram
# 
# 
# 
# data(aniso)
# va <- svariog(X=aniso, plot=FALSE)
# fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,10), nugget=0.2, max.dist=30, plot = TRUE)
# fit
# 
# data(sim03)
# va <- svariog(X=sim03, plot=TRUE)
# fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5, max.dist=200)
# fit
# 
# # graphical display
# fit <- fitsvariog(vario=va, ini.cov.pars=c(0.05,4.5), nugget=0.5, max.dist=200, plot=FALSE)
# plot(va$svario$u, va$svario$v)
# lines(fit$fit)
# 
# ## fit model to empirical variograms from field data and see how the maximum distance
# ## to be used can change the results
# data(crypho)
# va <- svariog(X=crypho, plot=TRUE, messages=FALSE)
# fit1 <- fitsvariog(vario=va, ini.cov.pars=c(0.03,100), nugget=0.1, max.dist=300, plot = TRUE)
# fit2 <- fitsvariog(vario=va, ini.cov.pars=c(0.03,100), nugget=0.1, max.dist=600, plot = TRUE)
# 
# # plot results
# plot(va$svario$u, va$svario$v)
# lines(fit1$fit, col="blue")
# lines(fit2$fit, col="red")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
