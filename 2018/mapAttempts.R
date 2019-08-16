library(sp)
library(lattice) # required for trellis.par.set():
trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()

# prepare nc sids data set:
library(maptools)
nc <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], proj4string=CRS("+proj=longlat +datum=NAD27"))
arrow = list("SpatialPolygonsRescale", layout.north.arrow(),
             offset = c(-76,34), scale = 0.5, which = 2)
#scale = list("SpatialPolygonsRescale", layout.scale.bar(),
#    offset = c(-77.5,34), scale = 1, fill=c("transparent","black"), which = 2)
#text1 = list("sp.text", c(-77.5,34.15), "0", which = 2)
#text2 = list("sp.text", c(-76.5,34.15), "1 degree", which = 2)
## multi-panel plot with filled polygons: North Carolina SIDS
spplot(nc, c("SID74", "SID79"), names.attr = c("1974","1979"),
       colorkey=list(space="bottom"), scales = list(draw = TRUE),
       main = "SIDS (sudden infant death syndrome) in North Carolina",
       sp.layout = list(arrow), as.table = TRUE)

#    sp.layout = list(arrow, scale, text1, text2), as.table = TRUE)



spplot(VAbasins,c('BASIN'),names.attr = c("Major Basins"),
       colorkey=list(space="bottom"), scales = list(draw = TRUE),
       main = "Virginia Probmon sites",
       sp.layout = list(arrow), as.table = TRUE)
       
spplot(VAbasins,'BASIN',col.regions = brewer.pal(n = length(unique(VAbasins@data$BASIN)), name = "Paired"), sp.layout = sitesLayer)

sitesLayer <- list("sp.points", IRstations2018, col = "grey21", pch=19,cex=0.5)

plot(IRstations2018,col='grey21',pch=19,cex=0.5,add=T)
legend('left',legend=c('Major Basins','Probabilistic Sites (Wadeable)')
       ,title='Legend',bty='n',inset=0.05,pch=c(15,15,19),cex=0.8
       ,col=c('wheat1','grey21'))
