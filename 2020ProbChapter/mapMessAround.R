# Figure 2.4-2 wadeable sites in IR window

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(VAoutline,border='wheat3',col='wheat1',ylim=c(37.4,38.5))
plot(IRstations2018,col='grey21',pch=19,cex=0.5,add = T)
legend('left',legend=c('Probabilistic Sites (Wadeable)')
       ,title='Legend',bty='n',inset=0.05,pch=c(15,15,19),cex=0.8
       ,col=c('wheat1','grey21'))

# set mapping theme to super basic
#theme_set(theme_bw() +
#            theme(rect = element_blank(),
#                  line = element_blank(), 
#                  panel.grid = element_blank(),
#                  text = element_blank()))

ggplot() +
  geom_sf(data = VAoutline, color = 'wheat3', fill = 'wheat1') +
  geom_sf(data = filter(allProb, IR2018 == 2018), 
          aes(color = 'Probabilistic Sites\n (Wadeable)'),
          size = 1, show.legend = "point") +
  coord_sf(datum = NA) + # take away grid
  scale_colour_manual(values = c('Probabilistic Sites\n (Wadeable)' = 'black'),
                      guide = guide_legend(override.aes = list(linetype = 'blank',
                                                               shape = 19))) +
  theme_minimal() +
  theme(legend.title=element_blank(), # no legend title
        legend.position =  c(0.2,0.52))
  



# plot basins

col.regions = brewer.pal(n = 12,#length(unique(VAbasins@data$BASIN)), 
                         name = "Paired")[c(1:5,7:12)]

# close but dots still in legend boxes

#ggplot() +
#  geom_sf(data = mutate(basinssmooth, Basin = BASIN), 
#          aes(fill = Basin) , show.legend = 'polygon') +
#  geom_sf(data = filter(allProb, IR2018 == 2018), 
#          aes(color = 'Probabilistic Sites\n (Wadeable)'),
#          size = 1, show.legend = 'point') +
#  coord_sf(datum = NA) + # take away grid
#  scale_fill_manual(values=col.regions) +
#  scale_colour_manual(values = c('Probabilistic Sites\n (Wadeable)' = 'black'),
#                      guide = guide_legend(override.aes = list(linetype = "blank", 
#                                                               shape = 19))) +
#  theme_minimal() +
#  theme(legend.title=element_blank(), # no legend title
#        legend.position =  'bottom')#c(0.2,0.52))


ggplot() +
  geom_sf(data = mutate(basinssmooth, Basin = BASIN), 
          aes(fill = Basin) , show.legend = 'polygon') +
  geom_sf(data = filter(allProb, IR2018 == 2018), 
          aes(color = 'Probabilistic Sites\n (Wadeable)'),
          size = 1, show.legend = 'point') +
  coord_sf(datum = NA) + # take away grid
  scale_fill_manual(values= col.regions, name = NULL,
                    guide = guide_legend(override.aes = list(linetype = 'solid',
                                                             shape = NA))) +
  scale_colour_manual(values = c('Probabilistic Sites\n (Wadeable)' = 'black'),
                      guide = guide_legend(override.aes = list(linetype = 0, 
                                                               shape = 19))) +
  theme_minimal() +
  theme(legend.title=element_blank(), # no legend title
        legend.position =  'bottom') + 
  guides(color=guide_legend(override.aes=list(fill=NA, shape = 19))) # still cant get box to go away on dots


plot(st_geometry(basinssmooth), color = 'BASIN')




        )
        #text = element_blank())
    rect = element_blank(),
  #                  line = element_blank(), 
  #                  panel.grid = element_blank(),
  #                  text = element_blank()))
  
    ggplot() +
      geom_sf(data = VAoutline, color = 'wheat3', fill = 'wheat1') +
      geom_sf(data = filter(allProb, IR2018 == 2018), 
              aes(color = 'Wadeable Sites'),
              size = 1, show.legend = 'point')  

nc <- st_read(system.file("shape/nc.shp", package="sf"))      
nc %>% ggplot(aes(fill=AREA)) +                                     
  geom_sf() +                                                         
  theme(legend.position = c(0.5,0))    




### Plot by ecoregions



ggplot() +
  geom_sf(data = mutate(VAecoregionsL3, `Ecoregion (Level III)` = US_L3NAME), 
          aes(fill = `Ecoregion (Level III)`) , show.legend = 'polygon') +
  geom_sf(data = filter(allProb, IR2018 == 2018), 
          aes(color = 'Probabilistic Sites\n (Wadeable)'),
          size = 1, show.legend = 'point') +
  coord_sf(datum = NA) + # take away grid
  scale_fill_manual(values= col.regions, name = NULL,
                    guide = guide_legend(override.aes = list(linetype = 'solid',
                                                             shape = NA))) +
  scale_colour_manual(values = c('Probabilistic Sites\n (Wadeable)' = 'black'),
                      guide = guide_legend(override.aes = list(linetype = 0, 
                                                               shape = 19))) +
  theme_minimal() +
  theme(legend.title=element_blank(), # no legend title
        legend.position =  'bottom') + 
  guides(color=guide_legend(override.aes=list(fill=NA, shape = 19)),
         fill = guide_legend(nrow=3,byrow=TRUE)) # still cant get box to go away on dots




#plot wadeable sites by IR cycle



ggplot() +
  geom_sf(data = VAoutline, color = 'wheat3', fill = 'wheat1') +
  geom_sf(data = filter(allProb, IR2008 == 2008), 
          aes(color = '2008 IR wadeable sites'), size = 1, show.legend = "point") +
  geom_sf(data = filter(allProb, IR2010 == 2010), 
          aes(color = '2010 IR wadeable sites'), size = 1, show.legend = "point") +
  geom_sf(data = filter(allProb, IR2012 == 2012), 
          aes(color = '2012 IR wadeable sites'), size = 1, show.legend = "point") +
  geom_sf(data = filter(allProb, IR2014 == 2014), 
          aes(color = '2014 IR wadeable sites'), size = 1, show.legend = "point") +
  geom_sf(data = filter(allProb, IR2016 == 2016), 
          aes(color = '2016 IR wadeable sites'), size = 1, show.legend = "point") +
  geom_sf(data = filter(allProb, IR2018 == 2018), 
          aes(color = '2018 IR wadeable sites'), size = 1, show.legend = "point") +
  geom_sf(data = filter(allProb, IR2020 == 2020), 
          aes(color = '2020 IR wadeable sites'), size = 1, show.legend = "point") +
  coord_sf(datum = NA) + # take away grid
  scale_colour_manual(values = c('2020 IR wadeable sites' = 'black',
                                 '2018 IR wadeable sites' = 'cyan',
                                 '2016 IR wadeable sites' = 'orange',
                                 '2014 IR wadeable sites' = 'purple',
                                 '2012 IR wadeable sites' = 'green',
                                 '2010 IR wadeable sites' = 'blue',
                                 '2008 IR wadeable sites' = 'red'),
                      guide = guide_legend(override.aes = list(linetype = c('blank','blank','blank','blank','blank','blank','blank'),
                                                               shape = c(19, 19, 19, 19, 19, 19, 19)))) +
  theme_minimal() +
  theme(legend.title=element_blank(), # no legend title
        legend.position =  'bottom') +
  guides(color = guide_legend(nrow=3,byrow=TRUE))







### wadeable and boatable sites
ggplot() +
  geom_sf(data = VAoutline, color = 'wheat3', fill = 'wheat1') +
  geom_sf(data = allProb, aes(color = 'Probabilistic Wadeable\n Sites (2001 - 2018)'),
          size = 1, show.legend = "point") +
  geom_sf(data = boatableSites, aes(color = 'Probabilistic Boatable\n Sites (2008 - 2018)'),
          size = 1, show.legend = "point") +
  scale_colour_manual(values = c('Probabilistic Wadeable\n Sites (2001 - 2018)' = 'grey21',
                                 'Probabilistic Boatable\n Sites (2008 - 2018)' = 'magenta'),
                      guide = guide_legend(override.aes = list(linetype = c('blank','blank'),
                                                               shape = c(19, 19)))) +
  coord_sf(datum = NA) + # take away grid
  theme_minimal() +
  theme(legend.title=element_blank(), # no legend title
        legend.position =  'bottom')











if (requireNamespace("sf", quietly = TRUE)) {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  ggplot(nc) +
    geom_sf(aes(fill = AREA))
  
  # If not supplied, coord_sf() will take the CRS from the first layer
  # and automatically transform all other layers to use that CRS. This
  # ensures that all data will correctly line up
  nc_3857 <- sf::st_transform(nc, 3857)
  ggplot() +
    geom_sf(data = nc) +
    geom_sf(data = nc_3857, colour = "red", fill = NA)
  
  # Unfortunately if you plot other types of feature you'll need to use
  # show.legend to tell ggplot2 what type of legend to use
  nc_3857$mid <- sf::st_centroid(nc_3857$geometry)
  ggplot(nc_3857) +
    geom_sf(colour = "white") +
    geom_sf(aes(geometry = mid, size = AREA), show.legend = "point")
  
  # You can also use layers with x and y aesthetics: these are
  # assumed to already be in the common CRS.
  ggplot(nc) +
    geom_sf() +
    annotate("point", x = -80, y = 35, colour = "red", size = 4)
  
  # Thanks to the power of sf, a geom_sf nicely handles varying projections
  # setting the aspect ratio correctly.
  library(maps)
  world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))
  ggplot() + geom_sf(data = world1)
  
  world2 <- sf::st_transform(
    world1,
    "+proj=laea +y_0=0 +lon_0=155 +lat_0=-90 +ellps=WGS84 +no_defs"
  )
  ggplot() + geom_sf(data = world2)
  
  # To add labels, use geom_sf_label().
  ggplot(nc_3857[1:3, ]) +
    geom_sf(aes(fill = AREA)) +
    geom_sf_label(aes(label = NAME))
}

