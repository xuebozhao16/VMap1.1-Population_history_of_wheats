##https://zhuanlan.zhihu.com/p/38164684 这是学习的代码
rm(list=ls())
library(assertthat)
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)
library(ggraph)
library(ggmap)
country_coords_txt <- "
1     3.00000  28.00000       Algeria
2    54.00000  24.00000           UAE
3   139.75309  35.68536         Japan
4    45.00000  25.00000 'Saudi Arabia'
5     9.00000  34.00000       Tunisia
6     5.75000  52.50000   Netherlands
7   103.80000   1.36667     Singapore
8   124.10000  -8.36667         Korea
9    -2.69531  54.75844            UK
10    34.91155  39.05901        Turkey
11  -113.64258  60.10867        Canada
12    77.00000  20.00000         India
13    25.00000  46.00000       Romania
14   135.00000 -25.00000     Australia
15    10.00000  62.00000        Norway"
# nodes come from the above table and contain geo-coordinates for some
# randomly picked countries
nodes <- read.delim(text = country_coords_txt, header = FALSE,
                    quote = "'", sep = "",
                    col.names = c('id', 'lon', 'lat', 'name'))
set.seed(123)  # set random generator state for the same output

N_EDGES_PER_NODE_MIN <- 1
N_EDGES_PER_NODE_MAX <- 4
N_CATEGORIES <- 4
# edges: create random connections between countries (nodes)
edges <- map_dfr(nodes$id, function(id) {
  n <- floor(runif(1, N_EDGES_PER_NODE_MIN, N_EDGES_PER_NODE_MAX+1))
  to <- sample(1:max(nodes$id), n, replace = FALSE)
  to <- to[to != id]
  categories <- sample(1:N_CATEGORIES, length(to), replace = TRUE)
  weights <- runif(length(to))
  data_frame(from = id, to = to, weight = weights, category = categories)
})

edges <- edges %>% mutate(category = as.factor(category))
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
edges_for_plot <- edges %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('from' = 'id')) %>%
  rename(x = lon, y = lat) %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('to' = 'id')) %>%
  rename(xend = lon, yend = lat)
assert_that(nrow(edges_for_plot) == nrow(edges))
nodes$weight = degree(g)
maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))
country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
                               data = map_data('world'),
                               fill = "#CECECE", color = "#515151",
                               size = 0.15)
mapcoords <- coord_fixed(xlim = c(-150, 180), ylim = c(-55, 80))
###
ggplot(nodes) + country_shapes +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     # draw edges as arcs
                 color = category, size = weight),
             data = edges_for_plot, curvature = 0.33,
             alpha = 0.5) +
  scale_size_continuous(guide = FALSE, range = c(0.25, 2)) + # scale for edge widths
  geom_point(aes(x = lon, y = lat, size = weight),           # draw nodes
             shape = 21, fill = 'white',
             color = 'black', stroke = 0.5) +
  scale_size_continuous(guide = FALSE, range = c(1, 6)) +    # scale for node size
  geom_text(aes(x = lon, y = lat, label = name),             # draw text labels
            hjust = 0, nudge_x = 1, nudge_y = 4,
            size = 3, color = "white", fontface = "bold") +
  mapcoords + maptheme

########################以上是测试的代码，现在是自己的数据啦
nodes <- read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/map/mig_map1.txt",
                                 header = F,sep = "\t",col.names = c('id', 'lon', 'lat', 'name'))

N_EDGES_PER_NODE_MIN <- 1
N_EDGES_PER_NODE_MAX <- 4
N_CATEGORIES <- 4
# edges: create random connections between countries (nodes)
edges <- map_dfr(nodes$id, function(id) {
  n <- floor(runif(1, N_EDGES_PER_NODE_MIN, N_EDGES_PER_NODE_MAX+1))
  to <- sample(1:max(nodes$id), n, replace = FALSE)
  to <- to[to != id]
  categories <- sample(1:N_CATEGORIES, length(to), replace = TRUE)
  weights <- runif(length(to))
  data_frame(from = id, to = to, weight = weights, category = categories)
})

edges <- edges %>% mutate(category = as.factor(category))
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
edges_for_plot <- edges %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('from' = 'id')) %>%
  rename(x = lon, y = lat) %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('to' = 'id')) %>%
  rename(xend = lon, yend = lat)
assert_that(nrow(edges_for_plot) == nrow(edges))
nodes$weight = degree(g)
maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))
country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
                               data = map_data('world'),
                               fill = "#CECECE", color = "#515151",
                               size = 0.15)
mapcoords <- coord_fixed(xlim = c(-150, 180), ylim = c(-55, 80))
###
ggplot(nodes) + country_shapes +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     # draw edges as arcs
                 color = category, size = weight),
             data = edges_for_plot, curvature = 0.33,
             alpha = 0.5) +
  scale_size_continuous(guide = FALSE, range = c(0.25, 2)) + # scale for edge widths
  geom_point(aes(x = lon, y = lat, size = weight),           # draw nodes
             shape = 21, fill = 'white',
             color = 'black', stroke = 0.5) +
  scale_size_continuous(guide = FALSE, range = c(1, 6)) +    # scale for node size
  geom_text(aes(x = lon, y = lat, label = name),             # draw text labels
            hjust = 0, nudge_x = 1, nudge_y = 4,
            size = 3, color = "white", fontface = "bold") +
  mapcoords + maptheme

#######
library(OpenStreetMap)
library(ggplot2)
map <- openmap(c(55,-30), c(4,130),type='esri')
plot(map,raster=TRUE)

###########################################################现在使用的是另一种方法做eems
#######现在做的是migration的图---D
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out4")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out_plot4","EEMS_D8")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/D/D_ex_out_plot/out",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=longlat +ellps=sphere +no_defs"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
#map_Eurasia <- map_world[which(map_world@data$continent == "Eurasia"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, col.outline = "gray90",
           add.demes =TRUE,col.demes = "red",min.cex.demes = 0.3, max.cex.demes = 1,add.grid = TRUE,
           m.plot.xy = { plot(map_world, col = NA, add = TRUE) },q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_world <- spTransform(map_world, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_world, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_world, col = NA, add = TRUE) },out.png = FALSE)
dev.off()



##现在看一下选择的闭合回路是不是
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/")
cood = read.table("map_file.txt",head =F)
plot(cood$V1,cood$V2,pch=19)
plot(cood$V1,cood$V2,type="l")



####现在做的是Structu re的rejection
require(reshape)
require (rworldmap)
require(rworldxtra)
col1 = rgb(28/255, 93/255, 127/255,1) #深蓝
col2 = rgb(155/255, 204/255, 227/255,1) #浅蓝
col3 = rgb(232/255, 228/255, 174/255,1) #黄色
col4 = rgb(227/255, 88/255, 87/255,1) #红色
col5 = rgb(233/255, 189/255, 180/255,1) #浅红
wheat1 <- read.table("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/Piemap_data.txt", header=T,sep="\t")
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$k.5),]
wheat_reshape <- wheat_reshape <- cast(wheat,cluster1_70+cluster2_70~k.5) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-15,145),ylim=c(-5,40),nameX="cluster2_70",nameY="cluster1_70",nameZs=c('1','2','3','4','5'),
        symbolSize=1,maxZVal = 1,
        zColours=c(col3,col2,col1,col5,col4),barOrient='vert',borderCol = NA,oceanCol="#d1d2d4",landCol="white")
dev.off()















