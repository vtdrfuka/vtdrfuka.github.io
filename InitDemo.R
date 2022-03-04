if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB,classInt,SWATmodel)
dir.create("~/MINT_SWATdemo")
setwd("~/MINT_SWATdemo/")

#### Global Elevation Retrieval. 
url="https://maps.google.com"
browseURL(url)

parseloc="https://www.google.com/maps/@11.7558581,37.710779,67574m/"
coords=as.numeric(strsplit(sub("https://www.google.com/maps/@(.*)m.*","\\1",parseloc),",")[[1]])
zoom=round(log(35200000/coords[3])/log(2))+1

proj4_utm = paste0("+proj=utm +zone=", trunc((180+coords[2])/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
latlon <- cbind(coords[2],coords[1])
cbasin_ll <- SpatialPoints(latlon)


proj4string(cbasin_ll)=proj4_ll
cbasin_utm=spTransform(cbasin_ll,crs_utm)
url=paste0("https://www.google.com/maps/@",
           coords[1],",",coords[2],",",zoom,"z")
browseURL(url)
searchlength=.75*coords[3]
bboxpts=cbasin_utm@coords+c(-searchlength,-searchlength)
bboxpts=rbind(bboxpts,cbasin_utm@coords+c(searchlength,searchlength))
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# grabbing and projecting latlon DEM to UTM
mydem_utm=get_aws_terrain(locations=bboxpts@coords, 
                      z = zoom, prj = proj4_utm,src ="aws")

writeRaster(mydem_utm,filename = "mydem.tif",overwrite=T)
plot(mydem_utm)
plot(bboxpts,add=T)
# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
slp=raster("mydemslp.tif")

# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
# Threshold 1/overguesstimate
# mythresh=.01*(raster::cellStats(sca,stat = max))
#
threshold=1000
system(paste0("mpiexec -n 8 threshold -ssa mydemsca.tif -src mydemsrc.tif -thresh ",threshold))
src=raster("mydemsrc.tif")
plot(src)
plot(cbasin_utm,add=T)
zoomext=cbasin_utm@coords
zoomext=rbind(zoomext,zoomext+res(src)*100)
zoomext=rbind(zoomext,zoomext-res(src)*100)
zoomext=SpatialPoints(zoomext,proj4string = crs_utm)
zoom(src,ext=zoomext)
plot(cbasin_utm,add=T)
# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)

## a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1){
  xy=locator(n=n)
  points(xy)
#Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

zoom(src)
makeshape.r("approxoutlets")


#outlet=SpatialPointsDataFrame(myflowgage$gagepoint_utm,data.frame(Id=c(1),outlet=paste("outlet",1,sep="")))
#writeOGR(outlet,dsn=".",layer="approxoutlets",driver="ESRI Shapefile", overwrite_layer=TRUE)

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o approxoutlets.shp -om outlet.shp")
outpt=readOGR("outlet.shp")
approxpt=readOGR("approxoutlets.shp")
zoom(src,ext=zoomext)
plot(approxpt, add=T,col="red")
plot(outpt, add=T,col="blue")

# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 
#threshold=myflowgage$area*1000*1000/(res(mybasindem)[1]^2)/10
# Threshold
system(paste0("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh ",threshold))
src1=raster("mydemsrc1.tif")
plot(src1)
plot(outpt, add=T,col="blue")
# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemssa.tif -src mydemsrc1.tif -o outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")

# Build a mask, trim, and crop the basin files
mydemw=raster("mydemw.tif")
plot(mydemw)
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem_utm,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
slp=raster("mydemslp.tif")
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )

nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
plot(TIC)
# A series of plots to show all of the components
#
# DEM - Filled DEM
mybasinfel=crop(fel,mybasinmask)
mybasinfel=mask(mybasinfel,mybasinmask)

par(mfrow = c(2, 2))
plot(mybasinfel)
plot(mybasinslp)
plot(TI)
plot(TIC)
streams=readOGR("mydemnet.shp")
plot(streams,add=TRUE)

# Lets look closer at our "subbasins"!
mybasindemw=crop(mydemw,mybasinmask)
mybasindemw=mask(mybasindemw,mybasinmask)
plot(mybasindemw)



### Vector Soils Sections
url="http://www.fao.org/geonetwork/srv/en/resources.get?id=14116&fname=DSMW.zip&access=private"
download.file(url,"DSMW.zip")
unzip("DSMW.zip")
DSMW_ll=readOGR("DSMW.shp")
cbasin_ll_null=cbasin_ll
bbox_ll=spTransform(bboxpts,crs_ll)
proj4string(bbox_ll)=""

DSMW_ll_crop=crop(DSMW_ll,bbox_ll)
proj4string(DSMW_ll_crop)=crs_ll
plot(DSMW_ll_crop,add=TRUE)
  
DSMW_utm=spTransform(DSMW_ll_crop,crs_utm)
raster::plot(mydemw)
plot(DSMW_utm,add=T)
rmysoil_utm=rasterize(mysoil_utm,TIC,field=as.numeric(mysoil_utm$mukey))

# GHCN vs CFSR Weather data for Basin


