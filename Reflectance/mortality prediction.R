# Prediction of Vineyard mortality

install.packages("shapefile")
library(rgdal) # for vector data processing
library(sp) # for spatial data processing
library(corrplot) # helps to many kind of plots
library(dplyr) # To manipulate dataframes
library(raster) # for raster data processing

library(caret) # for sorting data sets
library(randomForest) # for random forest classifiers functions
library(e1071) # for  svm classifier functions
library(factoextra) 
library(Factoshiny) # for PCA functions

# Create a CRS object
proj.4326=crs("+init=epsg:4326") # create crs object 4326

# Open vector datas 
rpg=readOGR("src","RPG_2018_VRC")
vm.rpg=readOGR("src","VM_RPG_2018")

# Reproject vector datas
rpg.wgs84=spTransform(rpg,proj.4326)
vm.rpg.wgs84=spTransform(vm.rpg,proj.4326)

plot(rpg.wgs84)
plot(vm.rpg.wgs84)
rpg.wgs84@polygons
## Request the type of data
class(rpg)
length(rpg)
extent(rpg)
crs(rpg)

# reading vector attribute tables
tbl.rpg.wgs84=rpg.wgs84@data
tbl.vm.rpg.wgs84=vm.rpg.wgs84@data


# Open raster and extent vector datas | raster
blue=raster("src/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9_FRE_B2.tif")
green=raster("src/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9_FRE_B3.tif")
red=raster("src/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9_FRE_B4.tif")
nir=raster("src/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9/SENTINEL2B_20180920-104712-964_L2A_T31TEH_D_V1-9_FRE_B8.tif")
stack.band=stack(blue,green,red,nir)

study.area=readOGR("src/mygeodata","Study_area_JRL_remote_sensing-polygon")

## Change CRS and crop raster to the extent of the study area | projectRaster and Crop
bands.wgs84=projectRaster(stack.band,crs=proj.4326)
bands.study.erea=crop(bands.wgs84,study.area)

blue.band=bands.study.erea$SENTINEL2B_20180920.104712.964_L2A_T31TEH_D_V1.9_FRE_B2
green.band=bands.study.erea$SENTINEL2B_20180920.104712.964_L2A_T31TEH_D_V1.9_FRE_B3
red.band=bands.study.erea$SENTINEL2B_20180920.104712.964_L2A_T31TEH_D_V1.9_FRE_B4
nir.band=bands.study.erea$SENTINEL2B_20180920.104712.964_L2A_T31TEH_D_V1.9_FRE_B8

## Generate indexes 
L=0.5
NDVI= (nir.band-red.band)/(nir.band + red.band)
SAVI= (1+L)*((nir.band - red.band)/(L+nir.band+red.band))

par(mfrow=c(1,2))
plot(NDVI)
plot(SAVI)

stack.indices.bands=stack(bands.study.erea,NDVI,SAVI)

# Extract NDVI, SAVI and reflectance from vector datas | extract
# VM.RPG
extract.mean.1=raster::extract(stack.indices.bands,vm.rpg.wgs84, fun=mean)
extract.sd.1=raster::extract(stack.indices.bands,vm.rpg.wgs84, fun=sd)
indices.vm.rpg=cbind(extract.mean.1,extract.sd.1)
join.indices.vm.rpg=cbind(vm.rpg.wgs84[,6],indices.vm.rpg)
#join.indices.vm.rpg=join.indices.vm.rpg %>% rename(mean_ndvi=X1,mean_savi=X2,sd_ndvi=X3,sd_savi=X4)
colnames(join.indices.vm.rpg@data)[2:13]=c("mean_b2","mean_b3","mean_b4","mean_b5",
                                          "mean_ndvi","mean_savi","sd_b2","sd_b3",
                                          "sd_b4","sd_b8", "sd_ndvi","sd_savi")

# RPG
extract.mean.2=raster::extract(stack.indices.bands,rpg.wgs84, fun=mean)
extract.sd.2=raster::extract(stack.indices.bands,rpg.wgs84, fun=sd)
indices.rpg=cbind(extract.mean.2,extract.sd.2)
join.indices.rpg=cbind(rpg.wgs84[,8],indices.rpg)
#join.indices.rpg=join.indices.rpg %>% rename(mean_ndvi=X1,mean_savi=X2,sd_ndvi=X3,sd_savi=X4)
colnames(join.indices.rpg@data)[2:13]=c("mean_b2","mean_b3","mean_b4","mean_b8",
                                        "mean_ndvi","mean_savi","sd_b2","sd_b3",
                                        "sd_b4","sd_b8", "sd_ndvi","sd_savi")
#join.indices.rpg=cbind(rpg.wgs84,indices.rpg)
#join.indices.rpg=join.indices.rpg@data %>% rename(mean_ndvi=X1,mean_savi=X2,sd_ndvi=X3,sd_savi=X4)

# Extracting training and testing data. 
training=subset(join.indices.vm.rpg, COD_INT %in% c(0:3))
#training$COD_INT[training$COD_INT==0]="no_mortality"
#training$COD_INT[training$INT==1]="low_mortality"
#training$COD_INT[training$COD_INT==2]="med_mortality"
#training$COD_INT[training$COD_INT==3]="high_mortality"
training@data[1]=lapply(training@data[1], as.factor)

# Sort train (70%) and test (30%) data from traning data
train.data=na.omit(training@data)

# Train Random forest  classifier
## RF
set.seed(123)
rf=randomForest(COD_INT~., data=train.data) 
varImpPlot(rf)
rf$importance

test_pred = predict(rf, test)
mat_conf=confusionMatrix(test_pred, test$COD_INT)
kappa= t(mat_conf$overall[1:2])
precision.prod_user = mat_conf$byClass[,c(1,3)]

prec_glob=data.frame(kappa)
prec_glob=matrix(prec_glob, ncol = 2, byrow = T)
write.csv2(prec_glob, 'Tables/précisions.csv')

prec_prod_util=as.data.frame(prod_util)
write.csv2(prec_prod_util, 'Tables/confusion.xls')
write.csv2(mat_conf[[17]]$table, 'Tables/mat_conf.csv')

# Trying PCA

pca= PCA(train.data,quali.sup=1,ncp=5)
pca$ind$coord
pca$eig
fviz_eig(pca)

# sort train and test data
new.train.data=na.omit(cbind(train.data[1],pca$ind$coord))
#set.seed(123)
#intrain = createDataPartition(y=new.train.data$COD_PAT, p=0.7, list = FALSE)
#train2= train.data[intrain,]
#test2= train.data[-intrain,]

# Train Random forest  classifier
## RF
set.seed(123)
rf2=randomForest(COD_INT~., data=new.train.data)

## prediction
## VM RPG
vm.rpg.pca= PCA(join.indices.vm.rpg@data,quali.sup=1,ncp=5)
data.vm.rpg.pca=data.frame(vm.rpg.pca$ind$coord)

pred.vm.rpg=predict(rf2,data.vm.rpg.pca)

final.class.1=cbind(join.indices.vm.rpg,pred.vm.rpg)

# RPG 
rpg.pca= PCA(join.indices.rpg@data[-1],quali.sup=1,ncp=5)
data.rpg.pca=data.frame(rpg.pca$ind$coord)

pred.rpg=predict(rf2,data.rpg.pca)

final.class.1=cbind(join.indices.rpg,pred.rpg)

### Export shapefile
writeOGR(final.class.1, dsn = "classified.vm.rpg", layer = "vm.rpg",driver = "ESRI Shapefile" )
writeOGR(final.class.2, dsn = "classified.rpg", layer = "rpg",driver = "ESRI Shapefile" )

