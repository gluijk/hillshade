# Construyendo mapas de sombras hillshade con R
# www.overfitting.net
# https://www.overfitting.net/2022/07/construyendo-mapas-de-sombras-hillshade.html

library(data.table)  # fread()
library(tiff)


# Centro de Descargas del Centro Nacional de Información Geográfica
# Modelos de elevaciones en formato raster MDT25
# Cotas en m, resolución rejilla=25m
# URL: http://centrodedescargas.cnig.es/CentroDescargas/index.jsp

# Leemos y procesamos datos raster
# 4 cuadrantes Sierra Norte de Madrid (Sierra de Guadarrama, Valle del Lozoya)
sierra_11=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0483_LID.txt", sep=" ", dec="."))
sierra_12=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0484_LID.txt", sep=" ", dec="."))
sierra_21=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0508_LID.txt", sep=" ", dec="."))
sierra_22=data.matrix(
    fread("PNOA_MDT25_ETRS89_HU30_0509_LID.txt", sep=" ", dec="."))

# Eliminate overlaps and final crop (values obtained manually)
DEM=matrix(0, nrow=1508, ncol=2269)
DEM[1:759, 11:1148]=sierra_11
DEM[14:768, 1136:2269]=sierra_12
DEM[741:1499, 1:1141]=sierra_21
DEM[754:1508, 1129:2265]=sierra_22
DEM=DEM[14:1499, 11:2265]
rm(sierra_11, sierra_12, sierra_21, sierra_22)

dx=25  # DEM resolution (m)
DIMX=nrow(DEM)
DIMY=ncol(DEM)

# Save DEM
DEMnorm=DEM-min(DEM)
DEMnorm=DEMnorm/max(DEMnorm)
writeTIFF(DEMnorm^(1/1.8),"DEM.tif",
          bits.per.sample=16, compression="LZW")

# Display DEM
image(t(DEM[nrow(DEM):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=1.8)),
      asp=nrow(DEM)/ncol(DEM), axes=FALSE)


# Lighting direction
for (dlight in list(c(0, 2, 3), c(0, 0, 1), c(0, -2, 3))) {
    # dlight=c(0, 2, 3)  # sunrise
    # dlight=c(0, 0, 1)  # midday
    # dlight=c(0,-2, 3)  # sunset
    
    dlightM=sum(dlight^2)^0.5
    
    # Calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMX-2), 2:(DIMY-1)] - DEM[3:DIMX,     2:(DIMY-1)])
    ny = 2*dx*(DEM[2:(DIMX-1), 1:(DIMY-2)] - DEM[2:(DIMX-1), 3:DIMY])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMX-2)x(DIMY-2) matrix
    
    # Reflectance
    hillshade=dn/(dlightM*nM)  # hillshade = cos(theta)
    hillshade[hillshade<0]=0  # clip negative values
    
    # Save hillshade
    writeTIFF(hillshade^(1/0.5),
              paste0("hillshade",dlight[1],dlight[2],dlight[3],".tif"),
              bits.per.sample=16, compression="LZW")
    }

# Display hillshade
image(t(hillshade[nrow(hillshade):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=0.5)),
      asp=nrow(hillshade)/ncol(hillshade), axes=FALSE)
