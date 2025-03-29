# Hillshade map
# www.overfitting.net
# https://www.overfitting.net/

library(terra)  # resample function
library(tiff)  # save 16-bit TIFF's


# Hillshade calculation
hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction. It can be defined in three ways:
    #   a) 1D value indicating light source Azimuth in degrees (0-360)
    #      (0=North, 90=East, 180=South, 270=West)
    #   b) 2D vector indicating light source (X,Y) coordinates
    #   c) 3D vector indicating light source (X,Y,Z) coordinates:
    #      (X=South, Y=East, Z=Up)
    #      dlight=c(0, 2, 3)  # sunrise
    #      dlight=c(0, 0, 1)  # midday
    #      dlight=c(0,-2, 3)  # sunset
    #   NOTE: both in a) and b) a 45º Elevation angle is applied
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    
    # Deal with lighting direction
    if (length(dlight)==1) dlight=c(-cos(dlight*pi/180), sin(dlight*pi/180))
    if (length(dlight)==2) dlight=c(dlight, (dlight[1]^2+dlight[2]^2)^0.5)
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}

# Generic array resample function
# works both for matrix (grayscale images) or 3-channel arrays (colour images)
arrayresample=function(img, DIMX, DIMY, method='bilinear') {
    require(terra)
    
    raster=rast(img)
    rasterrs=rast(nrows=DIMY, ncols=DIMX, extent=ext(raster))
    rasterrs=resample(raster, rasterrs, method=method)
    
    if (is.matrix(img)) return (matrix(as.array(rasterrs), nrow=nrow(rasterrs)))
    else return (as.array(rasterrs))  # convert back to matrix/array
}


#################################################

# 1. READ DEM

# Put all 7 images together manually in Photoshop -> "tenerifecomposite.tif"
RESOLUTION=25
MAXIMO=3710.062
DEM=readTIFF("tenerifecomposite.tif")*MAXIMO
hist(DEM[DEM>0], breaks=800)


#################################################

# 2. RESCALE DEM TO FULL HD (1920 X 1080)

f=1
fscale=0.347758887171561*f  # downsampling scaling (chosen to fit in Full HD)
RESOLUTION=RESOLUTION/fscale  # downsampling also reduces RESOLUTION
DIMY=round(nrow(DEM)*fscale)
DIMX=round(ncol(DEM)*fscale)
DEMrs=arrayresample(DEM, DIMX, DIMY)

DIMYfullHD=1080*f
DIMXfullHD=1920*f
DEM=matrix(0, nrow=DIMYfullHD, ncol=DIMXfullHD)
DEM[(DIMYfullHD/2-DIMY/2):(DIMYfullHD/2+DIMY/2-1),
    (DIMXfullHD/2-DIMX/2):(DIMXfullHD/2+DIMX/2-1)]=DEMrs
rm(DEMrs)

writeTIFF(DEM/max(DEM), "tenerifecomposite_fullHD.tif",
          bits.per.sample=16, compression="LZW")


#################################################

# 3. GENERATE BASIC HILLSHADE

# 3 ways to define the same light source
hillshade=hillshademap(DEM, dx=RESOLUTION, dlight=270)  # Azimuth
hillshade=hillshademap(DEM, dx=RESOLUTION, dlight=c(0,-1))  # (X,Y) coords
hillshade=hillshademap(DEM, dx=RESOLUTION, dlight=c(0,-1,1))  # (X,Y,Z) coords
writeTIFF(hillshade, "hillshade.tif", bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hillshade[nrow(hillshade):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=1)),
      asp=nrow(hillshade)/ncol(hillshade), axes=FALSE)


#################################################

# 4. GENERATE FALSE COLOUR HILLSHADE (IDEA FROM JOHN NELSON)

hillshadeW =hillshademap(DEM, dx=RESOLUTION, dlight=270)
hillshadeNW=hillshademap(DEM, dx=RESOLUTION, dlight=-45)
hillshadeN =hillshademap(DEM, dx=RESOLUTION, dlight=0)
img=replicate(3, hillshadeW)
img[,,2]=hillshadeNW
img[,,3]=hillshadeN
writeTIFF(img, "hillshadecolour.tif", bits.per.sample=16, compression="LZW")


