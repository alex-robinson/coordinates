


# Load functions
source("~/models/sicopolis/analysis/functions.r")

grid_stats <- function(name,var1,var2,mask2)
{
    err = var2 - var1
    err[mask2==0] = NA 

    n = length(which(!is.na(err)))

    fld_ave   = mean(var2[mask2==1],na.rm=TRUE)
    fld_range = range(var2[mask2==1],na.rm=TRUE)

    MAE       = sum(abs(as.vector(err)),na.rm=TRUE) / n
    AE_SD     = sd(abs(as.vector(err)),na.rm=TRUE)
    RRD       = MAE / diff(fld_range) * 100 

    table = data.frame(name=name, min=fld_range[1],max=fld_range[2], mean=fld_ave, 
                       MAE=MAE,AE_SD2=AE_SD*2.0, RRD=RRD)

    return(list(err=err,table=table))
}


## Load grids and maps
if (FALSE) {
    g1 = get.nc("maps/grid_CCSM3-T42a.nc")
    # g2 = get.nc("maps/grid_HIM-20KM.nc")
    g3 = get.nc("maps/grid_CCSM3-T42b.nc")

    # m1 = get.nc("maps/map_CCSM3-T42_HIM-20KM_20.nc")
    # m2 = get.nc("maps/map_HIM-20KM_CCSM3-T42_20.nc")

    stats = grid_stats("Hs",g1$Hs,g3$Hs,g3$mask)

}


## Check neighbors and quadrants
if (FALSE) {
    i = 20
    j = 22

    now  = data.frame(lon=g2$lon2D[i,j],lat=g2$lat2D[i,j])
    now1 = data.frame(lon=g1$lon2D[m1$i[i,j,]],lat=g1$lat2D[m1$i[i,j,]], 
                     dist=round(m1$dist[i,j,]*1e-3,1),quadrant=m1$quadrant[i,j,],weight=m1$weight[i,j,]*1e10)
    # now  = data.frame(lon=m2$lon2D[i,j],lat=m2$lat2D[i,j])
    # now1 = data.frame(lon=g1$lon2D[m2$i[i,j,]],lat=g1$lat2D[m2$i[i,j,]], 
    #                  dist=round(m2$dist[i,j,]*1e-3,1),quadrant=m2$quadrant[i,j,],weight=m2$weight[i,j,]*1e10)
    # now  = data.frame(lon=m3$x2D[i,j],lat=m3$y2D[i,j])
    # now1 = data.frame(lon=g2$x2D[m3$i[i,j,]],lat=g2$y2D[m3$i[i,j,]], 
    #                  dist=round(m3$dist[i,j,]*1e-3,1),quadrant=m3$quadrant[i,j,],weight=m3$weight[i,j,]*1e10)
    
    xlim = range(now1$lon)
    ylim = range(now1$lat)

    cols = now1$quadrant
    cols[cols==0] = 8 

    pchs = paste(now1$quadrant)

    plot(now$lon,now$lat,type="n",xlim=xlim,ylim=ylim)
    abline(v=now$lon,h=now$lat)
    points(now1$lon,now1$lat,col=cols,pch=pchs)
    legend("bottomleft",paste(c(1:4)),col=c(1:4),pch=1)
}

