
library(myr)




# Load data
if (TRUE) {

    dat = read.table("fort.13",header=FALSE,skip=2,
                col.names=c("xc","yc","F","Z","ZX","ZY"))
    
    #new = read.table("fort.13_natnew",header=FALSE,skip=0,
    #                 col.names=c("xc","yc","Z"))

    xlim = range(dat$xc)
    ylim = range(dat$yc)
    zlim = range(dat$F,dat$Z)
    zlim2 = c(-0.15,0.15)

    breaks = pretty(zlim,30)
    col    = colorRampPalette(jet.colors)(length(breaks)-1)

    nx = 25 
    ny = 25 

    par(mfrow=c(2,3))
    par(plt=c(0.1,0.8,0.05,0.95))
    quilt.plot(dat$xc,dat$yc,dat$F,xlim=xlim,ylim=ylim,zlim=zlim,nx=nx,ny=ny,breaks=breaks,col=col)
    quilt.plot(dat$xc,dat$yc,dat$Z,xlim=xlim,ylim=ylim,zlim=zlim,nx=nx,ny=ny,breaks=breaks,col=col)
    quilt.plot(dat$xc,dat$yc,dat$Z-dat$F,xlim=xlim,ylim=ylim,zlim=zlim2,nx=nx,ny=ny)
    
    #quilt.plot(new$xc,new$yc,new$Z,xlim=xlim,ylim=ylim,nx=nx,ny=ny,breaks=breaks,col=col)

}
