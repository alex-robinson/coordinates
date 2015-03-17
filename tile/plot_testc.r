
library(myr)




# Load data
if (TRUE) {

    dat = read.table("../fort.13",header=FALSE,skip=2,
                col.names=c("xc","yc","F","Z","ZX","ZY"))
    
    par(mfrow=c(1,3))
    quilt.plot(dat$xc,dat$yc,dat$F,nx=20,ny=20)
    quilt.plot(dat$xc,dat$yc,dat$Z,nx=20,ny=20)
    quilt.plot(dat$xc,dat$yc,dat$Z-dat$F,nx=20,ny=20)
       
}