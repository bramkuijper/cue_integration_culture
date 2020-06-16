# simple plot to give an overview of varying nc versus np
library("RColorBrewer")


# get the directory name of this script
# we use this to get the directory name of the
# data folder
script.dir <- dirname(sys.frame(1)$ofile)

# the filename of the data file 
filename <- "summary_vary_n_max_eta_.csv"

full_filename = file.path(script.dir,"../../data",filename)

dat <- read.table(full_filename, sep=";",header=T)

cue.intervals <- seq(-1,9,1)
colors <- brewer.pal(n=length(cue.intervals),name="Set3")

pdf("contour_vary_n.pdf")
print(
        levelplot(eta2_max ~ qjuv * sd_vc_noise | nch
                #                ,strip=myStrip
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,data=dat
                ,at=cue.intervals
                ,xlim=c(0.5,1)
                ,ylim=c(0,1)
                #                ,xlab=list(label=expression(paste("Probability environment changes, ",1-{},italic("p"))),cex=1.5)
                #                            ,ylab=as.list("Noise in vertical versus horizontal",expression(paste("social learning, ",sigma["v"],"=",1-sigma["h"],sep="")))


                ,ylab=list(label=expression(atop("Noise in vertical versus horizontal",
               paste("social learning, ",sigma["v"],"=",1-sigma["h"],sep=""))),cex=1.5)
                ,colorkey=F
                ,fontsize=16
                ,col.regions=colors
                ,scales=list(
                    y=list(
                      at=c(0,0.2,0.4,0.6,0.8,1.0)
                      ),
                    x=list(
                      at=c(0,0.2,0.4,0.6,0.8,1.0)
                      )
                    )
                ))
dev.off()
