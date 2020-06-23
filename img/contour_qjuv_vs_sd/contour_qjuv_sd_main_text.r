library("lattice")
library("grid")
library("RColorBrewer")
library("colorRamps")
source("/Users/bram/R/src/bramlib.r")

type <- "svg"

tickcex <- 0.75
labelcex <- 1
tick.lwd <- 0.5
lines.lwd <- 0.5

level.dividers <- seq(-1,9,1)

single.level <- function(
        row,
        col,
        dataset,
        y.ticks=T,
        x.ticks=T,
        y.label=F,
        x.label=F)
{
    lplot <- levelplot(
            eta2_max ~ qjuv * sd_vc_noise
            ,data=dataset
            ,at=level.dividers
            ,col.regions=brewer.pal(
                    n=length(level.dividers)
                    ,name="Set3"))
    
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col,
                            xscale=lplot$x.limits,
                            yscale=lplot$y.limits
                            ))

        do.call("panel.levelplot",trellis.panelArgs(lplot,1))
        grid.rect(gp=gpar(lwd=lines.lwd,fill="transparent"))

    upViewport()
    
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col-1,
                            xscale=lplot$x.limits,
                            yscale=lplot$y.limits
                            ))

    text <- ""

    if (y.label)
    {
        text  <- list(label=expression(atop("Noise in vertical versus horizontal",
       paste("social learning, ",sigma["v"],"=",1-sigma["h"],sep=""))),cex=labelcex)
    }

    single.axis(
            range=lplot$y.limits
            ,labels=y.ticks
            ,side="right"
            ,nmain=5
            ,nsub=6
            ,text=text
            ,x.text.off=0.1
            ,cex=tickcex
            ,labelcex=labelcex
            ,tck=0.5
            ,lwd=tick.lwd)

    upViewport()
    
    pushViewport(viewport(layout.pos.row=row+1,
                            layout.pos.col=col,
                            xscale=lplot$x.limits,
                            yscale=lplot$y.limits
                            ))

    text <- ""

    if (x.label)
    {
        text <- list(label=expression(paste("Noise in maternal cue versus individual learning, ",italic("q")["mat"]," = ",1-{},italic("q")["ind"])),cex=labelcex)
    }

    single.axis(
            range=lplot$x.limits
            ,side="top"
            ,labels=x.ticks
            ,nmain=5
            ,nsub=6
            ,text=text
            ,x.text.off=0.1
            ,cex=tickcex
            ,labelcex=labelcex
            ,tck=0.5
            ,lwd=tick.lwd)

    upViewport(1)
}


# get the directory name of this script
# we use this to get the directory name of the
# data folder
script.dir <- dirname(sys.frame(1)$ofile)

filename <- "summary_vary_sd_vs_qjuv_max_eta_.csv"

full_filename = file.path(script.dir,"../../data",filename)

the.data <- read.table(full_filename, sep=";",header=T)


init.plot("levelplot_qjuv_vs_sdsoc", 
                type=type,
                width=330,
                height=300,
                font="times")

widths <- c(0.3,1,0.1,1,0.3)
heights <- c(0.1,1,0.1,1,0.3)

# initial viewport
pushViewport(
            viewport(name="vp_head",
                just="center",
            height=.95,
            width=.95,
            layout=grid.layout(nrow=length(heights),
                                ncol=length(widths),
                                widths=widths,
                                heights=heights)))

    subs.1 <- subset(the.data, p == 0.8 & envt_change_at_birth == 0)

    plot_x <- single.level(
        row=2
        ,col=2
        ,dataset=subs.1
        ,y.label=F
        ,y.ticks=T
        ,x.label=F
        ,x.ticks=F)

    subs.2 <- subset(the.data, p == 0.8 & envt_change_at_birth == 1)
    
    plot_x <- single.level(
        row=4
        ,col=2
        ,dataset=subs.2
        ,y.label=F
        ,y.ticks=T
        ,x.label=F
        ,x.ticks=T)

    subs.3 <- subset(the.data, p == 0.2 & envt_change_at_birth == 0)
    
    plot_x <- single.level(
        row=2
        ,col=4
        ,dataset=subs.3
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F)

    subs.4 <- subset(the.data, p == 0.2 & envt_change_at_birth == 1)
    
    plot_x <- single.level(
        row=4
        ,col=4
        ,dataset=subs.4
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=T)

upViewport()
exit.plot()

