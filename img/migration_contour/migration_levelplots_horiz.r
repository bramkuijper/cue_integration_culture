library("lattice")
library("grid")
library("RColorBrewer")
library("colorRamps")
library("gtools")
source("/Users/bram/R/src/bramlib.r")

type <- "pdf"

tickcex <- 0.75
labelcex <- 1
tick.lwd <- 0.5
lines.lwd <- 0.5

level.dividers <- seq(-1,9,1)

single.level <- function(
        row
        ,col
        ,dataset
        ,y.ticks=T
        ,x.ticks=T
        ,y.label=F
        ,x.label=F
        ,ind.label=""
        ,title="")
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

        if (ind.label != "")
        {
            grid.text(x=unit(units="native",x=0.55)
                    ,y=unit(units="native",x=0.95)
                    ,just="center"
                    ,label=ind.label)
        }

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
    
    # plot title

    if (typeof(title) == typeof(expression("")))
    {
        pushViewport(viewport(layout.pos.row=row-1,
                                layout.pos.col=col,
                                xscale=lplot$x.limits,
                                yscale=lplot$y.limits
                                ))
            grid.text(x=0.5
                    ,y=0.5
                    ,just="center"
                    ,label=title)
        upViewport()
    }

    
    pushViewport(viewport(layout.pos.row=row+1,
                            layout.pos.col=col,
                            xscale=lplot$x.limits,
                            yscale=lplot$y.limits
                            ))

    text <- ""

    if (x.label)
    {
        text <- expression(paste("Noise in maternal cue versus individual learning, ",italic("q")["mat"]," = ",1-{},italic("q")["ind"]))
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


legend <- function(row,col)
{
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col,
                            ))

    labels <- c("Horizontal, prestige"
            ,"Horizontal, conformity"
            ,"Vertical, prestige"
            ,"Vertical, conformity"
            ,"Individual learning"
            ,"Maternal phenotype"
            ,"Maternal environment")

    colors <- c("#c865a1"
            ,"#d9d9d9"
            ,"#fcbbcf"
            ,"#b3de68"
            ,"#bdb8d9"
            ,"#7fb0d1"
            ,"#f97f70")

    y.i = 0.95
    line.height <- 0.1
    box.hw <- 0.07

    for (i in 1:length(labels))
    {
        grid.rect(x=unit(units="native",x=0.1)
                ,y=unit(units="native",x=y.i)
                ,width=unit(units="native",x=box.hw)
                ,height=unit(units="native",x=box.hw)
                ,gp=gpar(lwd=lines.lwd,fill=colors[[i]])
                )

        grid.text(
                ,x=unit(units="native",x=0.18)
                ,y=unit(units="native",x=y.i)
                ,just="left"
                ,label=labels[[i]])

        y.i <- y.i - line.height
    }

    upViewport()
}


# get the directory name of this script
# we use this to get the directory name of the
# data folder
script.dir <- dirname(sys.frame(1)$ofile)


filename.horiz_pre <- "summary_horiz_pre_disperse_max_eta_.csv"
full_filename.horiz_pre = file.path(script.dir,"../../data",filename.horiz_pre)
the.data.horiz_pre <- read.table(full_filename.horiz_pre, sep=";",header=T)

filename.ind_pre_m <- "summary_ind_learn_post_m_max_eta_.csv"
full_filename.ind_pre_m = file.path(script.dir,"../../data",filename.ind_pre_m)
the.data.ind_pre_m <- read.table(full_filename.ind_pre_m, sep=";",header=T)

the.data.horiz_pre <- the.data.horiz_pre[the.data.horiz_pre$envt_change_at_birth == 0,]
the.data.ind_pre_m <- the.data.ind_pre_m[the.data.ind_pre_m$envt_change_at_birth == 0,]

init.plot("levelplot_vary_migration_horiz", 
                type=type,
                width=820,
                height=800,
                font="times")

plot.dim <- 0.75

widths <- c(0.5,plot.dim,0.1,plot.dim,0.1,plot.dim,1)

heights <- c(0.1,plot.dim,0.1,plot.dim,0.1,plot.dim,0.1,plot.dim,0.2)

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

############################################

# HORIZONTAL LEARNING PRIOR TO DISPERSAL

############################################


    subs.1 <- subset(the.data.horiz_pre, m == 0.3 & autocorr == 0.6)

    plot_x <- single.level(
        row=2
        ,col=2
        ,dataset=subs.1
        ,y.label=F
        ,y.ticks=T
        ,x.label=F
        ,x.ticks=F
        ,ind.label="A"
        ,title=expression(paste("Migration probability, ",italic("d")," = ",0.3))
        )

    subs.2 <- subset(the.data.horiz_pre, m == 0.5 & autocorr == 0.6)
    
    plot_x <- single.level(
        row=2
        ,col=4
        ,dataset=subs.2
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F
        ,ind.label="B"
        ,title=expression(paste("Migration probability, ",italic("d")," = ",0.5))
        )

    subs.3 <- subset(the.data.horiz_pre, m == 0.8 & autocorr == 0.6)
    
    plot_x <- single.level(
        row=2
        ,col=6
        ,dataset=subs.3
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F
        ,ind.label="C"
        ,title=expression(paste("Migration probability, ",italic("d")," = ",0.8))
        )

    subs.4 <- subset(the.data.horiz_pre, m == 0.3 & autocorr == -0.6)
    
    plot_x <- single.level(
        row=4
        ,col=2
        ,dataset=subs.4
        ,y.label=F
        ,y.ticks=T
        ,x.label=F
        ,x.ticks=F
        ,ind.label="D"
        ,title=""
        )

    subs.5 <- subset(the.data.horiz_pre, m == 0.5 & autocorr == -0.6)
    
    plot_x <- single.level(
        row=4
        ,col=4
        ,dataset=subs.5
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F
        ,ind.label="E"
        ,title=""
        )

    subs.6 <- subset(the.data.horiz_pre, m == 0.8 & autocorr == -0.6)
    
    plot_x <- single.level(
        row=4
        ,col=6
        ,dataset=subs.6
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F
        ,ind.label="F"
        ,title=""
        )



##### 
#
# VERTICAL POST MIGRATION
#
#####

row <- 6

ind.label.i <- 71

    subs.1 <- subset(the.data.ind_pre_m, m == 0.3 & autocorr == 0.6)

    plot_x <- single.level(
        row=row
        ,col=2
        ,dataset=subs.1
        ,y.label=F
        ,y.ticks=T
        ,x.label=F
        ,x.ticks=F
        ,ind.label=chr(ind.label.i)
        ,title=""
        )

    subs.2 <- subset(the.data.ind_pre_m, m == 0.5 & autocorr == 0.6)
    
ind.label.i <-  ind.label.i + 1
    plot_x <- single.level(
        row=row
        ,col=4
        ,dataset=subs.2
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F
        ,ind.label=chr(ind.label.i)
        ,title=""
        )

    subs.3 <- subset(the.data.ind_pre_m, m == 0.8 & autocorr == 0.6)
    
ind.label.i <-  ind.label.i + 1
    plot_x <- single.level(
        row=row
        ,col=6
        ,dataset=subs.3
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=F
        ,ind.label=chr(ind.label.i)
        ,title=""
        )

row <- 8

    subs.4 <- subset(the.data.ind_pre_m, m == 0.3 & autocorr == -0.6)
    
ind.label.i <-  ind.label.i + 1
    plot_x <- single.level(
        row=row
        ,col=2
        ,dataset=subs.4
        ,y.label=F
        ,y.ticks=T
        ,x.label=F
        ,x.ticks=T
        ,ind.label=chr(ind.label.i)
        ,title=""
        )

    subs.5 <- subset(the.data.ind_pre_m, m == 0.5 & autocorr == -0.6)
    
ind.label.i <-  ind.label.i + 1
    plot_x <- single.level(
        row=row
        ,col=4
        ,dataset=subs.5
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=T
        ,ind.label=chr(ind.label.i)
        ,title=""
        )

    subs.6 <- subset(the.data.ind_pre_m, m == 0.8 & autocorr == -0.6)
   
ind.label.i <-  ind.label.i + 1
    plot_x <- single.level(
        row=row
        ,col=6
        ,dataset=subs.6
        ,y.label=F
        ,y.ticks=F
        ,x.label=F
        ,x.ticks=T
        ,ind.label=chr(ind.label.i)
        ,title=""
        )

# indication that it is vertical pre, horizontal post
    grid.text(
            x=unit(units="native"
                    ,x=0.01)
            ,y=unit(units="native"
                    ,x=0.8)
            ,rot=90
            ,just="centre"
            ,hjust="centre"
            ,label=expression(paste({}%down%{}," Vertical + horizontal pre migration ",{}%down%{}))
           )

# indication that it is vertical post, horizontal post
    grid.text(
            x=unit(units="native"
                    ,x=0.01)
            ,y=unit(units="native"
                    ,x=0.23)
            ,rot=90
            ,just="centre"
            ,hjust="centre"
            ,label=expression(paste({}%down%{}," Individual learning post migration ",{}%down%{}))
           )

# y label

    grid.text(
            x=unit(units="native"
                    ,x=0.06)
            ,y=unit(units="native"
                    ,x=0.5)
            ,rot=90
            ,just="centre"
            ,hjust="centre"
            ,label=expression(atop("Noise in vertical versus horizontal",
                                   paste("social learning, ",sigma["v"]," = ",1-sigma["h"],sep="")))
           )

    grid.text(
            x=unit(units="native"
                    ,x=0.43)
            ,y=unit(units="native"
                    ,x=0.0)
            ,rot=0
            ,just="centre"
            ,hjust="centre"
            ,label=expression(paste("Fidelity of individually learned cues versus maternal cues, ",italic("q")["ind"]," = ",1.5-{},italic("q")["mat"]))
           )

   legend(row=2,col=7)
    
   
   x.pos.arrow.text <- 0.75
   grid.text(
            x=unit(units="native"
                    ,x=x.pos.arrow.text)
            ,y=unit(units="native"
                    ,x=0.8)
            ,just="left"
            ,label=expression(paste({}%<-%{},"Rate of change ",1 - {},italic("p")," = 0.2",sep=""))
           )
   
   grid.text(
            x=unit(units="native"
                    ,x=x.pos.arrow.text)
            ,y=unit(units="native"
                    ,x=0.55)
            ,just="left"
            ,label=expression(paste({}%<-%{},"Rate of change ",1 - {},italic("p")," = 0.8",sep=""))
           )
   grid.text(
            x=unit(units="native"
                    ,x=x.pos.arrow.text)
            ,y=unit(units="native"
                    ,x=0.32)
            ,just="left"
            ,label=expression(paste({}%<-%{},"Rate of change ",1 - {},italic("p")," = 0.2",sep=""))
           )
   
   grid.text(
            x=unit(units="native"
                    ,x=x.pos.arrow.text)
            ,y=unit(units="native"
                    ,x=0.08)
            ,just="left"
            ,label=expression(paste({}%<-%{},"Rate of change ",1 - {},italic("p")," = 0.8",sep=""))
           )
   
   upViewport()


exit.plot()

