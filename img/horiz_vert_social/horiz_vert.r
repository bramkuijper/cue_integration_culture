# script to generate a simple overview
# plot of the evolved cues when individuals rely on
# both horizontal and social learning
library("ggplot2")
library("reshape")


# read in the data
the.data <- read.table("../../data/summary_cue_int_big.csv",sep=";",header=T)

ylim <- c(-0.5,7)

pdf("overview_horiz_vert_a_weightings.pdf")
print(xyplot(mean_ajuv + mean_asoc_horiz + mean_asoc_vert ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()

dat.long <- melt(the.data
        ,id.vars=c("p","qjuv","sdsoc_horiz","sdsoc_vert")
        ,measure=c("mean_ajuv","mean_asoc_horiz","mean_asoc_vert","mean_agen","mean_amat"))

str(dat.long)

the.plot <- ggplot(dat.long
        ,aes(x=1 - p
                ,y=value
                ,colour=variable
        )
) + geom_line() + geom_point() + labs(x = "Rate of environmental change"
                            ,y="Cue weighting")

the.plot + facet_wrap(vars(qjuv,sdsoc_horiz,sdsoc_vert), labeller=label_both)

ggsave(filename="overview_horiz_vert_a_weightings_ggplot.pdf", width=40, height=40, units="cm")


pdf("overview_horiz_vert_a_weightings2.pdf")
print(xyplot(mean_agen + mean_amat ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()

pdf("overview_horiz_vert_v.pdf")
print(xyplot(mean_vp + mean_vc ~  (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()

pdf("overview_horiz_vert_h.pdf")
print(xyplot(mean_hp + mean_hc ~  (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()


dat.long <- melt(the.data
        ,id.vars=c("p","qjuv","sdsoc_horiz","sdsoc_vert")
        ,measure=c("mean_hp","mean_hc"))

str(dat.long)

the.plot <- ggplot(dat.long
        ,aes(x=1 - p
                ,y=value
                ,colour=variable
        )
) + geom_line() + geom_point() + labs(x = "Rate of environmental change"
                            ,y="Horizontal cue weighting")

the.plot + facet_wrap(vars(qjuv,sdsoc_horiz,sdsoc_vert), labeller=label_both)

ggsave(filename="overview_horiz_vert_h_weightings_ggplot.pdf", width=40, height=40, units="cm")

pdf("overview_horiz_vert_phenotypic_variance.pdf")
print(xyplot(var_phen_ad ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()
