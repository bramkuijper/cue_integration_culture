# script to generate a simple overview
# plot of the evolved cues when individuals rely on
# both horizontal and social learning

# read in the 
the.data <- read.table("../../data/summary_horiz_vert.csv",sep=";",header=T)


pdf("overview_horiz_vert_a_weightings.pdf")
print(xyplot(mean_ajuv + mean_asoc_horiz + mean_asoc_vert ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()

pdf("overview_horiz_vert_a_weightings2.pdf")
print(xyplot(mean_agen + mean_amat ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()

pdf("overview_horiz_vert_phenotypic_variance.pdf")
print(xyplot(var_phen_ad ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()
