# simple plot to give an overview of varying nc versus np


the.data <- read.table("../../data/summary_vary_nc_np.csv",sep=";",header=T)

print(table(the.data[,c("npv","ncv","nph","nch")]))

the.data.s <- subset(the.data, ncv == 6 - npv & p >= 0.5)

stopifnot(nrow(the.data.s) > 0)


pdf("mean_a_for_nc_vs_np.pdf")
print(
        xyplot(mean_ajuv + mean_amat + mean_asoc_vert + mean_asoc_horiz ~ npv | qjuv * qmat * p
        ,data=the.data.s
        ,auto.key=T
        ,strip=function(strip.levels,...) {strip.default(strip.levels=T,...) }
        )
        )
dev.off()

pdf("mean_a_for_nc_vs_np.pdf")
print(
        xyplot(mean_hp + mean_hc + mean_vp + mean_vc ~ npv | qjuv * qmat * p
        ,data=the.data.s
        ,auto.key=T
        ,strip=function(strip.levels,...) {strip.default(strip.levels=T,...) }
        )
        )
dev.off()
