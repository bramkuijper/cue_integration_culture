library("lattice")

dat <- read.table("../../data/summary_social_learn.csv"
                  ,sep=";"
                  ,header=T)

pdf("overview_b_and_d.pdf")
print(
      xyplot(mean_dp + mean_dc ~ p | qmat * qjuv * m * nc * np
             ,data=dat
             ,pch=21
             ,cex=0.5
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,auto.key=T)
      )
dev.off()

pdf("overview_a.pdf")
print(
      xyplot(mean_agen + mean_asoc + mean_ajuv + mean_amat ~ p | qmat * qjuv * m * nc * np
             ,data=dat
             ,pch=21
             ,cex=0.5
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,auto.key=T)
      )
dev.off()

dat <- read.table("../../data/summary_vary_qmat.csv"
                  ,sep=";"
                  ,header=T)

pdf("overview_a_vary_qmat.pdf")
print(
      xyplot(mean_agen + mean_asoc + mean_ajuv + mean_amat ~ qmat | p * qjuv * m * nc * np
             ,data=dat[dat$qmat >= 0.5,]
             ,pch=21
             ,cex=0.5
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,auto.key=T)
      )
dev.off()

str(dat)
