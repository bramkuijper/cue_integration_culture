library("lattice")
library("latex2exp")
library("ggplot2")
library("reshape2")

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

dat <- read.table("../../data/summary_social_learn_qmat.csv"
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

dat <- read.table("../../data/summary_vary_p_social.csv"
                  ,sep=";"
                  ,header=T)

pdf("overview_a_vary_p.pdf")
print(
      xyplot(mean_agen + mean_asoc + mean_ajuv + mean_amat ~ p| qmat * qjuv * m * nc * np
             ,data=dat
             ,pch=21
             ,cex=0.5
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,auto.key=T)
      )
dev.off()
pdf("overview_dp_dc_vary_p.pdf")
print(
      xyplot(mean_dc + mean_dp ~ p| qmat * qjuv * m * nc * np
             ,data=dat
             ,pch=21
             ,cex=0.5
             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
             ,auto.key=T)
      )
dev.off()

# now vary nc
dat <- read.table("../../data/summary_social_vary_nc.csv"
                  ,sep=";"
                  ,header=T)
print("ggplotting")

# try ggplot as we need to forwart it 2 Sasha
library("ggplot2")

dat.long <- melt(dat
        ,id.vars=c("p","qmat","qjuv","nc")
        ,measure=c("mean_agen", "mean_asoc", "mean_ajuv", "mean_amat"))

# some semi-fancy labels not really
xlab <- TeX("Rate of environmental change, $1-p$")
ylab <- TeX("Cue weight, $a$")

str(dat.long)

the.plot <- ggplot(dat.long
        ,aes(x=1 - p
                ,y=value
                ,colour=variable
        )
) + geom_smooth() + labs(x = "Rate of environmental change"
                            ,y="Cue weighting")

the.plot + facet_wrap(vars(qmat, qjuv, nc), labeller=label_both)

ggsave(filename="overview_a_vary_nc.pdf", width=40, height=40, units="cm")

#pdf("overview_a_vary_nc.pdf",width=13,height=13)
#print(
#      xyplot(mean_agen + mean_asoc + mean_ajuv + mean_amat ~ (1-p) | qmat * qjuv * m * np * nc
#             ,data=dat
#             ,pch=21
#             ,cex=0.8
#             ,xlab=xlab
#             ,ylab=ylab
#             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
#             ,auto.key=T)
#      )
#dev.off()
#pdf("overview_dp_dc_vary_nc.pdf",width=13,height=13)
#print(
#      xyplot(mean_dc + mean_dp ~ p | qmat * qjuv * m * np * nc
#             ,data=dat
#             ,pch=21
#             ,xlab=xlab
#             ,cex=0.8
#             ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
#             ,auto.key=T)
#      )
#dev.off()
