library("lattice")
library("RColorBrewer")
library("colorRamps")

dat <- read.table("../../data/summary_vary_sd_vs_qjuv_max_eta_.csv", sep=";",header=T)

divide.it <- seq(-1,9,1)
subset.cor <- dat[dat$sd_hc_noise >= 0.4 & dat$sd_hc_noise <= 0.43,]

f.abs <- function(x) {
    return(x)
}

# plot correlations across the qjuv gradient
pdf("corrs_qjuv_gradient.pdf")
print(xyplot(f.abs(envt_cor_phen_prestige_horiz) + f.abs(envt_cor_phen_prestige_vert) + f.abs(envt_cor_xconformist_horiz) + f.abs(envt_cor_xconformist_vert) + f.abs(envt_cor_cue_juv_envt_high) ~ qjuv | juvenile_learns_remote_envt * (1-p) * envt_change_at_birth * mu_hp * mu_vp
                ,data=subset.cor
                ,xlab="qjuv"
                ,ylab="corr"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                ))
dev.off()
print("correlations done.")

pdf("slopes_qjuv.pdf",width=20)
print(xyplot(
                + abs(mean_hp)
                + abs(mean_vp)
                + abs(mean_hc)
                + abs(mean_vc)
                + abs(mean_ajuv)
                ~ qjuv | juvenile_learns_remote_envt * (1-p) * envt_change_at_birth * mu_hp * mu_vp
                ,data=subset.cor
                ,xlab="qjuv"
                ,ylab="Mean cue sensitivities"
                #                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )

dev.off()


# overview of all levelplots
pdf("levelplots_overview.pdf")
print(
        levelplot(eta2_max ~ qjuv * sd_vc_noise | juvenile_learns_remote_envt * p * envt_change_at_birth * mu_hp * mu_vp
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,data=dat
                ,at=divide.it
                ,col.regions=brewer.pal(n=length(divide.it),name="Set3")
                ))
dev.off()

pdf("levelplots_overview_cloud.pdf")
print(
        cloud(eta2_max ~ qjuv * sd_vc_noise | juvenile_learns_remote_envt * p * envt_change_at_birth * mu_hp * mu_vp
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,data=dat
                ,default.scales=list(arrows=F)
                ,at=divide.it
                #    ,col.regions=brewer.pal(n=length(divide.it),name="Set3")
                ))
dev.off()


pdf("levelplots_overview_delta.pdf")
print(
        levelplot(delta_eta2 ~ qjuv * sd_vc_noise | juvenile_learns_remote_envt * p * envt_change_at_birth * mu_hp * mu_vp
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,data=dat
                ,col.regions=matlab.like
                ))
dev.off()

pdf("levelplots_overview_2nd_eta_max.pdf")
print(
        levelplot(eta2_2nd ~ qjuv * sd_vc_noise | juvenile_learns_remote_envt * p * envt_change_at_birth * mu_hp * mu_vp
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,data=dat
                ,at=divide.it
                ,col.regions=brewer.pal(n=length(divide.it),name="Set3")
                ))
dev.off()
stop()
juv_learn_u <- sort(unique(dat$juvenile_learns_remote_envt))
p_u <- sort(unique(dat$p))
envt_change_u <- sort(unique(dat$envt_change_at_birth))

for (juv_learn_i in juv_learn_u)
{
    for (p_i in p_u)
    {
        for (envt_change_i in envt_change_u)
        {
            s.dat <- dat[
                            dat$juvenile_learns_remote_envt == juv_learn_i
                            &
                            dat$p == p_i
                            &
                            dat$envt_change_at_birth == envt_change_i,]

            print(nrow(s.dat))


            svg(file=paste("single_strip_juv_learn_",juv_learn_i,"_p_",p_i,"_envt_change_birth_",envt_change_i,".svg",sep="")
                   ,family="Times" ,width=4,height=4
                    )
            print(
                    levelplot(eta2_max ~ qjuv * sd_vc_noise
                            ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                            ,data=s.dat
                            ,at=divide.it
                            ,xlim=c(0.5,1)
                            ,ylim=c(0,1)
                            ,xlab=""
                            ,ylab=""
                            ,colorkey=F
                            ,fontsize=16
                            ,col.regions=brewer.pal(n=length(divide.it),name="Set3")
                            ))
            dev.off()
            
        }
    }
}

