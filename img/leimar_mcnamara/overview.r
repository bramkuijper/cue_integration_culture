library("lattice")

dat <- read.table("../../data/summary_initial_runs.csv"
        ,sep=";"
        ,header=T)

pdf("overview_initial_runs_a.pdf")
print(
        xyplot(mean_ajuv + mean_amat + mean_agen ~ p | qmat * qjuv
                ,auto.key=T
                ,data=dat
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,xlab="Frequency high envt"
                ,ylab="Sensitivity"
                )
        )
dev.off()

pdf("overview_initial_runs_b.pdf")
print(
        xyplot(mean_bmat_phen + mean_bmat_envt ~ p | qmat * qjuv
                ,auto.key=T
                ,data=dat
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,xlab="Frequency high envt"
                ,ylab="Sensitivity"
                )
        )
dev.off()
