# script to generate a simple overview
# plot of the evolved cues when individuals rely on
# both horizontal and social learning
library("ggplot2")
library("reshape")

# I am currently plotting a scenario
# where we removed all the 'sub'-sigmoidals as
# in Leimar et al 


# read in the data
# I
if (!exists("the.data"))
{
    the.data <- read.table("../../data/summary_single_logistic.csv",sep=";",header=T)
}

# little function to find column names 
# (as there is a massive amount of columns)
findcol <- function(pattern)
{
    return(names(the.data)[grep(pattern, names(the.data))])
}

names_prop <- findcol("var_prop.*")

stopifnot(length(names_prop) > 0)

the.formula = paste(names_prop, collapse=" + ")

the.formula = paste(the.formula, " ~ (1.0 - p) | mu_bmat_envt * mu_hp * qjuv * qmat")

pdf("slopes.pdf")
print(xyplot(mean_bmat_phen 
                + mean_bmat_envt
                + mean_agen
                + mean_ajuv
                + mean_aintercept
                + mean_vc
                + mean_vp
                + mean_hp
                + mean_hc ~ (1-p) | mu_bmat_envt * mu_hp * qjuv * qmat
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,ylab="Variance components"
                #                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()

# print a sample file for further inspection
sample.series <- subset(the.data,qmat == 1.0 & qjuv == 0.5 & mu_hp == 0.01 & mu_bmat_envt == 0.01 & p > 0.05 & p < 0.07)

stopifnot(nrow(sample.series) > 0)

print('pass')

pdf("var_components.pdf")
print(xyplot(as.formula(the.formula)
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                ,ylab="Variance components"
                #                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()
stop()


# agen var comp
matching_columns.gen <- get_columns_pattern(
        pattern="^(var|cov)_agen_X_g.*"
        ,data=the.data)
sum_var_expression.gen <- paste(matching_columns.gen,collapse="+")
the.data[,"var_comp_gen"] <- with(the.data,eval(parse(text=sum_var_expression.gen)))



# juvenile cue
matching_columns.juvcue <- get_columns_pattern(
        pattern="^(var|cov)_ajuv_X_cue_juv_envt_high.*"
        ,data=the.data)

sum_var_expression.juvcue <- paste(matching_columns.juvcue,collapse="+")
the.data[,"var_comp_juvcue"] <- with(the.data,eval(parse(text=sum_var_expression.juvcue)))


# maternal phenotypic cue
matching_columns.matphen <- get_columns_pattern(
        pattern="^(var|cov)_bmat_phen_X_phen_mat_error.*"
        ,data=the.data)

print(the.data[the.data$p == 1.0,matching_columns.matphen])

sum_var_expression.matphen <- paste(matching_columns.matphen,collapse="+")
the.data[,"var_comp_matphen"] <- with(the.data,eval(parse(text=sum_var_expression.matphen)))

# maternal environmental cue
matching_columns.matenvt <- get_columns_pattern(
        pattern="^(var|cov)_bmat_envt_X_maternal_envt_cue_eror.*"
        ,data=the.data)

sum_var_expression.matenvt <- paste(matching_columns.matenvt,collapse="+")
the.data[,"var_comp_matenvt"] <- with(the.data,eval(parse(text=sum_var_expression.matenvt)))



pdf("nonsocial_variance_components.pdf")
print(xyplot(
                var_comp_gen
                + var_comp_juvcue
                + var_aintercept
                + var_phen_ad_logistic
                + var_comp_matphen 
                + var_comp_matenvt 
                ~ (1.0 - p) 
                ,data=the.data
                ,xlab="Probability environment changes, 1 - p"
                #                ,ylim=ylim
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,auto.key=T
                )
        )
dev.off()


stop()

dat.long <- melt(the.data
        ,id.vars=c("p","qjuv","sdsoc_horiz","sdsoc_vert", "sdmat","juvenile_survival")
        ,measure=c("mean_ajuv","mean_asoc_horiz","mean_asoc_vert","mean_agen","mean_amat"))

the.plot <- ggplot(dat.long
        ,aes(x=1 - p
                ,y=value
                ,colour=variable
        )
) + geom_point() + labs(x = "Rate of environmental change"
                            ,y="Cue weighting")

the.plot + facet_wrap(vars(qjuv,sdmat,sdsoc_horiz,sdsoc_vert,juvenile_survival), labeller=label_both)

ggsave(filename="overview_horiz_vert_a_weightings_ggplot.pdf", width=60, height=60, units="cm")


pdf("overview_horiz_vert_a_weightings2.pdf")
print(xyplot(mean_agen + mean_amat ~ (1.0 - p) | qjuv * sdsoc_horiz * sdsoc_vert *sdmat
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
