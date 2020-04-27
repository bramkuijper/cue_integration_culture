#!/usr/bin/env Rscript

library(heplots)
library(data.table)

args <- commandArgs(trailingOnly=T)

dir <- args[[1]]

flist <- list.files(path=dir
        ,pattern="sim_.*\\d_dist.csv$")

prepend.eta <- function(x) {
    return(paste("eta2_",x,sep=""))
}

for (file_i in flist)
{
    df = read.table(file_i,sep=";",header=T)

    df[,"bmat_phen_X_phen_mat_error"] = df[,"bmat_phen"] * (df[,"phen_mat_error"] - 0.5)

    df[,"bmat_envt_X_maternal_envt_cue_eror"] = df[,"bmat_envt"] * (df[,"maternal_envt_cue_error"] - 0.5)

    df[,"agen_X_g"] = df[,"agen"] * df[,"g"]
    df[,"ajuv_X_cue_juv_envt_high"] = df[,"ajuv"] * (df[,"cue_juv_envt_high"] - 0.5)

    df[,"vc_X_xconformist_vert_error"] = df[,"vc"] * (df[,"xconformist_vert_error"] - 0.5)
    df[,"vp_X_phen_prestige_vert_error"] = df[,"vp"] * (df[,"phen_prestige_vert_error"] - 0.5)

    df[,"hp_X_phen_prestige_horiz_error"] = df[,"hp"] * (df[,"phen_prestige_horiz_error"] - 0.5)
    df[,"hc_X_xconformist_horiz_error"] = df[,"hc"] * (df[,"xconformist_horiz_error"] - 0.5)

    columns <- c(
        "aintercept",
        "bmat_phen_X_phen_mat_error",
        "bmat_envt_X_maternal_envt_cue_eror",
        "agen_X_g",
        "ajuv_X_cue_juv_envt_high",
        "vc_X_xconformist_vert_error",
        "vp_X_phen_prestige_vert_error",
        "hp_X_phen_prestige_horiz_error",
        "hc_X_xconformist_horiz_error")

    ind.vars <- c()
    n.ind.vars <- c()
    for (col_i in columns)
    {
        if (var(df[,col_i]) > 0.0)
        {
            ind.vars <- c(ind.vars,col_i)
        } else
        {
            n.ind.vars <- c(n.ind.vars,col_i)
        }
    }

    formula <- paste("phen_ad_logistic ~ ",paste(ind.vars,collapse="+"))

    mod1 <- lm(as.formula(formula)
        ,data=df)

    etas <- etasq(mod1, partial=FALSE)

    eta.df <- as.data.frame(t(etas))

    # add the omitted columns
    eta.df[,n.ind.vars] <- 0

    eta.df.n <- unlist(lapply(X=names(eta.df), FUN=prepend.eta))

    names(eta.df) <- eta.df.n
    
    eta.df <- setcolorder(eta.df,sort(eta.df.n))

    write.table(x=eta.df
            ,sep=";"
            ,row.names=F
            ,file=paste(file_i,"etasq",sep=""))
}
