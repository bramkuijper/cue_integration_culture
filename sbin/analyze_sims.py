#!/usr/bin/env python3

# summarize the individual-based model
import summarizesims
import argparse
import pandas as pd
import os.path
import re
import statsmodels.api as sm
from statsmodels.formula.api import ols



parser = argparse.ArgumentParser(description="Process simulation files")

parser.add_argument("--pattern"
                    ,help="Regular expression pattern to match files"
                    ,default="(sim.*\d+$|iter.*\d+)"
                    )

parser.add_argument("--path"
                    ,help="The path over which to iterate"
                    ,default=".")

parser.add_argument("--ncores"
                    ,help="Number of cores to run this"
                    ,default=1)

parser.add_argument("--max_files"
                    ,help="Whether we should process a maximum number of files"
                    ,default=None)

parser.add_argument("--sep"
                    ,help="Column separator"
                    ,default=";")

args = vars(parser.parse_args())


# auxiliary function to generate a var-cov matrix of a 
# trait distribution
def process_dist(filename):
   
    # make a dictionary to contain all the 
    # variances we are going to calculate
    var_dict = {}

    splitname = os.path.splitext(filename)
    dist_filename = splitname[0] + "_dist" + splitname[1]

    dist_df = pd.read_csv(dist_filename, sep=";")

    # correlations between the selective 
    # environment and 
    # other variables
    envt_cor_vars = [
            "phen_ad" 
            ,"phen_ad_logistic" 
            ,"phen_juv" 
            ,"phen_juv_logistic" 
            ,"phen_mat" 
            ,"phen_mat_error" 
            ,"maternal_envt_cue_error" 
            ,"phen_prestige_vert" 
            ,"phen_prestige_vert_error" 
            ,"phen_prestige_horiz" 
            ,"phen_prestige_horiz_error" 
            ,"xconformist_vert" 
            ,"xconformist_vert_error" 
            ,"xconformist_horiz" 
            ,"xconformist_horiz_error" 
            ,"g" 
            ,"envt_sel" 
            ,"envt_prev" 
            ,"cue_ad_envt_high" 
            ,"cue_juv_envt_high"]

    correlation_matrix = dist_df[envt_cor_vars].corr()

    # get the column for the selective environment
    # out of the correlation matrix
    sel_envt_col = correlation_matrix["envt_sel"]

    row_names = sel_envt_col.index.values

    for row_name_i in row_names:
        var_dict["envt_cor_" + row_name_i] = sel_envt_col[row_name_i]

    # make columns of the relevant variables
    dist_df["bmat_phen_X_phen_mat_error"] = dist_df["bmat_phen"] * (dist_df["phen_mat_error"] - 0.5)
    dist_df["bmat_envt_X_maternal_envt_cue_eror"] = \
            dist_df["bmat_envt"] * (dist_df["maternal_envt_cue_error"] - 0.5)

    dist_df["agen_X_g"] = dist_df["agen"] * dist_df["g"]
    dist_df["ajuv_X_cue_juv_envt_high"] = dist_df["ajuv"] * (dist_df["cue_juv_envt_high"] - 0.5)

    dist_df["vc_X_xconformist_vert_error"] = dist_df["vc"] * (dist_df["xconformist_vert_error"] - 0.5)
    dist_df["vp_X_phen_prestige_vert_error"] = dist_df["vp"] * (dist_df["phen_prestige_vert_error"] - 0.5)

    dist_df["hp_X_phen_prestige_horiz_error"] = dist_df["hp"] * (dist_df["phen_prestige_horiz_error"] - 0.5)
    dist_df["hc_X_xconformist_horiz_error"] = dist_df["hc"] * (dist_df["xconformist_horiz_error"] - 0.5)

    # only want to know variances of the following columns
    columns_only_var = ["phen_ad_logistic","phen_juv_logistic","phen_ad","phen_juv"]

    for column_only_var_i in columns_only_var:
        var_dict["var_" + column_only_var_i] = dist_df[column_only_var_i].var()

    # all columns for the variance covariance matrix
    columns_eta = [
            "aintercept"
            ,"bmat_phen_X_phen_mat_error"
            ,"bmat_envt_X_maternal_envt_cue_eror"
            ,"agen_X_g"
            ,"ajuv_X_cue_juv_envt_high"
            ,"vc_X_xconformist_vert_error"
            ,"vp_X_phen_prestige_vert_error"
            ,"hp_X_phen_prestige_horiz_error"
            ,"hc_X_xconformist_horiz_error"]

    formula_start = "phen_ad_logistic ~ "

    formula = formula_start

    not_first = False

    for column in columns_eta:
        if dist_df[column].var() > 0.0:

            if not_first: 
                formula += " + " + column
            else:
                formula += column
                not_first = True
        else:
            # variable excluded from eta calculations
            # just add a 0 to the var dict
            var_dict["eta2_" + column] = 0.0

    if formula != formula_start:
        
        # create the linear model
        lm_model = ols(formula,data=dist_df).fit()

        # get the anova table
        anova_table = sm.stats.anova_lm(lm_model,typ=2)

        # get the eta squares
        # this is classical eta^2, not partial eta^2 as it is 
        # SSR/SST rather than SSR/(SST + SSE)
        # (i.e., it sums up to 1)
        eta_sq = anova_table[:-1]["sum_sq"]/sum(anova_table["sum_sq"])

        # make a eta squared dict
        eta_sq_dict = eta_sq.to_dict()

        # add to var_dict, step-by-step, allowing us to change names
        for key, value in eta_sq_dict.items():
            var_dict["eta2_" + key] = value

    # now make sure all keys have the same order between files
    # as python chokes on this otherwise
    var_dict_keys = list(var_dict.keys())
    var_dict_keys.sort()

    for key in var_dict_keys:
        var_dict[key] = var_dict.pop(key)

    # return a data frame
    return(pd.DataFrame(data=var_dict,index=[1]))


summarizesims.SummarizeSims(
    path=args["path"]
    ,sep=args["sep"]
    ,pattern=args["pattern"]
    ,testing=False
    ,max_number_files=None if args["max_files"] == None else int(args["max_files"])
    ,posthoc_function=process_dist
    ,n_process=int(args["ncores"])
    )
