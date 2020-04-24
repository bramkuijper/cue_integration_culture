#!/usr/bin/env python3

# summarize the individual-based model
import summarizesims
import argparse
import pandas as pd
import os.path
import re



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

    splitname = os.path.splitext(filename)
    dist_filename = splitname[0] + "_dist" + splitname[1]

    dist_df = pd.read_csv(dist_filename, sep=";")

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

    names = dist_df.columns.values
#
    # find whether there are any 'Unnamed...' columns
    unnamed_cols = [ i for i in names if re.search("Unnamed",i) != None ]
#

#    columns_covmat = list(set(names) - set(["id","patch_id"]) - set(unnamed_cols))
    
    var_dict = {}

    # only want to know variances of the following columns
    columns_only_var = ["phen_ad_logistic","phen_juv_logistic","phen_ad","phen_juv"]

    for column_only_var_i in columns_only_var:
        var_dict["var_" + column_only_var_i] = dist_df[column_only_var_i].var()

    # all columns for the variance covariance matrix
    columns_covmat = [
            "aintercept"
            ,"bmat_phen_X_phen_mat_error"
            ,"bmat_envt_X_maternal_envt_cue_eror"
            ,"agen_X_g"
            ,"ajuv_X_cue_juv_envt_high"
            ,"vc_X_xconformist_vert_error"
            ,"vp_X_phen_prestige_vert_error"
            ,"hp_X_phen_prestige_horiz_error"
            ,"hc_X_xconformist_horiz_error"]


#    # variance covariance matrix
    cov_mat = dist_df[columns_covmat].cov()


    dict_cov_names = {}

    # get the column_names for the entries along the 
    # upper diagonal. To prevent names occurring twice
    # we first generate the names combs, sort them 
    # and write them to a dict
    for name_i in columns_covmat:
        for name_j in columns_covmat:
            comb_list = [name_i,name_j]
            dict_cov_names["".join(comb_list)] = comb_list

    var_total = dist_df["phen_ad_logistic"].var()
    dist_df["phen_ad_logistic2"] = 0

    # calculate total variance minus variance in one element #####
    # I do this, as sometimes covariances turn out to be negative 
    # and then the whole variance is negative, so we must be missing
    # something here
    for name_i in columns_covmat:

        # generate column name
        colname = "total_ad_minus_" + name_i

        # make a column in the original data frame where we 
        # subtract the column in question from the total sum
        dist_df[colname] = dist_df["phen_ad_logistic"] - dist_df[name_i]

        # calculate a varianace of this column
        var_dict["var_total_ad_minus_" + name_i] = dist_df[colname].var()
        var_dict["var_prop_ad_minus_" + name_i] = (var_total - var_dict["var_total_ad_minus_" + name_i])/var_total

        dist_df["phen_ad_logistic2"] += dist_df[name_i]

    var_dict["var_phen_ad_logistic2"] = dist_df["phen_ad_logistic2"].var()

    # get the data out of the covariance matrix and into the final dataframe
    for key, value in dict_cov_names.items():

        # give names to the cov/var columns
        entry_name = "cov_" + value[0] + "_" + value[1] if value[0] != value[1] else "var_" + value[0]
           
        # store the covariance in the dictionary
        var_dict[entry_name] = cov_mat.loc[value[0],value[1]]

    # return a data frame
    return(pd.DataFrame(data=var_dict,index=[1]))


summarizesims.SummarizeSims(
    path=args["path"]
    ,sep=args["sep"]
    ,pattern=args["pattern"]
    ,testing=False
    ,max_number_files=None if args["max_files"] == None else int(args["max_files"])
    ,posthoc_function=process_dist
    )
