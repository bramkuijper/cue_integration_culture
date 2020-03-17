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

    # add some additional columns
    dist_df["bmat_phen_X_phen_mat_error"] = dist_df["bmat_phen"] * dist_df["phen_mat_error"]
    dist_df["bmat_envt_X_maternal_envt_cue_erorr"] = \
            dist_df["bmat_envt"] * dist_df["maternal_envt_cue_error"]

    dist_df["agen_X_g"] = dist_df["agen"] * dist_df["g"]
    dist_df["ajuv_X_cue_juv_envt_high"] = dist_df["ajuv"] * dist_df["cue_juv_envt_high"]

    dist_df["vc_X_xconformist_vert_error"] = dist_df["vc"] * dist_df["xconformist_vert_error"]
    dist_df["vp_X_phen_prestige_vert_error"] = dist_df["vp"] * dist_df["phen_prestige_vert_error"]

    dist_df["hp_X_phen_prestige_horiz_error"] = dist_df["hp"] * dist_df["phen_prestige_horiz_error"]
    dist_df["hc_X_xconformist_horiz_error"] = dist_df["hc"] * dist_df["xconformist_horiz_error"]

    names = dist_df.columns.values

    # find whether there are any 'Unnamed...' columns
    unnamed_cols = [ i for i in names if re.search("Unnamed",i) != None ]

    names_no_id = list(set(names) - set(["id","patch_id"]) - set(unnamed_cols))

    # variance covariance matrix
    cov_mat = dist_df[names_no_id].cov()

    var_dict = {}

    dict_cov_names = {}

    # get the column_names for the entries along the 
    # upper diagonal. To prevent names occurring twice
    # we first generate the names combs, sort them 
    # and write them to a dict
    for name_i in names_no_id:
        for name_j in names_no_id:
            comb_list = [name_i,name_j]
            comb_list.sort()

            dict_cov_names["".join(comb_list)] = comb_list

    for key, value in dict_cov_names.items():

        entry_name = "cov_" + value[0] + "_" + value[1] if value[0] != value[1] else "var_" + value[0]
            
        var_dict[entry_name] = cov_mat.loc[value[0],value[1]]
            
    return(pd.DataFrame(data=var_dict,index=[1]))


summarizesims.SummarizeSims(
    path=args["path"]
    ,sep=args["sep"]
    ,pattern=args["pattern"]
    ,testing=False
    ,max_number_files=None if args["max_files"] == None else int(args["max_files"])
    ,posthoc_function=process_dist
    )
