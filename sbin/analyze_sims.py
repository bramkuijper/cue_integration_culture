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


parser.add_argument("--sep"
                    ,help="Column separator"
                    ,default=";")

args = vars(parser.parse_args())

def process_dist(filename):

    splitname = os.path.splitext(filename)
    dist_filename = splitname[0] + "_dist" + splitname[1]

    dist_df = pd.read_csv(dist_filename, sep=";")

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
    ,max_number_files=5
    ,posthoc_function=process_dist
    )
