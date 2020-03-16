#!/usr/bin/env python3

# summarize the individual-based model
import summarizesims
import argparse
import pandas as pd


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

    names_no_id = set(names) - set(["id","patch_id"])

    # variance covariance matrix
    cov_mat = dist_df[names_no_id].cov()

    var_dict = {}

    # now write diagonial to dataframe
    for names_no_id as name_i:
        for names_no_id as name_j:
            if name_i == name_j:
                var_dict["var_" + name_i + "_" + name_j] = cov_mat.loc[name_i,name_j]
            else:
                var_dict["cov_" + name_i + "_" + name_j] = cov_mat.loc[name_i,name_j]
            
    return(pd.Dataframe(data=var_dict,index=[1]))


summarizesims.SummarizeSims(
    path=args["path"]
    ,sep=args["sep"]
    ,pattern=args["pattern"]
    ,posthoc_function=process_dist
    )
