#!/usr/bin/env python3

import os.path
import pandas as pd
import sys

def find_max_eta(row, traits_n_labels, lab=""):

    # get the names we need to sort on
    keys = list(traits_n_labels.keys())
    
    # find those columns in the row
    series_eta = row.loc[keys]

    # sort them on magnitude
    series_eta_sorted = series_eta.sort_values(ascending=False)

    # get the relevant keys, corresponding to the max
    list_keys = list(series_eta_sorted.iloc[[0,1]].index[[0,1]])

    # also sort the second one
    delta_eta2 = series_eta_sorted.iloc[0] - series_eta_sorted.iloc[1]

#    # find the maximum eta value
#    max_eta_idx = row[keys].astype("float64").idxmax()
#
#    # second biggest key 
#    max2_eta_idx = row[keys].astype("float64").idxmax()



    return(pd.Series({"eta2_" + lab + "max" : keys.index(list_keys[0]), "delta_eta2_" + lab + "max" : delta_eta2,"eta2_" + lab + "2nd": keys.index(list_keys[1]) }))


#data_file_name = "summary_merged.csv"
#data_file_name = "summary_contour_vary_sdsoc_vert.csv"



#data_dir = os.path.join(os.path.expanduser("~"),"Projects/cue_integration_culture/data/")
## read in the data
#data = pd.read_csv(os.path.join(data_dir,data_file_name)
#        ,sep=";")

filename = sys.argv[1]

data = pd.read_csv(filename,sep=";")

traits_n_labels = {
        "eta2_aintercept":r"$a_{\text{0}}$",
        "eta2_agen_X_g":r"$a_{\mathrm{g}}$",
        "eta2_ajuv_X_cue_juv_envt_high":r"$a_{\mathrm{juv}}$",
        "eta2_bmat_envt_X_maternal_envt_cue_eror":r"$m_{\mathrm{e}}$",
        "eta2_bmat_phen_X_phen_mat_error":r"$m_{\mathrm{m}}$",
        "eta2_vc_X_xconformist_vert_error":r"$v_{\mathrm{c}}$",
        "eta2_vp_X_phen_prestige_vert_error":r"$v_{\mathrm{p}}$",
        "eta2_hc_X_xconformist_horiz_error":r"$h_{\mathrm{c}}$",
        "eta2_hp_X_phen_prestige_horiz_error":r"$h_{\mathrm{p}}$",
        }

traits_n_labels_juv = {
        "eta2_jv_aintercept":r"$a_{\text{0}}$",
        "eta2_jv_agen_X_g":r"$a_{\mathrm{g}}$",
        "eta2_jv_ajuv_X_cue_juv_envt_high":r"$a_{\mathrm{juv}}$",
        "eta2_jv_bmat_envt_X_maternal_envt_cue_eror":r"$m_{\mathrm{e}}$",
        "eta2_jv_bmat_phen_X_phen_mat_error":r"$m_{\mathrm{m}}$",
        "eta2_jv_vc_X_xconformist_vert_error":r"$v_{\mathrm{c}}$",
        "eta2_jv_vp_X_phen_prestige_vert_error":r"$v_{\mathrm{p}}$",
        }

data[["eta2_max","delta_eta2","eta2_2nd"]] = data.apply(func=find_max_eta
            ,axis=1 # apply to each row
            ,args=(traits_n_labels,)
            )

data[["eta2_jv_max","delta_jv_eta2","eta2_jv_2nd"]] = data.apply(func=find_max_eta
            ,axis=1 # apply to each row
            ,args=(traits_n_labels_juv,"_jv_")
            )

file_path, ext = os.path.splitext(filename)

data.to_csv(file_path + "_max_eta_" + ext,sep=";")
