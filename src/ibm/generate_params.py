#!/usr/bin/env python3

# ptyhon script that 
# generates a batch file with all the different 
# simulations that need to run for all parameter
# combinations

import numpy as np
import socket
import datetime

# whether the survival curvive is 
# sigmoidal or not
sigmoidal_survival = [ 0 ]

envt_change_birth = [0]
juv_learns_remote = [1]

# frequency of the high environment
p = list(np.linspace(0,1,30))

survival_scalar_sig = [-2.5,3.5]
survival_scalar_quad = [0.8,0.0]

# combinations of maternal and juvenile cues
qjuv_mat_combinations = [[0.5,0.5]]

nloci_g = [ 3 ]

exe = "./xcue_integration.exe"

laplace = 1

nrep = 3

# for now we just need 10 zeros, which covers all the traits
#
#1.  g 
#2.  intercept 
#3.  ajuv 
#4.  agen 
#5.  bmat_phen
#6.  bmat_envt
#7.  hp
#8.  hc
#9.  vp
#10. vc
initval_list = [str(0.0) for x in range(0,10)]

initvals = " ".join(initval_list)

# ranges
aminmax = "-10.0 10.0"
gminmax = "-1.0 1.0"
bminmax = "-10.0 10.0"

#m = list(np.linspace(0, 1.0, 20))
m = [0.1]
#mu_g, mu_aintercept, mu_ajuv, mu_agen, mu_bmat_phen, mu_bmat_envt, mu_hp, mu_hc, mu_vp, mu_vc
#mu_combis = [[ 0.01 for x in range(0,10) ]]
zeros = [ 0 for i in range(0,10) ]

mu_combis = []

# mutation rate only genes 
mu_only_g = zeros[:]
mu_only_g[0] = 0.01
mu_combis.append(mu_only_g)

# mutation rate only genes and intercept
mu_only_ai = zeros[:]
mu_only_ai[1] = 0.01
mu_combis.append(mu_only_ai)

# mutation rate for everything except horiz and vert social learning
mu_no_social = [ 0.01 for i in range(0,6)]
mu_no_social = mu_no_social + [ 0 for i in range(0,4) ]

assert(len(mu_no_social) == 10)
mu_combis.append(mu_no_social)

mu_all = [ 0.01 for i in range(0,10) ]
mu_combis.append(mu_all)


# choose what consideration you want. For now only ai
mu_combis = [ mu_all ]
                        
#sd_h_noise_combs = [ [0.0,0.0],[0.05,0.05],[0.1,0.1],[0.2,0.2],[0.3,0.3],[0.4,0.4],[0.5,0.5] ]
sd_h_noise_combs = [[0.5,0.5]]
sd_vc_noise = [ 0.0 ]
sd_vp_noise = [ 0.0 ]
sd_mat_phen_noise = [ 0.0 ]

#mu_combis = [[ 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001 ]]
sdmu = "0.02 0.25 0.25"

# sampling sizes for social learning
# nph, nch (performance and conformity for horizontal)
# npv, ncv (performance and conformity for vertical)
nx = [[5,5,5,5]]

#nx = []
#nc_max = 5
#
## vary number of conformist models versus prestige models
#for nc_i in range(1,nc_max + 1):
#
#    np_i = nc_max + 1 - nc_i 
#    nx += [[np_i, nc_i, np_i, nc_i]]
#
## vary overall numbers of models
#for n_i in range(1, nc_max + 1):
#    nx += [[n_i, n_i, n_i, n_i ]]

# counter for the number of jobs
ctr = 1

# whether jobs should be run in the background
run_in_background = True

juvenile_survival = [ 0 ]
adult_survival = [ 1 ]

# never run background jobs on cluster
hostname = socket.gethostname()
if "athena" in hostname:
    run_in_background = False

# add ampersand to each job command if job needs
# to be run in background
bg = "&" if run_in_background else ""

date = datetime.datetime.now()
base_name = "sim_cue_integration_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

for rep_i in range(0,nrep):
    for sigmoidal_survival_i in sigmoidal_survival:

        if sigmoidal_survival_i == 0:
            survival_scalar_i = survival_scalar_quad
        else:
            survival_scalar_i = survival_scalar_sig

        survival_scalar_i_str = " ".join(str(x) for x in survival_scalar_i)

        for envt_change_birth_i in envt_change_birth:
            for p_i in p:
                p_i = round(p_i,3)
                for qjuv_mat_i in qjuv_mat_combinations:

                    qjuv_i = qjuv_mat_i[0]
                    qmat_i = qjuv_mat_i[1]

                    for nloci_g_i in nloci_g:
                        for sd_h_noise_i in sd_h_noise_combs:

                            sd_hc_noise_i = sd_h_noise_i[0]
                            sd_hp_noise_i = sd_h_noise_i[1]

                            for sd_vc_noise_i in sd_vc_noise:
                                for sd_vp_noise_i in sd_vp_noise:
                                    for sd_mat_phen_noise_i in sd_mat_phen_noise:
                                        for juv_learns_remote_i in juv_learns_remote:
                                            for m_i in m:
                                                m_i = round(m_i,3)
                                                for mu_combi_i in mu_combis:
                                                    
                                                    mu_combi_i_str = " ".join(
                                                            str(x) for x in mu_combi_i)

                                                    for nx_i in nx:
                                                        
                                                        nxstr = " ".join(
                                                                str(x) for x in nx_i)

                                                        for juvenile_survival_i in juvenile_survival:
                                                            for adult_survival_i in adult_survival:

                                                                print("echo " + str(ctr))


                                                                base_name_i = base_name + "_" + str(ctr)


                                                                print(exe + " \t"
                                                                        + str(sigmoidal_survival_i) + " \t"
                                                                        + str(laplace) + " "
                                                                        + str(p_i) + " \t"
                                                                        + survival_scalar_i_str + " \t"
                                                                        + str(qmat_i) + " "
                                                                        + str(qjuv_i) + " "
                                                                        + str(sd_hc_noise_i) + " "
                                                                        + str(sd_hp_noise_i) + " "
                                                                        + str(sd_vc_noise_i) + " "
                                                                        + str(sd_vp_noise_i) + " "
                                                                        + str(sd_mat_phen_noise_i) + " "
                                                                        + str(nloci_g_i) + " \t"
                                                                        + initvals + " \t"
                                                                        + gminmax + " "
                                                                        + aminmax + " "
                                                                        + bminmax + " \t"
                                                                        + mu_combi_i_str + " "
                                                                        + sdmu + " \t" 
                                                                        + str(m_i) + " "
                                                                        + str(nxstr) + " "
                                                                        + str(juvenile_survival_i) + " "
                                                                        + str(adult_survival_i) + " "
                                                                        + str(envt_change_birth_i) + " "
                                                                        + str(juv_learns_remote_i) + " "
                                                                        + base_name_i + " "
                                                                        + bg
                                                                        )

                                                                ctr += 1
