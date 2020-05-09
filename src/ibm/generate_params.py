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

# frequency of the high environment
p = list(np.linspace(0,1,50))
#p = [0.1]

survival_scalar_sig = [-2.5,3.5]
survival_scalar_quad = [0.8,0.0]

qmat = [0.5,1.0]
qjuv = [0.5,1.0]

# we go from qmat = 1.0 to 0.5
# while qjuv goes from 0.5 to 1.0
qmat = list(np.linspace(0.5,1.0,50))


nloci_g = [ 3 ]

exe = "./xcue_integration.exe"

laplace = 1

nrep = 1

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

sdmat = [ 0.05 ]
sdsoc_horiz = [ 0.5, 1.0 ]
sdsoc_vert = [ 0.05 ]

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
#mu_combis = [ mu_no_social, mu_all ]
mu_combis = [ mu_all ]
                        
max_error_conform_horiz = [ 0.0 ]
max_error_prestige_horiz = [ 0.0 ]
max_error_conform_vert = [ 0.0 ]
max_error_prestige_vert = [ 0.0 ]
max_error_mat_phen = [ 0.0 ]
max_error_mat_envt = [ 0.0 ]

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
run_in_background = False

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
                for qmat_i in qmat:
                    qjuv = [1.5 - qmat_i]
                    for qjuv_i in qjuv:
                        for nloci_g_i in nloci_g:
                            for max_error_conform_horiz_i in max_error_conform_horiz:
                                for max_error_prestige_horiz_i in max_error_prestige_horiz:
                                    for max_error_conform_vert_i in max_error_conform_vert:
                                        for max_error_prestige_vert_i in max_error_prestige_vert:
                                            for max_error_mat_phen_i in max_error_mat_phen:
                                                for max_error_mat_envt_i in max_error_mat_envt:
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
                                                                        ctr += 1


                                                                        base_name_i = base_name + "_" + str(ctr)


                                                                        print(exe + " \t"
                                                                                + str(sigmoidal_survival_i) + " \t"
                                                                                + str(laplace) + " "
                                                                                + str(p_i) + " \t"
                                                                                + survival_scalar_i_str + " \t"
                                                                                + str(qmat_i) + " "
                                                                                + str(qjuv_i) + " "
                                                                                + str(max_error_conform_horiz_i) + " "
                                                                                + str(max_error_prestige_horiz_i) + " "
                                                                                + str(max_error_conform_vert_i) + " "
                                                                                + str(max_error_prestige_vert_i) + " "
                                                                                + str(max_error_mat_phen_i) + " "
                                                                                + str(max_error_mat_envt_i) + " "
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
                                                                                + base_name_i + " "
                                                                                + bg
                                                                                )
