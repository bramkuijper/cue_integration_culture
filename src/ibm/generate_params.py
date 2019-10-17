#!/usr/bin/env python3

# ptyhon script that 
# generates a batch file with all the different 
# simulations that need to run for all parameter
# combinations

import numpy as np
import socket

# whether the survival curvive is 
# sigmoidal or not
sigmoidal_survival = [ 0 ]

# frequency of the high environment
#p = list(np.linspace(0,1,11))
p = list(np.linspace(0,1,5))
survival_scalar_sig = [-2.5,3.5]
survival_scalar_quad = [0.8,0.0]
qmat = [ 0.5, 0.9, 1.0 ]
qjuv = [ 0.5, 0.9, 1.0 ]

nloci_g = [ 3 ]

exe = "./xcue_integration"

laplace = 1

nrep = 5

# for now we just need 12 zeros, which covers all the traits
initvals = " ".join([str(0.0) for x in range(0,12)])

aminmax = "0.0 8.0"
gminmax = "-1.0 1.0"
bminmax = "-10.0 10.0"

sdmat = [ 0.05 ]
sdsoc_horiz = [ 1.0, 0.1, 0.05 ]
sdsoc_vert = [ 1.0, 0.1, 0.05 ]
#m = list(np.linspace(0, 1.0, 11))
m = [ 0.1]

# mu_g, mu_amat, mu_ajuv, mu_agen, mu_asoc_horiz, mu_asoc_vert, mu_bmat_phen, mu_bmat_envt, mu_hp, mu_hc, mu_vp, mu_vc
mu_combis = [ [ 0.01 for x in range(0,12) ]]
#mu_combis = [[ 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001 ]]
sdmu = "0.02 0.25 0.25"

# sampling sizes for social learning
# nph, nch (performance and conformity for horizontal)
# npv, ncv (performance and conformity for vertical)
nx = [[5,5,5,5]]

# counter for the number of jobs
ctr = 1

# whether jobs should be run in the background
run_in_background = False

# never run background jobs on cluster
hostname = socket.gethostname()
if "carson" in hostname:
    run_in_background = False

# add ampersand to each job command if job needs
# to be run in background
bg = "&" if run_in_background else ""

for rep_i in range(0,nrep):
    for sigmoidal_survival_i in sigmoidal_survival:

        if sigmoidal_survival_i == 0:
            survival_scalar_i = survival_scalar_quad
        else:
            survival_scalar_i = survival_scalar_sig

        survival_scalar_i_str = " ".join(str(x) for x in survival_scalar_i)

        for p_i in p:
            p_i = round(p_i,3)
            for qmat_i in qmat:
                for qjuv_i in qjuv:
                    for nloci_g_i in nloci_g:
                        for sdmat_i in sdmat:
                            for sdsoc_vert_i in sdsoc_vert:
                                for sdsoc_horiz_i in sdsoc_horiz:
                                    for m_i in m:
                                        m_i = round(m_i,3)
                                        for mu_combi_i in mu_combis:
                                            
                                            mu_combi_i_str = " ".join(
                                                    str(x) for x in mu_combi_i)

                                            for nx_i in nx:

                                                nxstr = " ".join(
                                                        str(x) for x in nx_i)

                                                print("echo " + str(ctr))

                                                ctr += 1

                                                print(exe + " "
                                                        + str(sigmoidal_survival_i) + " "
                                                        + str(laplace) + " "
                                                        + str(p_i) + " "
                                                        + survival_scalar_i_str + " "
                                                        + str(qmat_i) + " "
                                                        + str(qjuv_i) + " "
                                                        + str(nloci_g_i) + " "
                                                        + initvals + " "
                                                        + gminmax + " "
                                                        + aminmax + " "
                                                        + bminmax + " "
                                                        + str(sdmat_i) + " "
                                                        + str(sdsoc_vert_i) + " "
                                                        + str(sdsoc_horiz_i) + " "
                                                        + mu_combi_i_str + " "
                                                        + sdmu + " " 
                                                        + str(m_i) + " "
                                                        + str(m_i) + " "
                                                        + str(nxstr) + " "
                                                        + bg
                                                        )
