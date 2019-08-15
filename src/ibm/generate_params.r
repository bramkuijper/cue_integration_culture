
# r code to produce different permutations of 
# possible parameters
# to run the simulations

# whether the survival curve is simply decelerating
# or decelerating in a sigmoidal fashion (see Leimar &
# McNamara 2015)
sigmoidal_survival = c(0)

# frequency of the high 
p = c(0.9)

# survival scalar values when 
survival_scalar_sig = c(-2.5,3.5)
survival_scalar_quad = c(0.8, 0.0)
qmat = c(0.5,0.9)
qjuv = c(0.5,0.9)

# whether one uses laplace type of mutations or not
# yes
laplace = 1

nloci_g = c(3)

sdmat = 0.05

m = c(0.1,0.5)

# initial values for all the evolving traits
# minus the genetic cue loci
initvals = paste(rep("0.0 ",times=6),collapse=" ")

# mutation rates for all the evolving traits
mu_combis = c(
        paste(rep("0.0001",times=6),collapse=" ")
                )

# mutation standard devations for a, b and g
sdmuvals = " 0.02 0.25 0.25 "

gminmax = "-1.0 1.0"
aminmax = "0.0 8.0"
bminmax = "-10.0 10.0"

# sure, I could use something like expand.grid to get
# all parameter permutations. However, in case you get constraints
# where certain combinations only occur conditional upon others,
# then this does not work anymore
# Hence, basic for loops...

batch_str = ""

# the name of the exe file
# change this at will
executable_name = "./xcue_integration"

ctr = 0

# number of replicate simulations for each parameter permutation
nrep = 1

# loop through all possible combinations
for (rep_i in 1:nrep)
{
    for (sigmoidal_survival_i in sigmoidal_survival)
    {
        survival_scalar_i = ""

        if (sigmoidal_survival_i)
        {
            survival_scalar_i = paste(survival_scalar_sig, collapse=" ")
        } else {
            survival_scalar_i = paste(survival_scalar_sig, collapse=" ")
        }

        for (p_i in p)
        {
            for (qmat_i in qmat)
            {
                for (qjuv_i in qjuv)
                {
                    for (nloci_g_i in nloci_g)
                    {
                        for (sdmat_i in sdmat)
                        {
                            for (mu_combi_i in mu_combis)
                            {
                                for (m_i in m)
                                {
                                    # generate the command calling the exe with
                                    # all parameters
                                    str_command = paste(
                                            executable_name
                                            ,sigmoidal_survival_i
                                            ,laplace
                                            ,p_i
                                            ,survival_scalar_i
                                            ,qmat_i
                                            ,qjuv_i
                                            ,nloci_g_i
                                            ,initvals
                                            ,gminmax
                                            ,aminmax
                                            ,bminmax
                                            ,sdmat_i
                                            ,mu_combi_i
                                            ,sdmuvals
                                            ,m_i
                                            ,"\n"
                                            ,sep=" ")

                                    # generate a string saying echo 
                                    # so that DOS prints a number for each run
                                    echo_str = paste(
                                            "echo "
                                            ,ctr
                                            ,"\n")

                                    ctr = ctr + 1

                                    # add it to master string
                                    # which will be written to file later
                                    batch_str = paste(
                                            batch_str
                                            ,echo_str
                                            ,str_command)
                                }
                            }
                        }
                    }
                }
            } # end for (qjuv_i in qjuv)
        } # end for (p_i in p)
    } # for (sigmoidal_survival_i in sigmoidal_survival)
} # for (rep_i in 1:nrep)


filename = "parameter_batch_file.bat"
write(x=batch_str, filename)
