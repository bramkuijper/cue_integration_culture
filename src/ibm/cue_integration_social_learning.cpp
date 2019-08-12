// Extending Leimar & McNamara's cue integration model
// to include cultural transmission
// Bram Kuijper 
// 2019
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>


// various functions, such as unique filename creation
#include "auxiliary.hpp"
#include "individual.hpp"

#define DEBUG

// standard namespace
using namespace std;

int seed = get_nanoseconds();
mt19937 rng_r{static_cast<long unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);

// parameters & variables:

// number of individuals in population
const int NPatches = 400;
const int NBreeder = 100;
const int NClutch = 10;

// number of generations
int number_generations = 50000;

// environmental switch rate
//
// parameters below will be changed on the command line

// frequency of high state patches
double p = 0.0;

// quality (beyond 0.5 error rate) of juvenile cue
double qjuv = 0.5;

// quality (beyond 0.5 error rate) of maternal cue
double qmat = 0.5;

// whether survival selection has a sigmoidal or a quadratic function
bool sigmoidal_survival = true;

// number of loci underlying the genetic cue
int nloci_g = 3;

// initial values
double init_g = 0.0;
double init_amat = 0.0;
double init_ajuv = 0.0;
double init_agen = 0.0;
double init_bmat_phen = 0.0;
double init_bmat_envt = 0.0;

// ranges for traits
double gmin = 0.0;
double gmax = 0.0;
double amin = 0.0;
double amax = 0.0;
double bmin = 0.0;
double bmax = 0.0;

// mutation rates
double mu_g = 0.0;
double mu_amat = 0.0;
double mu_ajuv = 0.0;
double mu_agen = 0.0;
double mu_bmat_phen = 0.0;
double mu_bmat_envt = 0.0;

double sdmu_a = 0.0;
double sdmu_b = 0.0;
double sdmu_g = 0.0;

//write out the data every nth generation
int data_nth_generation = 10;

struct Patch
{
    Individual breeders[NBreeder];
    int n_breeders; // number of breeders currently in patch

    vector <Individual> juveniles;

    bool envt_high;  // high-state patch yes/no
};

Patch Pop[NPatches];

// survival function
double survival_probability(double const ad_phen, bool const state_high)
{
    if (sigmoidal_survival)
    {
        double scalar = state_high ? 3.5 : -2.5;

        return(1.0 / (1.0 + exp(scalar + 6.0 * ad_phen)));
    }

    return(state_high ? 
            1.0 - 0.8 * (1.0 - ad_phen) * (1.0 - ad_phen)
            :
            1.0 - 0.8 * ad_phen * ad_phen
            );
}


// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{
    sigmoidal_survival = atoi(argv[1]);
    p = atof(argv[2]);
    qmat = atof(argv[2]);
    qjuv = atof(argv[2]);
    n_loci_g = atoi(argv[2]);
    init_g = atof(argv[2]);
    init_amat = atof(argv[2]);
    init_ajuv = atof(argv[2]);
    init_agen = atof(argv[2]);
    init_bmat_phen = atof(argv[2]);
    init_bmat_envt = atof(argv[2]);
    gmin = atof(argv[2]);
    gmax = atof(argv[2]);
    amin = atof(argv[2]);
    amax = atof(argv[2]);
    bmin = atof(argv[2]);
    bmax = atof(argv[2]);
    sdmat = atof(argv[2]);

    mu_g = atof(argv[2]);
    sdmu_g = atof(argv[2]);
    mu_amat = atof(argv[2]);
    mu_ajuv = atof(argv[2]);
    mu_agen = atof(argv[2]);
    mu_bmat_phen = atof(argv[2]);
    mu_bmat_envt = atof(argv[2]);
    sdmu_a = atof(argv[2]);
    sdmu_b = atof(argv[2]);

}

// write down all parameters in the file
void write_parameters(ofstream &DataFile)
{
    DataFile << endl << endl
        << "sigmoidal_survival;" << sigmoidal_survival << ";"
        << "p;" << p << ";"
        << "qmat;" << qmat << ";"
        << "qjuv;" << qjuv << ";"
        << "n_loci_g;" << n_loci_g << ";"
        << "init_g;" << init_g << ";"
        << "init_amat;" << init_amat << ";"
        << "init_ajuv;" << init_ajuv << ";"
        << "init_agen;" << init_agen << ";"
        << "init_bmat_phen;" << init_bmat_phen << ";"
        << "init_bmat_envt;" << init_bmat_envt << ";"
        << "gmin;" << gmin << ";"
        << "gmax;" << gmax << ";"
        << "amin;" << amin << ";"
        << "amax;" << amax << ";"
        << "bmin;" << bmin << ";"
        << "bmax;" << bmax << ";"
        << "sdmat;" << sdmat << ";"
        << "mu_g;" << mu_g << ";"
        << "sdmu_g;" << sdmu_g << ";"
        << "mu_amat;" << mu_amat << ";"
        << "mu_ajuv;" << mu_ajuv << ";"
        << "mu_agen;" << mu_agen << ";"
        << "mu_bmat_phen;" << mu_bmat_phen << ";"
        << "mu_bmat_envt;" << mu_bmat_envt << ";"
        << "sdmu_a;" << sdmu_a << ";"
        << "sdmu_b;" << sdmu_b << ";"
        << "seed;" << seed << ";" << endl;
}

// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile << "generation;time;mean_theta_a;mean_theta_b;mean_phi_a;mean_phi_b;mean_resources;var_theta_a;var_theta_b;var_phi_a;var_phi_b;var_resources;nwinter;nstaging;nsummer;nkids;mean_flock_size_summer;mean_flock_size_winter;mean_staging_size_winter;mean_staging_size_summer;" << endl;
}

// write data both for winter and summer populations
void write_stats(ofstream &DataFile, int generation, int timestep)
{
}

// initialize the population at the start of the simulation
void init_population()
{
    // auxiliary variable whether mothers perceived a cue
    // that the environment is in a high state
    bool cue_juv_envt_high;

    // loop through all individuals 
    // and assign them values for the cue loci
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // patch in a high state or not
        Pop[patch_i].envt_high = uniform(rng_r) < p;
       
        // cue given to adults
        cue_juv_envt_high = uniform(rng_r) < qmat ?
            Pop[patch_i].envt_high 
            :
            !Pop[patch_i].envt_high;
        
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            for (int allele_i = 0; allele_i < 2; ++allele_i)
            {
                for (int g_loc_i = 0; g_loc_i < nloci_g; ++g_loc_i)
                {
                    // genetic cue values
                    Pop[patch_i].breeders[breeder_i].g[g_loc_i][allele_i] = init_g;
                }

                // maternal cue weighting
                Pop[patch_i].breeders[breeder_i].amat[allele_i] = init_amat;

                // juvenile cue weighting
                Pop[patch_i].breeders[breeder_i].ajuv[allele_i] = init_ajuv;

                // genetic cue weighting
                Pop[patch_i].breeders[breeder_i].agen[allele_i] = init_agen;
                
                // maternal phenotypic cue weighting
                Pop[patch_i].breeders[breeder_i].bmat_phen[allele_i] = init_bmat_phen;
                
                // maternal phenotypic cue weighting
                Pop[patch_i].breeders[breeder_i].bmat_envt[allele_i] = init_bmat_envt;
            } //end for allele_i
            
            Pop[patch_i].n_breeders = NBreeder;
            
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_juv_envt_high;

            Pop[patch_i].breeders[breeder_i].cue_juv_envt_high = 
                cue_juv_envt_high;

            Pop[patch_i].breeders[breeder_i].u = 1.0 /
                (1.0 + exp(-init_agen * init_g - ajuv * cue_juv_envt_high));
        } // end for breeder_i
    } // end for patch_i
} // end void init_population()

// mutation of a certain allele with value val
// given mutation rate mu and mutational distribution stdev sdmu
double mutation(double val, double mu, double sdmu)
{
    if (gsl_rng_uniform(rng_r) < mu)
    {
        val += gsl_ran_gaussian(rng_r, sdmu);
    }

    return(val);
}

// auxiliary function to bound values between [min,max]
void clamp(double &val, double min, double max)
{
    val = val > max ? max : val < min ? min : val;
}

// create a new offspring
void create_offspring(Individual &mother
        ,Individual &father
        ,Individual &offspring
        ,bool offspring_envt)
{
    // set up a bernoulli distribution that returns 0s or 1s
    // at equal probability to sample alleles from the first or
    // second set of chromosomes of diploid individuals
    bernoulli_distribution allele_sample(0.5);

    double sum_genes = 0.0;

    for (int g_loc_i = 0; g_loc_i < nloci_g; ++g_loc_i)
    {
        // genetic cue values
        Kid.g[g_loc_i][0] = mutation(
                Mother.g[g_loc_i][allele_sample(rng_r)],
                mu_g,
                sdmu_g);

        clamp(Kid.g[g_loc_i][0], gmin, gmax);
        
        Kid.g[g_loc_i][1] = mutation(
                Father.g[g_loc_i][allele_sample(rng_r)],
                mu_g,
                sdmu_g);

        clamp(Kid.g[g_loc_i][1], gmin, gmax);

        sum_genes += 0.5 * (Kid.g[g_loc_i][0] + Kid.g[g_loc_i][1]);
    }

    // inheritance of maternal cue values 
    Kid.amat[0] = mutation(
            Mother.amat[allele_sample(rng_r)],
            mu_amat,
            sdmu_a);

    clamp(Kid.amat[0], amin, amax);

    Kid.amat[1] = mutation(
            Father.amat[allele_sample(rng_r)],
            mu_amat,
            sdmu_a);

    clamp(Kid.amat[1], amin, amax);

    double amat_phen = 0.5 * (Kid.amat[0] + Kid.amat[1]);
   

    // inheritance of juvenile cue values 
    Kid.ajuv[0] = mutation(
            Mother.ajuv[allele_sample(rng_r)],
            mu_ajuv,
            sdmu_a);

    clamp(Kid.ajuv[0], amin, amax);

    Kid.ajuv[1] = mutation(
            Father.ajuv[allele_sample(rng_r)],
            mu_ajuv,
            sdmu_a);

    clamp(Kid.ajuv[1], amin, amax);

    double ajuv_phen = 0.5 * (Kid.ajuv[0] + Kid.ajuv[1]);
   


    // inheritance of genetic cue values 
    Kid.agen[0] = mutation(
            Mother.agen[allele_sample(rng_r)],
            mu_agen,
            sdmu_a);

    clamp(Kid.agen[0], amin, amax);

    Kid.agen[1] = mutation(
            Father.agen[allele_sample(rng_r)],
            mu_agen,
            sdmu_a);

    clamp(Kid.agen[1], amin, amax);

    double agen_phen = 0.5 * (Kid.agen[0] + Kid.agen[1]);


    // inheritance of maternal phenotypic cue values 
    Kid.bmat_phen[0] = mutation(
            Mother.bmat_phen[allele_sample(rng_r)],
            mu_bmat_phen,
            sdmu_b);

    clamp(Kid.bmat_phen[0], bmin, bmax);

    Kid.bmat_phen[1] = mutation(
            Father.bmat_phen[allele_sample(rng_r)],
            mu_bmat_phen,
            sdmu_b);

    clamp(Kid.bmat_phen[1], bmin, bmax);

    double bmat_phen_phen = 0.5 * (Kid.bmat_phen[0] + Kid.bmat_phen[1]);
    
    // inheritance of maternal environmental cue values 
    Kid.bmat_envt[0] = mutation(
            Mother.bmat_envt[allele_sample(rng_r)],
            mu_bmat_envt,
            sdmu_b);

    clamp(Kid.bmat_envt[0], bmin, bmax);

    Kid.bmat_envt[1] = mutation(
            Father.bmat_envt[allele_sample(rng_r)],
            mu_bmat_envt,
            sdmu_b);

    clamp(Kid.bmat_envt[1], bmin, bmax);

    double bmat_envt_phen = 0.5 * (Kid.bmat_envt[0] + Kid.bmat_envt[1]);


    // kid receives juvenile cue
    Kid.cue_juv_envt_high = uniform(rng_r) < qjuv ? 
        offspring_envt : !offspring_envt;

    // kid receives adult cue
    Kid.cue_ad_envt_high = uniform(rng_r) < qmat ? 
        offspring_envt : !offspring_envt;

    double dmat_weighting = Kid.cue_ad_envt_high ? 1.0 : -1.0;

    // generate maternal cue
    Kid.xoff = 1.0 /
        (1.0 + exp(
                   -Kid.bmat_phen_phen * (Mother.ad_phen - 0.5) + 
                   dmat_weighting * Kid.bmat_envt_phen));

    // noise in the maternal cue
    normal_distribution<> maternal_noise(0,sdmat);

    xmat = Mother.xoff + maternal_noise(rng_r);

    Kid.ad_phen = 1.0 / 
        (1.0 + exp(
                   -amat_phen * Mother.xmat +
                   -agen_phen * sum_genes +
                   -ajuv_phen * Kid.cue_juv_envt_high));

} // end create_offspring()

void survive_reproduce()
{
    // set up a random number generator to 
    // sample from the remaining breeders
    uniform_int_distribution<> random_patch(
            0,
            NPatches - 1);


    // how many offspring need to be produced in each patch
    int n_offspring_per_patch = Clutch * NBreeders;
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        assert(Pop[patch_i].n_breeders > 0);
        assert(Pop[patch_i].n_breeders <= NBreeder);

        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            // individual dies if random uniform number is larger 
            // than the survival probability
            if (uniform(rng_r) > 
                    survival_probability(
                        Pop[patch_i].breeders[breeder_i].u
                        , Pop[patch_i].envt_high)
            )
            {
                // delete individual
                Pop[patch_i].breeders[breeder_i] = Pop[patch_i].breeders[Nbreeder - 1];
                --breeder_i;
                --Pop[patch_i].n_breeders;
            }
        } // end for (int breeder_i
        
        if (Pop[patch_i].n_breeders > 0)
        {
            // set up a random number generator to 
            // sample from the remaining breeders
            uniform_int_distribution<> random_breeder(
                    0,
                    Pop[patch_i].n_breeders - 1);

            // let breeders produce offspring
            for (int offspring_i = 0;
                    offspring_i < n_offspring_per_patch;
                    ++offspring_i)
            {
                Individual Kid;

                create_offspring(
                        Pop[patch_i].breeders[random_breeder(rng_r)]
                        ,Pop[patch_i].breeders[random_breeder(rng_r)]
                        ,Kid
                        ,Pop[patch_i].envt_high
                );

                // dispersal or not
                if (uniform(rng_r) > m)
                {
                    Pop[patch_i].juveniles.push_back(Kid);
                }
                else
                {
                    Pop[random_patch(rng_r)].juveniles.push_back(Kid);
                }
            } // for (int offspring_i = 0;
        } // end if (Pop[patch_i].n_breeders > 0)

        // change the envt
        if (p < uniform(rng_r))
        {
            Pop[patch_i].envt_high = !Pop[patch_i].envt_high;
        }
    } // end for int patch_i
} // end survive_reproduce()

void replace()
{
    int rand_juv;

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {

        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {

            assert(Pop[patch_i].juveniles.size() >= NBreeder);

            // set up a random number generator to 
            // sample from the remaining breeders
            uniform_int_distribution<> random_juveniles(
                    0,
                    Pop[patch_i].juveniles.size() - 1);

            rand_juv = random_juveniles(rng_r);

            Pop[patch_i].breeders[breeder_i] = 
                Pop[patch_i].juveniles[rand_juv];

            Pop[patch_i].juveniles.erase(
                    Pop[patch_i].juveniles.begin() 
                    + rand_juv);

        } // for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)

        Pop[patch_i].n_breeders = NBreeder;
        
        // clear out all the juveniles
        Pop[patch_i].juveniles.clear();

    } // end for (int patch_i = 0
}

// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    string filename = "sim_cue_integration";
    create_filename(filename);
    ofstream DataFile(filename.c_str());  // output file 

    // get command line arguments
    init_arguments(argc, argv);

    // write headers to the datafile
    write_data_headers(DataFile);

    // initialize the population
    init_population();

    for (int generation = 0; generation < number_generations; ++generation)
    {
        survive_reproduce();

        replace();

        if (generation % data_nth_generation == 0)
        {
            write_stats(DataFile, generation, 2);
        }
    }

    write_parameters(DataFile);
}
