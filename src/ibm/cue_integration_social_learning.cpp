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

// C++ random number generation
int seed = get_nanoseconds();
mt19937 rng_r{static_cast<long unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);

// parameters & variables:

// number of individuals in population
const int NPatches = 400;
const int NBreeder = 100;

// number of generations
int number_generations = 50000;

// environmental switch rate
//
// parameters below will be changed on the command line

// frequency of high state patches
double p = 0.0;

// survival function in low vs high patches
double survival_scalar[2] = {0.0, 0.0};

// quality (beyond 0.5 error rate) of juvenile cue
double qjuv = 0.5;

// quality (beyond 0.5 error rate) of maternal cue
double qmat = 0.5;

// noise in maternal cue 
double sdmat = 0.0;

// whether survival selection has a sigmoidal or a quadratic function
bool sigmoidal_survival = true;

bool laplace = true;

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


// migration rates
double m = 0.0;

//write out the data every nth generation
int data_nth_generation = 10;

// survival statistics which I obtain in the survive()
// function and then use in the write_data() function
double mean_survival[2] = {0.0,0.0};
double var_survival[2] = {0.0,0.0};


struct Patch
{
    // number of hermaphroditic breeder
    Individual breeders[NBreeder];

    // next generation breeders
    Individual breeders_t1[NBreeder];
    
    int n_breeders; // number of breeders currently in patch
    bool envt_high;  // high-state patch yes/no
};

Patch Pop[NPatches];

// survival function
double survival_probability(double const ad_phen, bool const state_high)
{
    // eqns 3,4 Leimar & McNamara (2015) Amnat
    if (sigmoidal_survival)
    {
        double phen = state_high ? -ad_phen : ad_phen;
        return(1.0 / (1.0 + exp(survival_scalar[state_high] + 6.0 * phen)));
    }

    // eqns 1,2 Leimar & McNamara (2015) Amnat
    return(state_high ? 
            1.0 - survival_scalar[0] * (1.0 - ad_phen) * (1.0 - ad_phen)
            :
            1.0 - survival_scalar[0] * ad_phen * ad_phen
            );
}


// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{
    sigmoidal_survival = atoi(argv[1]);
    laplace = atoi(argv[2]);
    p = atof(argv[3]);
    survival_scalar[0] = atof(argv[4]);
    survival_scalar[1] = atof(argv[5]);
    qmat = atof(argv[6]);
    qjuv = atof(argv[7]);
    nloci_g = atoi(argv[8]);
    init_g = atof(argv[9]);
    init_amat = atof(argv[10]);
    init_ajuv = atof(argv[11]);
    init_agen = atof(argv[12]);
    init_bmat_phen = atof(argv[13]);
    init_bmat_envt = atof(argv[14]);
    gmin = atof(argv[15]);
    gmax = atof(argv[16]);
    amin = atof(argv[17]);
    amax = atof(argv[18]);
    bmin = atof(argv[19]);
    bmax = atof(argv[20]);
    sdmat = atof(argv[21]);

    mu_g = atof(argv[22]);
    mu_amat = atof(argv[23]);
    mu_ajuv = atof(argv[24]);
    mu_agen = atof(argv[25]);
    mu_bmat_phen = atof(argv[26]);
    mu_bmat_envt = atof(argv[27]);
    sdmu_a = atof(argv[28]);
    sdmu_b = atof(argv[29]);
    sdmu_g = atof(argv[30]);
    m = atof(argv[31]);

}

// write down all parameters in the file
void write_parameters(ofstream &DataFile)
{
    DataFile << endl << endl
        << "sigmoidal_survival;" << sigmoidal_survival << ";"<< endl
        << "laplace;" << laplace << ";"<< endl
        << "p;" << p << ";"<< endl
        << "qmat;" << qmat << ";"<< endl
        << "qjuv;" << qjuv << ";"<< endl
        << "nloci_g;" << nloci_g << ";"<< endl
        << "init_g;" << init_g << ";"<< endl
        << "init_amat;" << init_amat << ";"<< endl
        << "init_ajuv;" << init_ajuv << ";"<< endl
        << "init_agen;" << init_agen << ";"<< endl
        << "init_bmat_phen;" << init_bmat_phen << ";"<< endl
        << "init_bmat_envt;" << init_bmat_envt << ";"<< endl
        << "gmin;" << gmin << ";"<< endl
        << "gmax;" << gmax << ";"<< endl
        << "amin;" << amin << ";"<< endl
        << "amax;" << amax << ";"<< endl
        << "bmin;" << bmin << ";"<< endl
        << "bmax;" << bmax << ";"<< endl
        << "sdmat;" << sdmat << ";"<< endl
        << "mu_g;" << mu_g << ";"<< endl
        << "sdmu_g;" << sdmu_g << ";"<< endl
        << "mu_amat;" << mu_amat << ";"<< endl
        << "mu_ajuv;" << mu_ajuv << ";"<< endl
        << "mu_agen;" << mu_agen << ";"<< endl
        << "mu_bmat_phen;" << mu_bmat_phen << ";"<< endl
        << "mu_bmat_envt;" << mu_bmat_envt << ";"<< endl
        << "sdmu_a;" << sdmu_a << ";"<< endl
        << "sdmu_b;" << sdmu_b << ";"<< endl
        << "m;" << m << ";"<< endl
        << "survival_scalar0;" << survival_scalar[0] << ";"<< endl
        << "survival_scalar1;" << survival_scalar[1] << ";"<< endl
        << "seed;" << seed << ";" << endl;

} // void write_parameters(ofstream &DataFile)

void write_dist(ofstream &DataFile)
{
    double g;
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            DataFile << patch_i << ";" 
                << breeder_i << ";"
                << Pop[patch_i].breeders[breeder_i].ad_phen << ";"
                << Pop[patch_i].breeders[breeder_i].ad_mat << ";"
                << Pop[patch_i].breeders[breeder_i].xmat << ";"
                << 0.5 * (Pop[patch_i].breeders[breeder_i].agen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1]) << ";"
                << 0.5 * (Pop[patch_i].breeders[breeder_i].ajuv[0]
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1]) << ";"
                << 0.5 * (Pop[patch_i].breeders[breeder_i].amat[0]
                    +
                    Pop[patch_i].breeders[breeder_i].amat[1]) << ";"
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_phen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1]) << ";"
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_envt[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1]) << ";";

            g = 0.0;

            for (int g_i = 0; g_i < nloci_g; ++g_i)
            {
                g += 0.5 * (
                    Pop[patch_i].breeders[breeder_i].g[0][g_i]
                    +
                    Pop[patch_i].breeders[breeder_i].g[1][g_i]);
            }

            DataFile << g << ";"
                << Pop[patch_i].envt_high << ";"
                << Pop[patch_i].breeders[breeder_i].cue_ad_envt_high << ";"
                << Pop[patch_i].breeders[breeder_i].cue_juv_envt_high << ";"
                << endl;
        } // end for (int breeder_i = 0; breeder_i < NPatches; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)
} // end  void write_dist(ofstream &DataFile)


// list of the data headers at the start of the file
// in which the distribution of evolved trait values
// and states is written at the end of the file
void write_data_headers_dist(ofstream &DataFile)
{
    DataFile 
        << "patch_id;" 
        << "id;" 
        << "ad_phen;" 
        << "ad_mat;" 
        << "xmat;" 
        << "agen;" 
        << "ajuv;" 
        << "amat;" 
        << "bmat_phen;" 
        << "bmat_envt;" 
        << "g;" 
        << "envt;" 
        << "cue_ad_envt_high;" 
        << "cue_juv_envt_high;" 
        << endl; 
}

// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile 
        << "generation;" 
        << "mean_ad_phen;" 
        << "mean_agen;" 
        << "mean_ajuv;" 
        << "mean_amat;" 
        << "mean_bmat_phen;" 
        << "mean_bmat_envt;" 
        << "mean_g;" 
        << "var_ad_phen;" 
        << "var_agen;" 
        << "var_ajuv;" 
        << "var_amat;" 
        << "var_bmat_phen;" 
        << "var_bmat_envt;" 
        << "var_g;" 
        << "freq_high;" 
        << "mean_surv0;" 
        << "mean_surv1;" 
        << "var_surv0;" 
        << "var_surv1;" 
        << endl; 
}

// write data both for winter and summer populations
void write_stats(ofstream &DataFile, int generation, int timestep)
{
    // variables to store means and variances
    double mean_ad_phen = 0.0;
    double ss_ad_phen = 0.0;
   
    double mean_agen = 0.0;
    double ss_agen = 0.0;
   
    double mean_amat = 0.0;
    double ss_amat = 0.0;
    
    double mean_ajuv = 0.0;
    double ss_ajuv = 0.0;
    
    double mean_bmat_phen = 0.0;
    double ss_bmat_phen = 0.0;
    
    double mean_bmat_envt = 0.0;
    double ss_bmat_envt = 0.0;

    double mean_g = 0.0;
    double ss_g = 0.0;

    double freq_high = 0.0;

    // auxiliary variables to calculate an individual's phenotype
    double g, ad_phen, agen, amat, ajuv, bmat_phen, bmat_envt;

    // summing means and sums of squares over all patches and breeders
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        freq_high += Pop[patch_i].envt_high;

        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            for (int g_i = 0; g_i < nloci_g; ++g_i)
            {
                g =  0.5 * (
                        Pop[patch_i].breeders[breeder_i].g[0][g_i]
                        +
                        Pop[patch_i].breeders[breeder_i].g[1][g_i]
                        );

                mean_g += g;
                ss_g += g * g;
            } // end for (int g_i = 0; g_i < nloci_g; ++g_i)

            // adult phenotype
            ad_phen = Pop[patch_i].breeders[breeder_i].ad_phen;
            mean_ad_phen += ad_phen;
            ss_ad_phen += ad_phen * ad_phen;

            // sensitivity to genetic cues
            agen = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].agen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1] 
                    );

            mean_agen += agen;
            ss_agen += agen * agen;
            
            // sensitivity to maternal cues
            amat = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].amat[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].amat[1] 
                    );

            mean_amat += amat;
            ss_amat += amat * amat;
            
            // sensitivity to juvenile cues
            ajuv = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].ajuv[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1] 
                    );

            mean_ajuv += ajuv;
            ss_ajuv += ajuv * ajuv;
            
            // maternal sensitivity to phenotypic cues
            bmat_phen = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_phen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1] 
                    );

            mean_bmat_phen += bmat_phen;
            ss_bmat_phen += bmat_phen * bmat_phen;
            
            // maternal sensitivity to environmental cues
            bmat_envt = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_envt[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1] 
                    );

            mean_bmat_envt += bmat_envt;
            ss_bmat_envt += bmat_envt * bmat_envt;
        } // end for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)

    mean_ad_phen /= NPatches * NBreeder;
    double var_ad_phen = ss_ad_phen / NPatches * NBreeder - mean_ad_phen * mean_ad_phen;
   
    mean_agen /= NPatches * NBreeder;
    double var_agen = ss_agen / NPatches * NBreeder - mean_agen * mean_agen;
   
    mean_amat /= NPatches * NBreeder;
    double var_amat = NPatches * NBreeder - mean_amat * mean_amat;
    
    mean_ajuv /= NPatches * NBreeder;
    double var_ajuv = ss_ajuv / NPatches * NBreeder - mean_ajuv * mean_ajuv;
    
    mean_bmat_phen /= NPatches * NBreeder;
    double var_bmat_phen = ss_bmat_phen / NPatches * NBreeder - mean_bmat_phen * mean_bmat_phen;
    
    mean_bmat_envt /= NPatches * NBreeder;
    double var_bmat_envt = ss_bmat_envt / NPatches * NBreeder - mean_bmat_envt * mean_bmat_envt;

    mean_g /= NPatches * NBreeder;
    double var_g = ss_g / NPatches * NBreeder - mean_g * mean_g;

    freq_high /= NPatches;

    DataFile << generation << ";"
        << mean_ad_phen << ";"
        << mean_agen << ";"
        << mean_ajuv << ";"
        << mean_amat << ";"
        << mean_bmat_phen << ";"
        << mean_bmat_envt << ";"
        << mean_g << ";"
        << var_ad_phen << ";"
        << var_agen << ";"
        << var_ajuv << ";"
        << var_amat << ";"
        << var_bmat_phen << ";"
        << var_bmat_envt << ";"
        << var_g << ";" 
        << freq_high << ";" 
        << mean_survival[0] << ";" 
        << mean_survival[1] << ";" 
        << var_survival[0] << ";" 
        << var_survival[1] << ";" 
        << endl;
}

// initialize the population at the start of the simulation
void init_population()
{
    // auxiliary variable whether mothers perceived a cue
    // that the environment is in a high state
    bool cue_ad_envt_high;

    // loop through all individuals 
    // and assign them values for the cue loci
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // patch in a high state or not
        Pop[patch_i].envt_high = uniform(rng_r) < p;
       
        // cue given to adults
        cue_ad_envt_high = uniform(rng_r) < qmat ?
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
                    Pop[patch_i].breeders[breeder_i].g[allele_i].push_back(init_g);
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
            
            // initialize environmental cue values
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_ad_envt_high;

            Pop[patch_i].breeders[breeder_i].cue_juv_envt_high = 
                cue_ad_envt_high;

            // as this is generation t=0, forget about maternal cues for now
            Pop[patch_i].breeders[breeder_i].ad_phen = 1.0 /
                (1.0 + exp(-init_agen * init_g - init_ajuv * cue_ad_envt_high));
        } // end for breeder_i
    } // end for patch_i
} // end void init_population()

// mutation of a certain allele with value val
// given mutation rate mu and mutational distribution stdev sdmu
double mutation(double val, double mu, double sdmu)
{
    if (uniform(rng_r) < mu)
    {
        if (laplace)
        {
            // using a uniform dist to draw numbers from a Laplace dist
            // see https://en.wikipedia.org/wiki/Laplace_distribution 
            uniform_real_distribution<> mutational_effect(-0.5 + 0.0000001, 0.5);
            double U = mutational_effect(rng_r);

            double sgnU = U < 0.0 ? -1 : U > 0.0 ? 1.0 : 0.0;

            // effect size of Laplace
            double x = -sdmu * sgnU * log(1.0 - 2 * fabs(U));

            val += x;
        }
        else
        {
            normal_distribution<> mutational_effect(0.0, sdmu);
            val += mutational_effect(rng_r);
        }
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

    double allelic_val;

    assert((int)mother.g[0].size() == nloci_g);
    
    // reset arrays for the genetic cue values
    if (offspring.g[0].size() > 0)
    {
        offspring.g[0].clear();
        offspring.g[1].clear();
    }

    // inherit genetic cue values
    for (int g_loc_i = 0; g_loc_i < nloci_g; ++g_loc_i)
    {
        // maternal values 
        allelic_val = mutation(
                mother.g[allele_sample(rng_r)][g_loc_i],
                mu_g,
                sdmu_g);

        clamp(allelic_val, gmin, gmax);

        offspring.g[0].push_back(allelic_val);
        
        // paternal values 
        allelic_val = mutation(
                father.g[allele_sample(rng_r)][g_loc_i],
                mu_g,
                sdmu_g);

        clamp(allelic_val, gmin, gmax);

        offspring.g[1].push_back(allelic_val);

        sum_genes += 0.5 * (offspring.g[0][g_loc_i] + offspring.g[1][g_loc_i]);
    }
    
    assert((int)offspring.g[0].size() == nloci_g);

    // inheritance of maternal cue values 
    offspring.amat[0] = mutation(
            mother.amat[allele_sample(rng_r)],
            mu_amat,
            sdmu_a);

    clamp(offspring.amat[0], amin, amax);

    offspring.amat[1] = mutation(
            father.amat[allele_sample(rng_r)],
            mu_amat,
            sdmu_a);

    clamp(offspring.amat[1], amin, amax);

    double amat_phen = 0.5 * (offspring.amat[0] + offspring.amat[1]);
   

    // inheritance of juvenile cue values 
    offspring.ajuv[0] = mutation(
            mother.ajuv[allele_sample(rng_r)],
            mu_ajuv,
            sdmu_a);

    clamp(offspring.ajuv[0], amin, amax);

    offspring.ajuv[1] = mutation(
            father.ajuv[allele_sample(rng_r)],
            mu_ajuv,
            sdmu_a);

    clamp(offspring.ajuv[1], amin, amax);

    double ajuv_phen = 0.5 * (offspring.ajuv[0] + offspring.ajuv[1]);
   


    // inheritance of genetic cue values 
    offspring.agen[0] = mutation(
            mother.agen[allele_sample(rng_r)],
            mu_agen,
            sdmu_a);

    clamp(offspring.agen[0], amin, amax);

    offspring.agen[1] = mutation(
            father.agen[allele_sample(rng_r)],
            mu_agen,
            sdmu_a);

    clamp(offspring.agen[1], amin, amax);

    double agen_phen = 0.5 * (offspring.agen[0] + offspring.agen[1]);


    // inheritance of maternal phenotypic cue values 
    offspring.bmat_phen[0] = mutation(
            mother.bmat_phen[allele_sample(rng_r)],
            mu_bmat_phen,
            sdmu_b);

    clamp(offspring.bmat_phen[0], bmin, bmax);

    offspring.bmat_phen[1] = mutation(
            father.bmat_phen[allele_sample(rng_r)],
            mu_bmat_phen,
            sdmu_b);

    clamp(offspring.bmat_phen[1], bmin, bmax);

    // inheritance of maternal environmental cue values 
    offspring.bmat_envt[0] = mutation(
            mother.bmat_envt[allele_sample(rng_r)],
            mu_bmat_envt,
            sdmu_b);

    clamp(offspring.bmat_envt[0], bmin, bmax);

    offspring.bmat_envt[1] = mutation(
            father.bmat_envt[allele_sample(rng_r)],
            mu_bmat_envt,
            sdmu_b);

    clamp(offspring.bmat_envt[1], bmin, bmax);

    // kid receives juvenile cue
    offspring.cue_juv_envt_high = uniform(rng_r) < qjuv ? 
        offspring_envt : !offspring_envt;

    // adult cue will be received after potential envt'al change
    //
    // has the mother observed a high cue or a low one?
    double dmat_weighting = mother.cue_ad_envt_high ? -1.0 : 1.0;

    offspring.ad_mat = mother.ad_phen;

    // generate maternal cue
    double xoff = 1.0 /
        (1.0 + exp(
                   -0.5 * (mother.bmat_phen[0] + mother.bmat_phen[1]) * (mother.ad_phen - 0.5) + 
                   dmat_weighting * 0.5 * (mother.bmat_envt[0] + mother.bmat_envt[1])));

    // noise in the maternal cue
    normal_distribution<> maternal_noise(0.0, sdmat);

    offspring.xmat = xoff + maternal_noise(rng_r);

    clamp(offspring.xmat, 0.0, 1.0);

    offspring.ad_phen = 1.0 / 
        (1.0 + exp(
                   -amat_phen * offspring.xmat +
                   -agen_phen * sum_genes +
                   -ajuv_phen * offspring.cue_juv_envt_high));

} // end create_offspring()

void survive()
{
    // set up a random number generator to 
    // sample from the remaining breeders
    uniform_int_distribution<> random_patch(
            0,
            NPatches - 1);


    // reset survival statistics
    mean_survival[0] = 0.0;
    var_survival[0] = 0.0;
    mean_survival[1] = 0.0;
    var_survival[1] = 0.0;

    // keep track of the number of breeders
    // in high patches to calculate survival stats
    int n_high_patches = 0.0;

    double surv = 0.0;

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        assert(Pop[patch_i].n_breeders > 0);
        assert(Pop[patch_i].n_breeders <= NBreeder);
        
        // keep track of the number of breeders in a
        // high envt
        if (Pop[patch_i].envt_high)
        {
            n_high_patches += NBreeder;
        }

        // breeders endure survival selection
        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            // calculate survival probability
            surv = survival_probability(
                        Pop[patch_i].breeders[breeder_i].ad_phen
                        ,Pop[patch_i].envt_high);

            // store the survival value in this envt
            // for sake of statistics
            mean_survival[Pop[patch_i].envt_high] += surv;

            // store sum of squares now, calc variance later
            var_survival[Pop[patch_i].envt_high] += surv * surv;

            // individual dies if random uniform number is larger 
            // than the survival probability
            if (uniform(rng_r) > surv)
            {
                // delete individual
                Pop[patch_i].breeders[breeder_i] = 
                    Pop[patch_i].breeders[Pop[patch_i].n_breeders - 1];

                --breeder_i;
                --Pop[patch_i].n_breeders;
            }
        } // end for (int breeder_i


        // environmental change
        if (p < uniform(rng_r))
        {
            Pop[patch_i].envt_high = !Pop[patch_i].envt_high;
        }
    } // end for int patch_i

    // finalize survival statistics
    mean_survival[0] /= NPatches * NBreeder - n_high_patches;
    mean_survival[1] /= n_high_patches;

    // var = E[xx] - E[x]E[x]
    var_survival[0] = var_survival[0] / (NPatches * NBreeder - n_high_patches)
        - mean_survival[0] * mean_survival[0];
    
    var_survival[1] = var_survival[1] / n_high_patches
        - mean_survival[1] * mean_survival[1];

    // finalize statistics
} // end survive_reproduce()

void replace()
{
    // randomly chosen remote patch to obtain
    // individuals from
    int random_remote_patch;
        
    // set up a random number generator to 
    // sample remote patches
    uniform_int_distribution<> patch_sampler(
            0,
            NPatches - 1);

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            if (uniform(rng_r) > m 
                    &&
                    Pop[patch_i].n_breeders > 0)
            {
                // set up a random number generator to 
                // sample from the remaining breeders
                uniform_int_distribution<> random_local_breeder(
                        0,
                        Pop[patch_i].n_breeders - 1);

                create_offspring(
                        Pop[patch_i].breeders[random_local_breeder(rng_r)]
                        ,Pop[patch_i].breeders[random_local_breeder(rng_r)]
                        ,Pop[patch_i].breeders_t1[breeder_i]
                        ,Pop[patch_i].envt_high
                );
                
                assert((int)Pop[patch_i].breeders_t1[breeder_i].g[0].size() == nloci_g);
            }
            else
            {
                do {

                    random_remote_patch = patch_sampler(rng_r);

                }
                while(Pop[random_remote_patch].n_breeders < 1);
        
                uniform_int_distribution<> random_remote_breeder(
                0,
                Pop[random_remote_patch].n_breeders - 1);

                create_offspring(
                        Pop[random_remote_patch].breeders[random_remote_breeder(rng_r)]
                        ,Pop[random_remote_patch].breeders[random_remote_breeder(rng_r)]
                        ,Pop[patch_i].breeders_t1[breeder_i]
                        ,Pop[random_remote_patch].envt_high
                );
            
                assert((int)Pop[patch_i].breeders_t1[breeder_i].g[0].size() == nloci_g);
            }
        } // for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
    } // end for (int patch_i = 0

    bool cue_ad_envt_high;

    // all new breeders established, copy them over
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // calculate adult cue value supplied to mothers
        // the cue value is the same for all mothers
        cue_ad_envt_high = uniform(rng_r) < qmat ? 
            Pop[patch_i].envt_high 
            : 
            !Pop[patch_i].envt_high;
        
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            Pop[patch_i].breeders[breeder_i] = 
                Pop[patch_i].breeders_t1[breeder_i];
    
            assert((int)Pop[patch_i].breeders_t1[breeder_i].g[0].size() == nloci_g);
            
            Pop[patch_i].n_breeders = NBreeder;
        
            // give breeder an environmental cue as adult
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_ad_envt_high;
        }
    } // end for (int patch_i = 0
}

// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    string filename = "sim_cue_integration";
    create_filename(filename);
    ofstream DataFile(filename.c_str());  // output file 

    // output file to write out the complete 
    // trait and state distribution of individuals
    // in the last generation
    string filename_final_dist = filename + "_dist";
    ofstream DataFileDist(filename_final_dist.c_str());  

    // get command line arguments
    init_arguments(argc, argv);

    // write headers to the datafile
    write_data_headers(DataFile);
    
    write_data_headers_dist(DataFileDist);

    // initialize the population
    init_population();

    for (int generation = 0; generation < number_generations; ++generation)
    {
        survive();

        replace();

        if (generation % data_nth_generation == 0)
        {
            write_stats(DataFile, generation, 2);
        }
    }

    write_dist(DataFileDist);
    
    write_parameters(DataFile);
}
