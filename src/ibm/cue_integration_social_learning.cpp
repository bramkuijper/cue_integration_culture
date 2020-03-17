// Extending Leimar & McNamara's cue integration model
// to include cultural transmission
// Bram Kuijper 
// 2019
//
#define DEBUG

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>


// the individual class, which defines properties of each
// individual used here
#include "individual.hpp"

// C++ random number generation unsigned int seed = get_nanoseconds();
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r{seed};
std::uniform_real_distribution<> uniform(0.0,1.0);

// file base name
std::string base_name;
// parameters & variables:
// a lot of the parameters declared here will be overridden
// in the init_arguments function which reads parameters from the command line

// number of individuals in population
const int NPatches = 400;
const int NBreeder = 100;

// number of generations
int number_generations = 75000;

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

double max_error_conform_horiz = 0;
double max_error_prestige_horiz = 0;
double max_error_conform_vert = 0;
double max_error_prestige_vert = 0;

double max_error_mat_phen = 0;
double max_error_mat_envt = 0;

// whether survival selection has a sigmoidal or a quadratic function
bool sigmoidal_survival = true;

bool laplace = true;

// number of loci underlying the genetic cue
int nloci_g = 3;

// initial values
// these will be overridden when the function
// init_arguments is called
double init_g = 0.0;
double init_aintercept = 0.0;
double init_ajuv = 0.0;
double init_agen = 0.0;
double init_bmat_phen = 0.0;
double init_bmat_envt = 0.0;
double init_hc = 0.0;
double init_hp = 0.0;
double init_vc = 0.0;
double init_vp = 0.0;

// ranges for traits
double gmin = 0.0;
double gmax = 0.0;
double amin = 0.0;
double amax = 0.0;
double bmin = 0.0;
double bmax = 0.0;

// mutation rates
double mu_g = 0.0;
double mu_aintercept = 0.0;
double mu_ajuv = 0.0;
double mu_agen = 0.0;
double mu_bmat_phen = 0.0;
double mu_bmat_envt = 0.0;
double mu_hp = 0.0;
double mu_hc = 0.0;
double mu_vp = 0.0;
double mu_vc = 0.0;

double sdmu_a = 0.0;
double sdmu_b = 0.0;
double sdmu_g = 0.0;

bool juvenile_survival = false;


// number of individuals to sample
// in case of horizontal social learning
int nph = 0; // performance bias
int nch = 0; // conformity bias

// number of individuals to sample
// in case of vertical social learning
int npv = 0; // performance bias
int ncv = 0; // conformity bias

// migration rates
double m = 0.0;

//write out the data 
//every nth generation
int data_nth_generation = 50;

// survival statistics which I obtain in the survive()
// function and then use in the write_data() function
double mean_survival[2] = {0.0,0.0};
double var_survival[2] = {0.0,0.0};


// the patch with its breeders
struct Patch
{
    // number of hermaphroditic breeder
    Individual breeders[NBreeder];

    // next generation breeders
    Individual breeders_t1[NBreeder];
    
    int n_breeders; // number of breeders currently in patch
    int n_recruits; // number of breeders currently in patch
    bool envt_high;  // high-state patch yes/no
};

Patch Pop[NPatches];

// survival function
double survival_probability(double const phen_ad, bool const state_high)
{
    // eqns 3,4 Leimar & McNamara (2015) Amnat
    if (sigmoidal_survival)
    {
        double phen = state_high ? -phen_ad : phen_ad;
        return(1.0 / (1.0 + exp(survival_scalar[state_high] + 6.0 * phen)));
    }

    // eqns 1,2 Leimar & McNamara (2015) Amnat
    return(state_high ? 
            1.0 - survival_scalar[0] * (1.0 - phen_ad) * (1.0 - phen_ad)
            :
            1.0 - survival_scalar[0] * phen_ad * phen_ad
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
    max_error_conform_horiz = atof(argv[8]);
    max_error_prestige_horiz = atof(argv[9]);
    max_error_conform_vert = atof(argv[10]);
    max_error_prestige_vert = atof(argv[11]);
    max_error_mat_phen = atof(argv[12]);
    max_error_mat_envt = atof(argv[13]);
    nloci_g = atoi(argv[14]);
    
    init_g = atof(argv[15]);
    init_aintercept = atof(argv[16]);
    init_ajuv = atof(argv[17]);

    init_agen = atof(argv[18]);
    init_bmat_phen = atof(argv[19]);
    init_bmat_envt = atof(argv[20]);

    init_hp = atof(argv[21]);
    init_hc = atof(argv[22]);
    init_vp = atof(argv[23]);

    init_vc = atof(argv[24]);

    gmin = atof(argv[25]);
    gmax = atof(argv[26]);
    amin = atof(argv[27]);
    amax = atof(argv[28]);
    bmin = atof(argv[29]);
    bmax = atof(argv[30]);

    mu_g = atof(argv[31]);
    mu_aintercept = atof(argv[32]);
    mu_ajuv = atof(argv[33]);
    mu_agen = atof(argv[34]);
    mu_bmat_phen = atof(argv[35]);
    mu_bmat_envt = atof(argv[36]);
    mu_hp = atof(argv[37]);
    mu_hc = atof(argv[38]);
    mu_vp = atof(argv[39]);
    mu_vc = atof(argv[40]);
    sdmu_a = atof(argv[41]);
    sdmu_b = atof(argv[42]);
    sdmu_g = atof(argv[43]);
    m = atof(argv[44]);
    nph = atoi(argv[45]);
    nch = atoi(argv[46]);
    npv = atoi(argv[47]);
    ncv = atoi(argv[48]);
    juvenile_survival = atoi(argv[49]);
    base_name = argv[50];
}

// write down all parameters to the file DataFile
void write_parameters(std::ofstream &DataFile)
{
    DataFile << std::endl << std::endl
        << "sigmoidal_survival;" << sigmoidal_survival << ";"<< std::endl
        << "laplace;" << laplace << ";"<< std::endl
        << "p;" << p << ";"<< std::endl
        << "qmat;" << qmat << ";"<< std::endl
        << "qjuv;" << qjuv << ";"<< std::endl
        << "nloci_g;" << nloci_g << ";"<< std::endl
        << "init_g;" << init_g << ";"<< std::endl
        << "init_aintercept;" << init_aintercept << ";"<< std::endl
        << "init_ajuv;" << init_ajuv << ";"<< std::endl
        << "init_agen;" << init_agen << ";"<< std::endl
        << "init_bmat_phen;" << init_bmat_phen << ";"<< std::endl
        << "init_bmat_envt;" << init_bmat_envt << ";"<< std::endl
        << "init_hp;" << init_hp << ";"<< std::endl
        << "init_hc;" << init_hc << ";"<< std::endl
        << "init_vp;" << init_vp << ";"<< std::endl
        << "init_vc;" << init_vc << ";"<< std::endl
        << "gmin;" << gmin << ";"<< std::endl
        << "gmax;" << gmax << ";"<< std::endl
        << "amin;" << amin << ";"<< std::endl
        << "amax;" << amax << ";"<< std::endl
        << "bmin;" << bmin << ";"<< std::endl
        << "bmax;" << bmax << ";"<< std::endl
        << "mu_g;" << mu_g << ";"<< std::endl
        << "sdmu_g;" << sdmu_g << ";"<< std::endl
        << "mu_ajuv;" << mu_ajuv << ";"<< std::endl
        << "mu_agen;" << mu_agen << ";"<< std::endl
        << "mu_aintercept;" << mu_aintercept << ";"<< std::endl
        << "juvenile_survival;" << juvenile_survival << ";"<< std::endl
        << "mu_bmat_phen;" << mu_bmat_phen << ";"<< std::endl
        << "mu_bmat_envt;" << mu_bmat_envt << ";"<< std::endl
        << "mu_vc;" << mu_vc << ";"<< std::endl
        << "mu_vp;" << mu_vp << ";"<< std::endl
        << "mu_hc;" << mu_hc << ";"<< std::endl
        << "mu_hp;" << mu_hp << ";"<< std::endl
        << "sdmu_a;" << sdmu_a << ";"<< std::endl
        << "sdmu_b;" << sdmu_b << ";"<< std::endl
        << "m;" << m << ";"<< std::endl
        << "nph;" << nph << ";"<< std::endl
        << "nch;" << nch << ";"<< std::endl
        << "npv;" << npv << ";"<< std::endl
        << "ncv;" << ncv << ";"<< std::endl
        << "survival_scalar0;" << survival_scalar[0] << ";"<< std::endl
        << "survival_scalar1;" << survival_scalar[1] << ";"<< std::endl
        << "seed;" << seed << ";" << std::endl;

} // void write_parameters(std::ofstream &DataFile)


// write all properties of all individuals
// to the file DataFile (to obtain information about
// the distribution of traits)
void write_dist(std::ofstream &DataFile)
{
    double g; // auxiliary variable to temporarily 
                // store trait expression

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            DataFile << patch_i << ";" 
                << breeder_i << ";"
                << Pop[patch_i].breeders[breeder_i].phen_ad << ";"
                << Pop[patch_i].breeders[breeder_i].phen_ad_logistic << ";"
                << Pop[patch_i].breeders[breeder_i].phen_juv << ";"
                << Pop[patch_i].breeders[breeder_i].phen_juv_logistic << ";"
                << Pop[patch_i].breeders[breeder_i].phen_mat << ";"
                << Pop[patch_i].breeders[breeder_i].phen_mat_error << ";"
                << Pop[patch_i].breeders[breeder_i].maternal_envt_cue_error << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_vert << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_vert_error << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_horiz << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_horiz_error << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_vert << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_vert_error << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_horiz << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_horiz_error << ";"
                
                // the intercept locus
                << 0.5 * (Pop[patch_i].breeders[breeder_i].aintercept[0]
                    +
                    Pop[patch_i].breeders[breeder_i].aintercept[1]) << ";"

                // agen
                << 0.5 * (Pop[patch_i].breeders[breeder_i].agen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1]) << ";"

                // ajuv
                << 0.5 * (Pop[patch_i].breeders[breeder_i].ajuv[0]
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1]) << ";"

                // bmat_phen 
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_phen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1]) << ";"
                
                // bmat_envt
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_envt[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1]) << ";"

                // hc
                << 0.5 * (Pop[patch_i].breeders[breeder_i].hc[0]
                    +
                    Pop[patch_i].breeders[breeder_i].hc[1]) << ";"

                // hp
                << 0.5 * (Pop[patch_i].breeders[breeder_i].hp[0]
                    +
                    Pop[patch_i].breeders[breeder_i].hp[1]) << ";"
                
                // vc
                << 0.5 * (Pop[patch_i].breeders[breeder_i].vc[0]
                    +
                    Pop[patch_i].breeders[breeder_i].vc[1]) << ";"

                // vp
                << 0.5 * (Pop[patch_i].breeders[breeder_i].vp[0]
                    +
                    Pop[patch_i].breeders[breeder_i].vp[1]) << ";";

            g = 0.0;

            // genetic cue
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
                << Pop[patch_i].breeders[breeder_i].cue_juv_envt_high << ";";
            
            // now output the whole phenotypic variance stuff
            DataFile 
                // agen X g
                << 0.5 * (Pop[patch_i].breeders[breeder_i].agen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1]) * g << ";"
                
                // bmat_phen X phen_mat_error
                << 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_phen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1]) * Pop[patch_i].breeders[breeder_i].phen_mat_error << ";"
                
                // bmat_envt X maternal_envt_cue_error
                << 0.5 * (Pop[patch_i].breeders[breeder_i].bmat_envt[0]
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1]) * Pop[patch_i].breeders[breeder_i].maternal_envt_cue_error << ";"

                // ajuv X cue_juv_envt_high
                << 0.5 * (Pop[patch_i].breeders[breeder_i].ajuv[0]
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1]) * Pop[patch_i].breeders[breeder_i].cue_juv_envt_high << ";"
               
                // vp_X_phen_prestige_vert_error
                << 0.5 * (Pop[patch_i].breeders[breeder_i].vp[0]
                    +
                    Pop[patch_i].breeders[breeder_i].vp[1]) * Pop[patch_i].breeders[breeder_i].phen_prestige_vert_error << ";"

                // vc X xconformist_vert_error
                << 0.5 * (Pop[patch_i].breeders[breeder_i].vc[0]
                    +
                    Pop[patch_i].breeders[breeder_i].vc[1]) * Pop[patch_i].breeders[breeder_i].xconformist_vert_error << ";"
                
                // hp X phen_prestige_horiz_error
                << 0.5 * (Pop[patch_i].breeders[breeder_i].hp[0]
                    +
                    Pop[patch_i].breeders[breeder_i].hp[1]) * Pop[patch_i].breeders[breeder_i].phen_prestige_horiz_error << ";"
                
                // hc X xconformist_horiz_error
                << 0.5 * (Pop[patch_i].breeders[breeder_i].hc[0]
                    +
                    Pop[patch_i].breeders[breeder_i].hc[1]) * Pop[patch_i].breeders[breeder_i].xconformist_horiz_error << ";"
                << std::endl;

        } // end for (int breeder_i = 0; breeder_i < NPatches; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)
} // end  void write_dist(std::ofstream &DataFile)


// list of the data headers at the start of the file
// in which the distribution of evolved trait values
// and states is written at the end of the file
void write_data_headers_dist(std::ofstream &DataFile)
{
    DataFile 
        << "patch_id;" 
        << "id;" 
        << "phen_ad;" 
        << "phen_ad_logistic;" 
        << "phen_juv;" 
        << "phen_juv_logistic;" 
        << "phen_mat;" 
        << "phen_mat_error;" 
        << "maternal_envt_cue_error;" 
        << "phen_prestige_vert;" 
        << "phen_prestige_vert_error;" 
        << "phen_prestige_horiz;" 
        << "phen_prestige_horiz_error;" 
        << "xconformist_vert;" 
        << "xconformist_vert_error;" 
        << "xconformist_horiz;" 
        << "xconformist_horiz_error;" 
        << "aintercept;" 
        << "agen;" 
        << "ajuv;" 
        << "bmat_phen;" 
        << "bmat_envt;" 
        << "hc;" 
        << "hp;" 
        << "vc;" 
        << "vp;" 
        << "g;" 
        << "envt;" 
        << "cue_ad_envt_high;" 
        << "cue_juv_envt_high;" 
        << "agen_X_g;" 
        << "bmat_phen_X_phen_mat_error;" 
        << "bmat_envt_X_maternal_envt_cue_error;" 
        << "ajuv_X_cue_juv_envt_high;" 
        << "vp_X_phen_prestige_vert_error;" 
        << "vc_X_xconformist_vert_error;"
        << "hp_X_phen_prestige_horiz_error;"
        << "hc_X_xconformist_horiz_error;"
        << std::endl; 
}

// list of the data headers at the start of the file
void write_data_headers(std::ofstream &DataFile)
{
    DataFile 
        << "generation;" 
        << "mean_phen_ad;" 
        << "mean_phen_juv;" 
        << "mean_phen_prestige_vert;" 
        << "mean_phen_prestige_horiz;" 
        << "mean_aintercept;" 
        << "mean_agen;" 
        << "mean_ajuv;" 
        << "mean_bmat_phen;" 
        << "mean_bmat_envt;" 
        << "mean_hc;" 
        << "mean_hp;" 
        << "mean_vc;" 
        << "mean_vp;" 
        << "mean_g;" 
        << "var_phen_ad;" 
        << "var_phen_juv;" 
        << "var_phen_prestige_vert;" 
        << "var_phen_prestige_horiz;" 
        << "var_aintercept;" 
        << "var_agen;" 
        << "var_ajuv;" 
        << "var_bmat_phen;" 
        << "var_bmat_envt;" 
        << "var_hc;" 
        << "var_hp;" 
        << "var_vc;" 
        << "var_vp;" 
        << "var_g;" 
        << "freq_high;" 
        << "mean_surv0;" 
        << "mean_surv1;" 
        << "var_surv0;" 
        << "var_surv1;" << std::endl;
}

// write data both for winter and summer populations
void write_stats(
        std::ofstream &DataFile, 
        int generation, 
        int timestep)
{
    // variables to store means and variances
    double mean_phen_ad = 0.0;
    double ss_phen_ad = 0.0;
    
    double mean_phen_juv = 0.0;
    double ss_phen_juv = 0.0;
    
    double mean_phen_prestige_vert = 0.0;
    double ss_phen_prestige_vert = 0.0;
    
    double mean_phen_prestige_horiz = 0.0;
    double ss_phen_prestige_horiz = 0.0;
    
    double mean_aintercept = 0.0;
    double ss_aintercept = 0.0;
   
    double mean_agen = 0.0;
    double ss_agen = 0.0;
   
    double mean_ajuv = 0.0;
    double ss_ajuv = 0.0;
    
    double mean_bmat_phen = 0.0;
    double ss_bmat_phen = 0.0;
    
    double mean_bmat_envt = 0.0;
    double ss_bmat_envt = 0.0;
    
    double mean_hc = 0.0;
    double ss_hc = 0.0;

    double mean_hp = 0.0;
    double ss_hp = 0.0;
    
    double mean_vc = 0.0;
    double ss_vc = 0.0;

    double mean_vp = 0.0;
    double ss_vp = 0.0;

    double mean_g = 0.0;
    double ss_g = 0.0;

    double freq_high = 0.0;
    
    // auxiliary variables to calculate an individual's phenotype
    double val;

    int total_individuals = 0;

    // summing means and sums of squares over all patches and breeders
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        freq_high += Pop[patch_i].envt_high;

        total_individuals += Pop[patch_i].n_breeders;

        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            val = 0.0;
            for (int g_i = 0; g_i < nloci_g; ++g_i)
            {
                val +=  0.5 * (
                        Pop[patch_i].breeders[breeder_i].g[0][g_i]
                        +
                        Pop[patch_i].breeders[breeder_i].g[1][g_i]
                        );

            } // end for (int g_i = 0; g_i < nloci_g; ++g_i)

            mean_g += val;
            ss_g += val * val;

            // adult phenotype
            val = Pop[patch_i].breeders[breeder_i].phen_ad;
            mean_phen_ad += val;
            ss_phen_ad += val * val;
            
            // juvenile phenotype
            val = Pop[patch_i].breeders[breeder_i].phen_juv;
            mean_phen_juv += val;
            ss_phen_juv += val * val;
            
            // prestige phenotype
            val = Pop[patch_i].breeders[breeder_i].phen_prestige_vert_error;
            mean_phen_prestige_vert += val;
            ss_phen_prestige_vert += val * val;
            
            // prestige phenotype
            val = Pop[patch_i].breeders[breeder_i].phen_prestige_horiz_error;
            mean_phen_prestige_horiz += val;
            ss_phen_prestige_horiz += val * val;

            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].aintercept[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].aintercept[1] 
                    );

            mean_aintercept += val;
            ss_aintercept += val * val;

            // sensitivity to genetic cues
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].agen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1] 
                    );

            mean_agen += val;
            ss_agen += val * val;

            // sensitivity to juvenile cues
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].ajuv[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1] 
                    );
            
            mean_ajuv += val;
            ss_ajuv += val * val;

            // maternal sensitivity to phenotypic cues
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_phen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1] 
                    );

            mean_bmat_phen += val;
            ss_bmat_phen += val * val;
            
            // maternal sensitivity to environmental cues
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_envt[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1] 
                    );

            mean_bmat_envt += val;
            ss_bmat_envt += val * val;
            
            // sensitivity to performance-based cues when learning
            // horizontally
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hp[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].hp[1] 
                    );

            mean_hp += val;
            ss_hp += val * val;
            
            // sensitivity to confirmity-based cues when learning
            // horizontally
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hc[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].hc[1] 
                    );

            mean_hc += val;
            ss_hc += val * val;
            
            
            // sensitivity to performance-based cues when learning
            // vertically 
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].vp[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].vp[1] 
                    );

            mean_vp += val;
            ss_vp += val * val;
            
            // sensitivity to confirmity-based cues when learning
            // vertically
            val = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].vc[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].vc[1] 
                    );

            mean_vc += val;
            ss_vc += val * val;
        } // end for (int breeder_i = 0; breeder_i < n_breeders; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)

    mean_phen_ad /= total_individuals;
    double var_phen_ad = ss_phen_ad / total_individuals - mean_phen_ad * mean_phen_ad;
    
    mean_phen_juv /= total_individuals;
    double var_phen_juv = ss_phen_juv / total_individuals - mean_phen_juv * mean_phen_juv;
    
    mean_phen_prestige_vert /= total_individuals;

    double var_phen_prestige_vert = ss_phen_prestige_vert / total_individuals 
        - mean_phen_prestige_vert * mean_phen_prestige_vert;
   
    mean_phen_prestige_horiz /= total_individuals;

    double var_phen_prestige_horiz = ss_phen_prestige_horiz / total_individuals 
        - mean_phen_prestige_horiz * mean_phen_prestige_horiz;


    mean_aintercept /= total_individuals;
    double var_aintercept = ss_aintercept / total_individuals - mean_aintercept * mean_aintercept;

    mean_agen /= total_individuals;
    double var_agen = ss_agen / total_individuals - mean_agen * mean_agen;
   
    mean_ajuv /= total_individuals;
    double var_ajuv = ss_ajuv / total_individuals -  mean_ajuv * mean_ajuv;
    
    mean_bmat_phen /= total_individuals;
    double var_bmat_phen = ss_bmat_phen / total_individuals - mean_bmat_phen * mean_bmat_phen;
    
    mean_bmat_envt /= total_individuals;
    double var_bmat_envt = ss_bmat_envt / total_individuals - mean_bmat_envt * mean_bmat_envt;
    
    mean_hp /= total_individuals;
    double var_hp = ss_hp / total_individuals - mean_hp * mean_hp;
    
    mean_hc /= total_individuals;
    double var_hc = ss_hc / total_individuals - mean_hc * mean_hc;
    
    mean_vp /= total_individuals;
    double var_vp = ss_vp / total_individuals - mean_vp * mean_vp;
    
    mean_vc /= total_individuals;
    double var_vc = ss_vc / total_individuals - mean_vc * mean_vc;

    mean_g /= total_individuals;
    double var_g = ss_g / total_individuals - mean_g * mean_g;
    
    freq_high /= NPatches;

    DataFile << generation << ";"
        << mean_phen_ad << ";"
        << mean_phen_juv << ";"
        << mean_phen_prestige_vert << ";"
        << mean_phen_prestige_horiz << ";"
        << mean_aintercept << ";"
        << mean_agen << ";"
        << mean_ajuv << ";"
        << mean_bmat_phen << ";"
        << mean_bmat_envt << ";"
        << mean_hc << ";"
        << mean_hp << ";"
        << mean_vc << ";"
        << mean_vp << ";"
        << mean_g << ";"
        << var_phen_ad << ";"
        << var_phen_juv << ";"
        << var_phen_prestige_vert << ";"
        << var_phen_prestige_horiz << ";"
        << var_aintercept << ";"
        << var_agen << ";"
        << var_ajuv << ";"
        << var_bmat_phen << ";"
        << var_bmat_envt << ";"
        << var_hc << ";"
        << var_hp << ";"
        << var_vc << ";"
        << var_vp << ";"
        << var_g << ";" 
        << freq_high << ";" 
        << mean_survival[0] << ";" 
        << mean_survival[1] << ";" 
        << var_survival[0] << ";" 
        << var_survival[1] << ";" << std::endl;
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
        Pop[patch_i].envt_high = uniform(rng_r) < 0.5;
       
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
                Pop[patch_i].breeders[breeder_i].aintercept[allele_i] = init_aintercept;

                // juvenile cue weighting
                Pop[patch_i].breeders[breeder_i].ajuv[allele_i] = init_ajuv;

                // genetic cue weighting
                Pop[patch_i].breeders[breeder_i].agen[allele_i] = init_agen;
                
                // maternal phenotypic cue weighting
                Pop[patch_i].breeders[breeder_i].bmat_phen[allele_i] = init_bmat_phen;
                
                // maternal phenotypic cue weighting
                Pop[patch_i].breeders[breeder_i].bmat_envt[allele_i] = init_bmat_envt;
                
                // weighting of a horizontally learnt prestige cue
                Pop[patch_i].breeders[breeder_i].hp[allele_i] = init_hp;
                
                // weighting of a horizontally learnt conformist cue
                Pop[patch_i].breeders[breeder_i].hc[allele_i] = init_hc;
                
                // weighting of a horizontally learnt prestige cue
                Pop[patch_i].breeders[breeder_i].vp[allele_i] = init_vp;
                
                // weighting of a horizontally learnt conformist cue
                Pop[patch_i].breeders[breeder_i].vc[allele_i] = init_vc;

            } //end for allele_i
            
            Pop[patch_i].n_breeders = NBreeder;
            
            // initialize environmental cue values
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_ad_envt_high;

            Pop[patch_i].breeders[breeder_i].cue_juv_envt_high = 
                cue_ad_envt_high;

            // give initial cue to individuals
            Pop[patch_i].breeders[breeder_i].phen_juv_logistic = 
                Pop[patch_i].breeders[breeder_i].phen_ad_logistic = 
                    init_aintercept
                    + init_agen * init_g
                    + init_ajuv * cue_ad_envt_high;

            // as this is generation t=0, forget about 
            // most cues when determining phenotype
            Pop[patch_i].breeders[breeder_i].phen_ad = 1.0 /
                (1.0 + exp(-Pop[patch_i].breeders[breeder_i].phen_ad_logistic));
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
            std::uniform_real_distribution<> mutational_effect(-0.5 + 0.0000001, 0.5);
            double U = mutational_effect(rng_r);

            double sgnU = U < 0.0 ? -1.0 : U > 0.0 ? 1.0 : 0.0;

            // effect size of Laplace
            double x = -sdmu/sqrt(2.0) * sgnU * log(1.0 - 2 * fabs(U));

            val += x;
        }
        else
        {
            std::normal_distribution<> mutational_effect(0.0, sdmu);
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
        ,bool const offspring_envt_high
        ,double const phen_prestige_vert
        ,double const xconformist_vert
        )
{
    // set up a bernoulli distribution that returns 0s or 1s
    // at equal probability to sample alleles from the first or
    // second set of chromosomes of diploid individuals
    std::bernoulli_distribution allele_sample(0.5);

    assert((int)mother.g[0].size() == nloci_g);
    
    // reset arrays for the genetic cue values
    if (offspring.g[0].size() > 0)
    {
        offspring.g[0].clear();
        offspring.g[1].clear();
    }

    // inherit genetic cue values
    // first auxiliary variables
    double sum_genes = 0.0;
    double allelic_val;

    // iterate over all gene loci and inherit
    for (int g_loc_i = 0; g_loc_i < nloci_g; ++g_loc_i)
    {
        // maternal values 
        allelic_val = mutation(
                mother.g[allele_sample(rng_r)][g_loc_i],
                mu_g,
                sdmu_g
                );

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
    
    // inheritance of intercept values 
    offspring.aintercept[0] = mutation(
            mother.aintercept[allele_sample(rng_r)],
            mu_aintercept,
            sdmu_a);

    clamp(offspring.aintercept[0], amin, amax);

    offspring.aintercept[1] = mutation(
            father.aintercept[allele_sample(rng_r)],
            mu_aintercept,
            sdmu_a);

    clamp(offspring.aintercept[1], amin, amax);

    double aintercept_phen = 0.5 * 
        (offspring.aintercept[0] + offspring.aintercept[1]);

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
   


    // inheritance of genetic cue sensitivity values 
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


    
    // inheritance of rel. sensitivity to 
    // vertical prestige biases
    offspring.vp[0] = mutation(
            mother.vp[allele_sample(rng_r)],
            mu_vp,
            sdmu_b);

    clamp(offspring.vp[0], bmin, bmax);

    offspring.vp[1] = mutation(
            father.vp[allele_sample(rng_r)],
            mu_vp,
            sdmu_b);

    clamp(offspring.vp[1], bmin, bmax);

    double vp_phen = 0.5 * (offspring.vp[0] + offspring.vp[1]);
    
    // inheritance of rel. sensitivity to 
    // vertical conformity biases
    offspring.vc[0] = mutation(
            mother.vc[allele_sample(rng_r)],
            mu_vc,
            sdmu_b);

    clamp(offspring.vc[0], bmin, bmax);

    offspring.vc[1] = mutation(
            father.vc[allele_sample(rng_r)],
            mu_vc,
            sdmu_b);

    clamp(offspring.vc[1], bmin, bmax);

    double vc_phen = 0.5 * (offspring.vc[0] + offspring.vc[1]);
    
    
    // inheritance of rel. sensitivity to 
    // horizontal prestige biases
    offspring.hp[0] = mutation(
            mother.hp[allele_sample(rng_r)],
            mu_hp,
            sdmu_b);

    clamp(offspring.hp[0], bmin, bmax);

    offspring.hp[1] = mutation(
            father.hp[allele_sample(rng_r)],
            mu_hp,
            sdmu_b);

    clamp(offspring.hp[1], bmin, bmax);

    // inheritance of rel. sensitivity to 
    // horizontal conformity biases
    offspring.hc[0] = mutation(
            mother.hc[allele_sample(rng_r)],
            mu_hc,
            sdmu_b);

    clamp(offspring.hc[0], bmin, bmax);

    offspring.hc[1] = mutation(
            father.hc[allele_sample(rng_r)],
            mu_hc,
            sdmu_b);

    clamp(offspring.hc[1], bmin, bmax);

    // kid receives juvenile cue
    offspring.cue_juv_envt_high = uniform(rng_r) < qjuv ? 
        offspring_envt_high : !offspring_envt_high;

    // get a cue of the maternal phenotype
    offspring.phen_mat = mother.phen_ad;

    // add noise to the maternal phenotypic cue
    offspring.phen_mat_error = offspring.phen_mat 
        + uniform(rng_r) * max_error_mat_phen;

    clamp(offspring.phen_mat_error, 0, 1.0);

    // get a cue of the maternal environment
    offspring.maternal_envt_cue = mother.cue_ad_envt_high;
    
    // add noise to the maternal phenotypic cue
    offspring.maternal_envt_cue_error = offspring.maternal_envt_cue;

    // express sensitivity to maternal phenotype
    double b_phen = 0.5 * (offspring.bmat_phen[0] + offspring.bmat_phen[1]);
    assert(b_phen >= bmin);
    assert(b_phen <= bmax);

    // express sensitivity to maternal environment 
    double b_envt = 0.5 * (offspring.bmat_envt[0] + offspring.bmat_envt[1]);
    assert(b_envt >= bmin);
    assert(b_phen <= bmax);

//    // generate maternal cue
//    double xoff = 1.0 /
//        (1.0 + exp(
//                   -b_phen * (mother.phen_ad - 0.5) 
//                   + 
//                   dmat_weighting * b_envt));
//
//    // generate maternal cue again 
//    // but then holding the maternal environment constant
//    // this is for stats purposes
//    double xoff_phen_only = 1.0 /
//        (1.0 + exp(
//                   -b_phen * (mother.phen_ad - 0.5)));
//
//    // generate maternal cue again 
//    // but then holding the maternal phenotype constant
//    // this is for stats purposes
//    double xoff_envt_only = 1.0 /
//        (1.0 + exp(dmat_weighting * b_envt));

//    // social learning
    offspring.xconformist_vert = xconformist_vert;
    offspring.xconformist_vert_error = xconformist_vert + 
        uniform(rng_r) * max_error_conform_vert;

    clamp(offspring.xconformist_vert_error, 0.0, 1.0);

    offspring.phen_prestige_vert = phen_prestige_vert;
    offspring.phen_prestige_vert_error = phen_prestige_vert +
        uniform(rng_r) * max_error_prestige_vert;
    
    clamp(offspring.phen_prestige_vert_error, 0.0, 1.0);
//
    offspring.phen_juv_logistic = 
                    aintercept_phen
                    + b_phen * (offspring.phen_mat_error - 0.5)
                    + b_envt * (offspring.maternal_envt_cue_error - 0.5)
                    + agen_phen * sum_genes 
                    + ajuv_phen * (offspring.cue_juv_envt_high - 0.5)
                    + vc_phen * (offspring.xconformist_vert_error - 0.5)
                    + vp_phen * (offspring.phen_prestige_vert_error - 0.5);
    
    // expressing a juvenile phenotype
    offspring.phen_juv = 1.0 / (1.0 + exp( - offspring.phen_juv_logistic ));

    offspring.phen_ad = NAN;
} // end create_offspring()

void adult_survival()
{
    // set up a random number generator to 
    // sample from the remaining breeders
    std::uniform_int_distribution<> random_patch(
            0,
            NPatches - 1);


    // reset survival statistics
    mean_survival[0] = 0.0;
    var_survival[0] = 0.0;
    mean_survival[1] = 0.0;
    var_survival[1] = 0.0;

    double surv = 0.0;
   
    // remember total number of breeders
    // in low and high patches
    // to calc averages of survival probabilities
    int n_total_individuals[2] = {0,0};

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // update counts of breeders in this enviroment
        n_total_individuals[Pop[patch_i].envt_high] += 
            Pop[patch_i].n_breeders;

        // breeders endure survival selection
        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            // check whether adult phenotypes indeed exist
            assert(
                    std::isnormal(
                        Pop[patch_i].breeders[breeder_i].phen_ad) > 0);

            // calculate survival probability
            surv = survival_probability(
                        Pop[patch_i].breeders[breeder_i].phen_ad
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

                // reduce for-loop counter
                --breeder_i;
                // reduce total number of adult breeders
                --Pop[patch_i].n_breeders;
            }
        } // end for (int breeder_i
    } // end for int patch_i

    // finalize survival statistics
    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        mean_survival[envt_i] /= n_total_individuals[envt_i];

        var_survival[envt_i] = 
            var_survival[envt_i] / n_total_individuals[envt_i]
            - mean_survival[envt_i] * mean_survival[envt_i];
    }
} // end survive_reproduce()

void social_learning(
        int const patch_i  // the patch in which social learning takes place
        ,bool const learning_is_horizontal // whether learning is horizontal or vertical
        ,double &prestige_phen // returning the prestige-bias phenotype
        ,double &xconformist) // returning the conformism-bias cue
{
    double surv, phen;

    int np_local, nc_local;

    int np = learning_is_horizontal ? nph : npv;
    int nc = learning_is_horizontal ? nch : ncv;

    // set up a random number generator to 
    // sample and socially learn from 
    // the surviving breeders
    std::uniform_int_distribution<> random_local_breeder(
            0,
            Pop[patch_i].n_breeders - 1);

    // check whether sampling np is feasible
    np_local = np > Pop[patch_i].n_breeders ? Pop[patch_i].n_breeders : np;

    // keep track on the current phenotype and survival value
    // of the individual with the highest 'prestige'
    //
    // if no individuals are sampled to learn about prestige
    // cue value is always 0.0. 
    prestige_phen = 0.0;
    double prestige_surv = 0.0;

    // compare survivorship values from np breeders and pick the highest
    for (int prest_i = 0; prest_i < np_local; ++prest_i)
    {
        // if learning is horizontal learn from other individual's juvenile 
        // phenotypes, otherwise learn from parental generation adult phenotypes
        phen = learning_is_horizontal ?
            Pop[patch_i].breeders[random_local_breeder(rng_r)].phen_juv
            :
            Pop[patch_i].breeders[random_local_breeder(rng_r)].phen_ad;

        assert(std::isnormal(phen) > 0);

        surv = survival_probability(phen,Pop[patch_i].envt_high);

        // update prestige bias
        if (surv > prestige_surv)
        {
            prestige_phen = phen;

            prestige_surv = surv;
        }
    }

    int nlo = 0;
    int nhi = 0;

    // check whether sampling nc is feasible
    nc_local = nc > Pop[patch_i].n_breeders ? Pop[patch_i].n_breeders : nc;

    // check which phenotype is most common, those below 0.5 or those above
    for (int conform_i = 0; conform_i < nc_local; ++conform_i)
    {
        phen = learning_is_horizontal ?
            Pop[patch_i].breeders[random_local_breeder(rng_r)].phen_juv
            :
            Pop[patch_i].breeders[random_local_breeder(rng_r)].phen_ad;

        assert(std::isnormal(phen) > 0);

        if (phen > 0.5)
        {
            ++nhi;
        }
        else if (phen <= 0.5)
        {
            ++nlo;
        }
    }

    // give conformist cue which is between 0 and 1
    xconformist = nhi/(nhi + nlo);
} // end void social_learning
 


// births of new offspring
// juvenile cue integration
// juvenile selection
// adult cue integration (aka horizontal social learning)
void replace()
{
    // randomly chosen remote patch to obtain
    // individuals from
    int random_remote_patch;

    // auxiliary variables to calculate
    // socially learnt cues through prestige-based social learning
    // and through conformism-based social learning
    double prestige_phen, xconformist;

    // set up a random number generator to 
    // sample remote patches
    std::uniform_int_distribution<> patch_sampler(
            0,
            NPatches - 1);


    // auxiliary variable to keep track 
    // of the survival probability
    double surv;

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        Pop[patch_i].n_recruits = 0;

        // now make NBreeder offspring and have them survive
        for (int offspring_i = 0; offspring_i < NBreeder; ++offspring_i)
        {
            Individual Kid;

            // offspring is born in local patch hence sample from local parents
            if (uniform(rng_r) < 1.0 - m 
                    &&
                    Pop[patch_i].n_breeders > 0)
            {
                // set up a random number generator to 
                // sample from the remaining breeders
                std::uniform_int_distribution<> random_local_breeder(
                        0,
                        Pop[patch_i].n_breeders - 1);

                // prepare cues in local patch from 
                // vertical social learning (i.e., from parents)
                social_learning(
                        patch_i
                        ,false
                        ,prestige_phen
                        ,xconformist);

                create_offspring(
                        Pop[patch_i].breeders[random_local_breeder(rng_r)]
                        ,Pop[patch_i].breeders[random_local_breeder(rng_r)]
                        ,Kid
                        ,Pop[patch_i].envt_high
                        ,prestige_phen
                        ,xconformist
                );
            }
            else // offspring born in remote patch
            {
                // sample a random remote pathc
                do {

                    random_remote_patch = patch_sampler(rng_r);

                }
                while(Pop[random_remote_patch].n_breeders < 1);
        
                std::uniform_int_distribution<> random_remote_breeder(
                0,
                Pop[random_remote_patch].n_breeders - 1);
               
                // prepare cues in remote patch (where parents are breeding)
                // from vertical social learning
                social_learning(
                        random_remote_patch
                        ,false
                        ,prestige_phen
                        ,xconformist);

                create_offspring(
                        Pop[random_remote_patch].breeders[random_remote_breeder(rng_r)]
                        ,Pop[random_remote_patch].breeders[random_remote_breeder(rng_r)]
                        ,Kid
                        ,Pop[random_remote_patch].envt_high
                        ,prestige_phen
                        ,xconformist
                );
            } // end else offsprign born in remote patch

            // in case juvenile selection acts,
            // and offspring dies, just continue on the next
            // iteration of this loop without adding the offspring
            // to the breeder population
            if (juvenile_survival)
            {
                surv = survival_probability(
                            Kid.phen_juv
                            ,Pop[patch_i].envt_high);

                // offspring dies 
                if (uniform(rng_r) > surv)
                {
                    continue;
                }
            }

            // offspring survives
            Pop[patch_i].breeders_t1[Pop[patch_i].n_recruits] = Kid;
            assert((int)Pop[patch_i].breeders_t1[Pop[patch_i].n_recruits].g[0].size() == nloci_g);

            ++Pop[patch_i].n_recruits;

            assert(Pop[patch_i].n_recruits <= NBreeder);

        } // for (int offspring_i = 0; offspring_i < NBreeder; ++offspring_i)
    } // end for (int patch_i = 0

    bool cue_ad_envt_high;

    // auxiliary variables for horizontal social learning
    double hp, hc;

    // all new breeders born etc, copy them over
    // and perform horiz social learning
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        // calculate adult cue value supplied to mothers
        // the cue value is the same for all mothers
        cue_ad_envt_high = uniform(rng_r) < qmat ? 
            Pop[patch_i].envt_high 
            : 
            !Pop[patch_i].envt_high;
        
        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_recruits; ++breeder_i)
        {
            Pop[patch_i].breeders[breeder_i] = 
                Pop[patch_i].breeders_t1[breeder_i];

            assert((int)Pop[patch_i].breeders[breeder_i].g[0].size() == nloci_g);
            
            // give breeder an environmental cue as adult
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_ad_envt_high;
        } // all juveniles are now assigned a breeding position
        
        Pop[patch_i].n_breeders = Pop[patch_i].n_recruits;

        // horizontal social learning 
        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            // ad_phen should be NaN as it is yet to be set
            assert(std::isnan(Pop[patch_i].breeders[breeder_i].phen_ad) > 0);

            // horizontal social learning takes place
            social_learning(
                    patch_i
                    ,true
                    ,prestige_phen
                    ,xconformist);

            hp = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hp[0]
                    + 
                    Pop[patch_i].breeders[breeder_i].hp[1]);
            
            hc = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hc[0]
                    + 
                    Pop[patch_i].breeders[breeder_i].hc[1]);

            Pop[patch_i].breeders[breeder_i].xconformist_horiz = xconformist;
            Pop[patch_i].breeders[breeder_i].xconformist_horiz_error = xconformist + 
                uniform(rng_r) * max_error_conform_horiz;

            clamp(Pop[patch_i].breeders[breeder_i].xconformist_horiz_error, 0.0, 1.0);

            Pop[patch_i].breeders[breeder_i].phen_prestige_horiz = prestige_phen;
            Pop[patch_i].breeders[breeder_i].phen_prestige_horiz_error = prestige_phen +
                uniform(rng_r) * max_error_prestige_horiz;

            clamp(Pop[patch_i].breeders[breeder_i].phen_prestige_horiz_error, 0.0, 1.0);

            Pop[patch_i].breeders[breeder_i].phen_ad_logistic = 
                Pop[patch_i].breeders[breeder_i].phen_juv_logistic
                + hp * (Pop[patch_i].breeders[breeder_i].phen_prestige_horiz_error - 0.5)
                + hc * (Pop[patch_i].breeders[breeder_i].xconformist_horiz_error - 0.5);

            // add social learning to the adult phenotype
            Pop[patch_i].breeders[breeder_i].phen_ad = 1.0 / 
                (1.0 + exp(-Pop[patch_i].breeders[breeder_i].phen_ad_logistic));
        } // for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        
        // envtal change after breeder establishment
        if (uniform(rng_r) < 1.0 - p)
        {
            Pop[patch_i].envt_high = !Pop[patch_i].envt_high;
        }

    } // end for (int patch_i = 0
}

// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    // get command line arguments
    init_arguments(argc, argv);

    std::string filename_stats = base_name + ".csv";
    std::ofstream DataFile(filename_stats.c_str());  // output file 

    // output file to write out the complete 
    // trait and state distribution of individuals
    // in the last generation
    std::string filename_dist = base_name + "_dist.csv";
    std::ofstream DataFileDist(filename_dist.c_str());  

    // write headers to the datafile
    write_data_headers(DataFile);
    
    write_data_headers_dist(DataFileDist);

    // initialize the population
    init_population();

    // auxiliary variable to store current generation
    int generation;

    for (generation = 0; generation < number_generations; ++generation)
    {
        // survival of adult breeders followed by reproduction
        adult_survival();

        // new breeder establishment, 
        // followed by horizontal learning
        // and finally enviromental change
        replace();

        if (generation % data_nth_generation == 0)
        {
            write_stats(DataFile, generation, 2);
        }
    }
            
    write_stats(DataFile, generation, 2);

    write_dist(DataFileDist);
    
    write_parameters(DataFile);
}
