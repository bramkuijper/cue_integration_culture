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

// the individual class, which defines properties of each
// individual used here
#include "individual.hpp"

#define DEBUG

// standard namespace
using namespace std;

// C++ random number generation
int seed = get_nanoseconds();
mt19937 rng_r{static_cast<long unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);

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

// noise in maternal cue 
double sdmat = 0.0;

// noise in social cue which is obtained 
// vertically and/or horizontally
double sdsoc_vert = 0.0;
double sdsoc_horiz = 0.0;

// whether survival selection has a sigmoidal or a quadratic function
bool sigmoidal_survival = true;

bool laplace = true;

// number of loci underlying the genetic cue
int nloci_g = 3;

// initial values
// these will be overridden when the function
// init_arguments is called
double init_g = 0.0;
double init_amat = 0.0;
double init_ajuv = 0.0;
double init_agen = 0.0;
double init_asoc_vert = 0.0;
double init_asoc_horiz = 0.0;
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
double mu_amat = 0.0;
double mu_ajuv = 0.0;
double mu_agen = 0.0;
double mu_asoc_horiz = 0.0;
double mu_asoc_vert = 0.0;
double mu_bmat_phen = 0.0;
double mu_bmat_envt = 0.0;
double mu_hp = 0.0;
double mu_hc = 0.0;
double mu_vp = 0.0;
double mu_vc = 0.0;

double sdmu_a = 0.0;
double sdmu_b = 0.0;
double sdmu_g = 0.0;

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
    nloci_g = atoi(argv[8]);
    init_g = atof(argv[9]);
    init_amat = atof(argv[10]);
    init_ajuv = atof(argv[11]);
    init_agen = atof(argv[12]);
    init_asoc_horiz = atof(argv[13]);
    init_asoc_vert = atof(argv[14]);
    init_bmat_phen = atof(argv[15]);
    init_bmat_envt = atof(argv[16]);
    init_hp = atof(argv[17]);
    init_hc = atof(argv[18]);
    init_vp = atof(argv[19]);
    init_vc = atof(argv[20]);

    gmin = atof(argv[21]);
    gmax = atof(argv[22]);
    amin = atof(argv[23]);
    amax = atof(argv[24]);
    bmin = atof(argv[25]);
    bmax = atof(argv[26]);
    sdmat = atof(argv[27]);
    sdsoc_vert = atof(argv[28]);
    sdsoc_horiz = atof(argv[29]);

    mu_g = atof(argv[30]);
    mu_amat = atof(argv[31]);
    mu_ajuv = atof(argv[32]);
    mu_agen = atof(argv[33]);
    mu_asoc_horiz = atof(argv[34]);
    mu_asoc_vert = atof(argv[35]);
    mu_bmat_phen = atof(argv[36]);
    mu_bmat_envt = atof(argv[37]);
    mu_hp = atof(argv[38]);
    mu_hc = atof(argv[39]);
    mu_vp = atof(argv[40]);
    mu_vc = atof(argv[41]);
    sdmu_a = atof(argv[42]);
    sdmu_b = atof(argv[43]);
    sdmu_g = atof(argv[44]);
    m = atof(argv[45]);
    nph = atoi(argv[46]);
    nch = atoi(argv[47]);
    npv = atoi(argv[48]);
    ncv = atoi(argv[49]);
}

// write down all parameters to the file DataFile
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
        << "init_asoc_vert;" << init_asoc_vert << ";"<< endl
        << "init_asoc_horiz;" << init_asoc_horiz << ";"<< endl
        << "init_bmat_phen;" << init_bmat_phen << ";"<< endl
        << "init_bmat_envt;" << init_bmat_envt << ";"<< endl
        << "init_hp;" << init_hp << ";"<< endl
        << "init_hc;" << init_hc << ";"<< endl
        << "init_vp;" << init_vp << ";"<< endl
        << "init_vc;" << init_vc << ";"<< endl
        << "gmin;" << gmin << ";"<< endl
        << "gmax;" << gmax << ";"<< endl
        << "amin;" << amin << ";"<< endl
        << "amax;" << amax << ";"<< endl
        << "bmin;" << bmin << ";"<< endl
        << "bmax;" << bmax << ";"<< endl
        << "sdmat;" << sdmat << ";"<< endl
        << "sdsoc_vert;" << sdsoc_vert << ";"<< endl
        << "sdsoc_horiz;" << sdsoc_horiz << ";"<< endl
        << "mu_g;" << mu_g << ";"<< endl
        << "sdmu_g;" << sdmu_g << ";"<< endl
        << "mu_amat;" << mu_amat << ";"<< endl
        << "mu_ajuv;" << mu_ajuv << ";"<< endl
        << "mu_agen;" << mu_agen << ";"<< endl
        << "mu_asoc_horiz;" << mu_asoc_horiz << ";"<< endl
        << "mu_asoc_vert;" << mu_asoc_vert << ";"<< endl
        << "mu_bmat_phen;" << mu_bmat_phen << ";"<< endl
        << "mu_bmat_envt;" << mu_bmat_envt << ";"<< endl
        << "mu_vc;" << mu_vc << ";"<< endl
        << "mu_vp;" << mu_vp << ";"<< endl
        << "mu_hc;" << mu_hc << ";"<< endl
        << "mu_hp;" << mu_hp << ";"<< endl
        << "sdmu_a;" << sdmu_a << ";"<< endl
        << "sdmu_b;" << sdmu_b << ";"<< endl
        << "m;" << m << ";"<< endl
        << "nph;" << nph << ";"<< endl
        << "nch;" << nch << ";"<< endl
        << "npv;" << npv << ";"<< endl
        << "ncv;" << ncv << ";"<< endl
        << "survival_scalar0;" << survival_scalar[0] << ";"<< endl
        << "survival_scalar1;" << survival_scalar[1] << ";"<< endl
        << "seed;" << seed << ";" << endl;

} // void write_parameters(ofstream &DataFile)


// write all properties of all individuals
// to the file DataFile (to obtain information about
// the distribution of traits)
void write_dist(ofstream &DataFile)
{
    double g; // auxiliary variable to temporarily 
                // store trait expression

    
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            DataFile << patch_i << ";" 
                << breeder_i << ";"
                << Pop[patch_i].breeders[breeder_i].phen_ad << ";"
                << Pop[patch_i].breeders[breeder_i].phen_mat << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_vert << ";"
                << Pop[patch_i].breeders[breeder_i].phen_prestige_horiz << ";"
                << Pop[patch_i].breeders[breeder_i].xmat << ";"
                << Pop[patch_i].breeders[breeder_i].xsoc_vert << ";"
                << Pop[patch_i].breeders[breeder_i].xsoc_horiz << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_vert << ";"
                << Pop[patch_i].breeders[breeder_i].xconformist_horiz << ";"

                // agen
                << 0.5 * (Pop[patch_i].breeders[breeder_i].agen[0]
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1]) << ";"

                // ajuv
                << 0.5 * (Pop[patch_i].breeders[breeder_i].ajuv[0]
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1]) << ";"

                // amat
                << 0.5 * (Pop[patch_i].breeders[breeder_i].amat[0]
                    +
                    Pop[patch_i].breeders[breeder_i].amat[1]) << ";"

                // asoc_vert 
                << 0.5 * (Pop[patch_i].breeders[breeder_i].asoc_vert[0]
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_vert[1]) << ";"
                
                // asoc_horiz 
                << 0.5 * (Pop[patch_i].breeders[breeder_i].asoc_horiz[0]
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[1]) << ";"

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
                << Pop[patch_i].breeders[breeder_i].maternal_cue << ";"
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
        << "phen_ad;" 
        << "phen_mat;" 
        << "phen_prestige_vert;" 
        << "phen_prestige_horiz;" 
        << "xmat;" 
        << "xsoc_vert;" 
        << "xsoc_horiz;" 
        << "xconformist_vert;" 
        << "xconformist_horiz;" 
        << "agen;" 
        << "ajuv;" 
        << "amat;" 
        << "asoc_vert;" 
        << "asoc_horiz;" 
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
        << "maternal_cue;" 
        << endl; 
}

// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile 
        << "generation;" 
        << "mean_phen_ad;" 
        << "mean_phen_juv;" 
        << "mean_phen_prestige_vert;" 
        << "mean_phen_prestige_horiz;" 
        << "mean_agen;" 
        << "mean_ajuv;" 
        << "mean_amat;" 
        << "mean_asoc_vert;" 
        << "mean_asoc_horiz;" 
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
        << "var_agen;" 
        << "var_ajuv;" 
        << "var_amat;" 
        << "var_asoc_vert;" 
        << "var_asoc_horiz;" 
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
        << "var_surv1;" 
        << "var_component_gen;"
        << "var_component_ajuv;"
        << "var_component_amat;"
        << "var_component_amat_envt;"
        << "var_component_amat_phen;"
        << "cov_amat_ajuv;"
        << "cov_amat_envt_ajuv;"
        << "cov_amat_phen_ajuv;"
        << "var_component_asoc_vert;"
        << "var_component_asoc_vert_p;"
        << "var_component_asoc_vert_c;"
        << "var_component_asoc_horiz;"
        << "var_component_asoc_horiz_p;"
        << "var_component_asoc_horiz_c;"
        << "cov_amat_asoc_vert;"
        << "cov_amat_asoc_vert_c;"
        << "cov_amat_asoc_vert_p;"
        << "cov_amat_asoc_horiz;"
        << "cov_amat_asoc_horiz_c;"
        << "cov_amat_asoc_horiz_p;"
        << "cov_amat_envt_asoc_vert;"
        << "cov_amat_envt_asoc_vert_c;"
        << "cov_amat_envt_asoc_vert_p;"
        << "cov_amat_envt_asoc_horiz;"
        << "cov_amat_envt_asoc_horiz_c;"
        << "cov_amat_envt_asoc_horiz_p;"
        << "cov_amat_phen_asoc_vert;"
        << "cov_amat_phen_asoc_vert_c;"
        << "cov_amat_phen_asoc_vert_p;"
        << "cov_amat_phen_asoc_horiz;"
        << "cov_amat_phen_asoc_horiz_c;"
        << "cov_amat_phen_asoc_horiz_p;"

        << "cov_ajuv_asoc_vert;"
        << "cov_ajuv_asoc_vert_c;"
        << "cov_ajuv_asoc_vert_p;"
        
        << "cov_ajuv_asoc_horiz;"
        << "cov_ajuv_asoc_horiz_c;"
        << "cov_ajuv_asoc_horiz_p;"

        << "cov_agen_asoc_vert;"
        << "cov_agen_asoc_vert_c;"
        << "cov_agen_asoc_vert_p;"

        << "cov_agen_asoc_horiz;"
        << "cov_agen_asoc_horiz_c;"
        << "cov_agen_asoc_horiz_p;"

        << "cov_agen_ajuv;"
        << endl; 
}

// write data both for winter and summer populations
void write_stats(ofstream &DataFile, int generation, int timestep)
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
   
    double mean_agen = 0.0;
    double ss_agen = 0.0;
   
    double mean_amat = 0.0;
    double ss_amat = 0.0;
    
    double mean_ajuv = 0.0;
    double ss_ajuv = 0.0;
    
    double mean_asoc_vert = 0.0;
    double ss_asoc_vert = 0.0;
    
    double mean_asoc_horiz = 0.0;
    double ss_asoc_horiz = 0.0;
    
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

    // now some stats at the logit scale
    // maternal variance relative to total variance
    double ss1_gen_component = 0.0;
    double ss2_gen_component = 0.0;

    double ss1_ajuv_component = 0.0;
    double ss2_ajuv_component = 0.0;

    double ss1_amat_component = 0.0;
    double ss2_amat_component = 0.0;

    double ss1_amat_envt_component = 0.0;
    double ss2_amat_envt_component = 0.0;

    double ss1_amat_phen_component = 0.0;
    double ss2_amat_phen_component = 0.0;

    double ss1_asoc_vert_component = 0.0;
    double ss2_asoc_vert_component = 0.0;
    
    double ss1_asoc_horiz_component = 0.0;
    double ss2_asoc_horiz_component = 0.0;

    double ss1_asoc_vert_c_component = 0.0;
    double ss2_asoc_vert_c_component = 0.0;
    double ss1_asoc_vert_p_component = 0.0;
    double ss2_asoc_vert_p_component = 0.0;
    
    double ss1_asoc_horiz_c_component = 0.0;
    double ss2_asoc_horiz_c_component = 0.0;
    double ss1_asoc_horiz_p_component = 0.0;
    double ss2_asoc_horiz_p_component = 0.0;

    double ss2_cov_agen_ajuv = 0.0;

    double ss2_cov_amat_ajuv = 0.0;
    double ss2_cov_amat_envt_ajuv = 0.0;
    double ss2_cov_amat_phen_ajuv = 0.0;

    double ss2_cov_agen_asoc_vert = 0.0;
    double ss2_cov_agen_asoc_vert_c = 0.0;
    double ss2_cov_agen_asoc_vert_p = 0.0;

    double ss2_cov_ajuv_asoc_vert = 0.0;
    double ss2_cov_ajuv_asoc_vert_c = 0.0;
    double ss2_cov_ajuv_asoc_vert_p = 0.0;
    
    double ss2_cov_ajuv_asoc_horiz = 0.0;
    double ss2_cov_ajuv_asoc_horiz_c = 0.0;
    double ss2_cov_ajuv_asoc_horiz_p = 0.0;

    double ss2_cov_amat_asoc_vert = 0.0;
    double ss2_cov_amat_asoc_vert_c = 0.0;
    double ss2_cov_amat_asoc_vert_p = 0.0;
    
    double ss2_cov_amat_asoc_horiz = 0.0;
    double ss2_cov_amat_asoc_horiz_c = 0.0;
    double ss2_cov_amat_asoc_horiz_p = 0.0;

    double ss2_cov_amat_phen_asoc_vert = 0.0;
    double ss2_cov_amat_phen_asoc_vert_c = 0.0;
    double ss2_cov_amat_phen_asoc_vert_p = 0.0;
    
    double ss2_cov_amat_envt_asoc_vert = 0.0;
    double ss2_cov_amat_envt_asoc_vert_c = 0.0;
    double ss2_cov_amat_envt_asoc_vert_p = 0.0;
    
    double ss2_cov_amat_phen_asoc_horiz = 0.0;
    double ss2_cov_amat_phen_asoc_horiz_c = 0.0;
    double ss2_cov_amat_phen_asoc_horiz_p = 0.0;
    
    double ss2_cov_amat_envt_asoc_horiz = 0.0;
    double ss2_cov_amat_envt_asoc_horiz_c = 0.0;
    double ss2_cov_amat_envt_asoc_horiz_p = 0.0;

    double ss2_cov_agen_asoc_horiz = 0.0;
    double ss2_cov_agen_asoc_horiz_c = 0.0;
    double ss2_cov_agen_asoc_horiz_p = 0.0;

    // auxiliary variables to calculate an individual's phenotype
    double g, z, agen, amat, ajuv, xmat_envt, xmat, xmat_phen, xsoc_vert, xsoc_vert_p, xsoc_vert_c, asoc_horiz, asoc_vert, cue_ad, xsoc_horiz_p, xsoc_horiz_c, xsoc_horiz;


    // summing means and sums of squares over all patches and breeders
    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        freq_high += Pop[patch_i].envt_high;

        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            g = 0.0;
            for (int g_i = 0; g_i < nloci_g; ++g_i)
            {
                g +=  0.5 * (
                        Pop[patch_i].breeders[breeder_i].g[0][g_i]
                        +
                        Pop[patch_i].breeders[breeder_i].g[1][g_i]
                        );

            } // end for (int g_i = 0; g_i < nloci_g; ++g_i)

            mean_g += g;
            ss_g += g * g;

            // adult phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_ad;
            mean_phen_ad += z;
            ss_phen_ad += z * z;
            
            // juvenile phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_juv;
            mean_phen_juv += z;
            ss_phen_juv += z * z;
            
            // prestige phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_prestige_vert;
            mean_phen_prestige_vert += z;
            ss_phen_prestige_vert += z * z;
            
            // prestige phenotype
            z = Pop[patch_i].breeders[breeder_i].phen_prestige_horiz;
            mean_phen_prestige_horiz += z;
            ss_phen_prestige_horiz += z * z;

            // sensitivity to genetic cues
            agen = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].agen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].agen[1] 
                    );

            mean_agen += agen;
            ss_agen += agen * agen;

            // variance components (liability scale)
            // sum of squares for the 
            // genetic variance component is
            // var(agen * sum(g)) = E[agen^2 * sum(g)^2] - E[agen * sum(g)]^2
            ss1_gen_component += agen * g;
            ss2_gen_component += agen * agen * g * g;

            // sensitivity to maternal cues
            amat = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].amat[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].amat[1] 
                    );

            mean_amat += amat;
            ss_amat += amat * amat;

            xmat_envt = Pop[patch_i].breeders[breeder_i].xmat_envt_only;
            xmat_phen  = Pop[patch_i].breeders[breeder_i].xmat_phen_only;
            xmat = Pop[patch_i].breeders[breeder_i].xmat;

            // variance components (liability scale)
            // sum of squares for the 
            // maternal environmental component
            ss1_amat_envt_component += amat * xmat_envt;
            ss2_amat_envt_component += amat * amat * xmat_envt * xmat_envt;

            // variance components (liability scale)
            // sum of squares for the 
            // maternal phenotypic component
            ss1_amat_phen_component += amat * xmat_phen;
            ss2_amat_phen_component += amat * amat * xmat_phen * xmat_phen;

            // variance components (liability scale)
            // sum of squares for the 
            // total maternal component
            ss1_amat_component += amat * xmat;
            ss2_amat_component += amat * amat * xmat * xmat;

            // sensitivity to juvenile cues
            ajuv = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].ajuv[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].ajuv[1] 
                    );
            
            mean_ajuv += ajuv;
            ss_ajuv += ajuv * ajuv;

            cue_ad = Pop[patch_i].breeders[breeder_i].cue_ad_envt_high;

            // variance components juvenile cue
            ss1_ajuv_component += ajuv * cue_ad;
            ss2_ajuv_component += ajuv * ajuv * cue_ad * cue_ad;

            // covariance between juvenile cue and maternal total effect
            // E[abcd] - E[ab]E[cd] (the latter product was what we already
            // calculated before)
            ss2_cov_amat_ajuv += amat * xmat * ajuv * cue_ad;
            ss2_cov_amat_envt_ajuv += amat * xmat_envt * ajuv * cue_ad;
            ss2_cov_amat_phen_ajuv += amat * xmat_phen * ajuv * cue_ad;

            ss2_cov_agen_ajuv += ajuv * cue_ad * agen * g;

            // sensitivity to socially learnt cues (vertically)
            asoc_vert = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].asoc_vert[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_vert[1] 
                    );

            mean_asoc_vert += asoc_vert;
            ss_asoc_vert += asoc_vert * asoc_vert;

            xsoc_vert = Pop[patch_i].breeders[breeder_i].xsoc_vert;
            xsoc_vert_c = Pop[patch_i].breeders[breeder_i].xsoc_vert_c;
            xsoc_vert_p = Pop[patch_i].breeders[breeder_i].xsoc_vert_p;

            // variance due to the total socially learnt component
            ss1_asoc_vert_component += asoc_vert * xsoc_vert;
            ss2_asoc_vert_component += asoc_vert * asoc_vert * xsoc_vert * xsoc_vert;

            // variance due to the conformity bias socially learnt component
            ss1_asoc_vert_c_component += asoc_vert * xsoc_vert_c;
            ss2_asoc_vert_c_component += asoc_vert * asoc_vert * xsoc_vert_c * xsoc_vert_c;
            
            // variance due to the prestige bias socially learnt component
            ss1_asoc_vert_p_component += asoc_vert * xsoc_vert_p;
            ss2_asoc_vert_p_component += asoc_vert * asoc_vert * xsoc_vert_p * xsoc_vert_p;
            
            // covariance between maternal and vertically learnt components
            ss2_cov_amat_asoc_vert += amat * xmat * asoc_vert * xsoc_vert;
            ss2_cov_amat_asoc_vert_c += amat * xmat * asoc_vert * xsoc_vert_c;
            ss2_cov_amat_asoc_vert_p += amat * xmat * asoc_vert * xsoc_vert_p;

            // covariance between maternal environment and vertically learnt components
            ss2_cov_amat_envt_asoc_vert += amat * xmat_envt * asoc_vert * xsoc_vert;
            ss2_cov_amat_envt_asoc_vert_c += amat * xmat_envt * asoc_vert * xsoc_vert_c;
            ss2_cov_amat_envt_asoc_vert_p += amat * xmat_envt * asoc_vert * xsoc_vert_p;

            // covariance between maternal phenotype and vertically learnt components
            ss2_cov_amat_phen_asoc_vert += amat * xmat_phen * asoc_vert * xsoc_vert;
            ss2_cov_amat_phen_asoc_vert_c += amat * xmat_phen * asoc_vert * xsoc_vert_c;
            ss2_cov_amat_phen_asoc_vert_p += amat * xmat_phen * asoc_vert * xsoc_vert_p;
            
            // covariance between juvenile cue and vertically learnt components
            ss2_cov_ajuv_asoc_vert += ajuv * cue_ad * asoc_vert * xsoc_vert;
            ss2_cov_ajuv_asoc_vert_c += ajuv * cue_ad * asoc_vert * xsoc_vert_c;
            ss2_cov_ajuv_asoc_vert_p += ajuv * cue_ad * asoc_vert * xsoc_vert_p;

            // covariance between genetic cue and vertically learnt components
            ss2_cov_agen_asoc_vert += agen * g * asoc_vert * xsoc_vert;
            ss2_cov_agen_asoc_vert_c += agen * g * asoc_vert * xsoc_vert_c;
            ss2_cov_agen_asoc_vert_p += agen * g * asoc_vert * xsoc_vert_p;

            // sensitivity to socially learnt cues (horizontally)
            asoc_horiz = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[1] 
                    );

            mean_asoc_horiz += asoc_horiz;
            ss_asoc_horiz += asoc_horiz * asoc_horiz;

            xsoc_horiz = Pop[patch_i].breeders[breeder_i].xsoc_horiz;
            xsoc_horiz_c = Pop[patch_i].breeders[breeder_i].xsoc_horiz_c;
            xsoc_horiz_p = Pop[patch_i].breeders[breeder_i].xsoc_horiz_p;

            ss1_asoc_horiz_component += asoc_horiz * xsoc_horiz;
            ss2_asoc_horiz_component += asoc_horiz * asoc_horiz * xsoc_horiz * xsoc_horiz;

            ss1_asoc_horiz_c_component += asoc_horiz * xsoc_horiz_c;
            ss2_asoc_horiz_c_component += asoc_horiz * asoc_horiz * xsoc_horiz_c * xsoc_horiz_c;

            ss1_asoc_horiz_p_component += asoc_horiz * xsoc_horiz_p;
            ss2_asoc_horiz_p_component += asoc_horiz * asoc_horiz * xsoc_horiz_p * xsoc_horiz_p;

            // covariance between maternal and vertically learnt components
            ss2_cov_amat_asoc_horiz += amat * xmat * asoc_horiz * xsoc_horiz;
            ss2_cov_amat_asoc_horiz_c += amat * xmat * asoc_horiz * xsoc_horiz_c;
            ss2_cov_amat_asoc_horiz_p += amat * xmat * asoc_horiz * xsoc_horiz_p;
            
            // covariance between maternal environmental cue and vertically learnt components
            ss2_cov_amat_envt_asoc_horiz += amat * xmat_envt * asoc_horiz * xsoc_horiz;
            ss2_cov_amat_envt_asoc_horiz_c += amat * xmat_envt * asoc_horiz * xsoc_horiz_c;
            ss2_cov_amat_envt_asoc_horiz_p += amat * xmat_envt * asoc_horiz * xsoc_horiz_p;
            
            // covariance between maternal phenotypic cue and vertically learnt components
            ss2_cov_amat_phen_asoc_horiz += amat * xmat_phen * asoc_horiz * xsoc_horiz;
            ss2_cov_amat_phen_asoc_horiz_c += amat * xmat_phen * asoc_horiz * xsoc_horiz_c;
            ss2_cov_amat_phen_asoc_horiz_p += amat * xmat_phen * asoc_horiz * xsoc_horiz_p;
            
            // covariance between juvenile cue and horizically learnt components
            ss2_cov_ajuv_asoc_horiz += ajuv * cue_ad * asoc_horiz * xsoc_horiz;
            ss2_cov_ajuv_asoc_horiz_c += ajuv * cue_ad * asoc_horiz * xsoc_horiz_c;
            ss2_cov_ajuv_asoc_horiz_p += ajuv * cue_ad * asoc_horiz * xsoc_horiz_p;

            // covariance between genetic cue and horizically learnt components
            ss2_cov_agen_asoc_horiz += agen * g * asoc_horiz * xsoc_horiz;
            ss2_cov_agen_asoc_horiz_c += agen * g * asoc_horiz * xsoc_horiz_c;
            ss2_cov_agen_asoc_horiz_p += agen * g * asoc_horiz * xsoc_horiz_p;
            
            // maternal sensitivity to phenotypic cues
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_phen[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_phen[1] 
                    );

            mean_bmat_phen += z;
            ss_bmat_phen += z * z;
            
            // maternal sensitivity to environmental cues
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].bmat_envt[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].bmat_envt[1] 
                    );

            mean_bmat_envt += z;
            ss_bmat_envt += z * z;
            
            // sensitivity to performance-based cues when learning
            // horizontally
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hp[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].hp[1] 
                    );

            mean_hp += z;
            ss_hp += z * z;
            
            // sensitivity to confirmity-based cues when learning
            // horizontally
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].hc[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].hc[1] 
                    );

            mean_hc += z;
            ss_hc += z * z;
            
            
            // sensitivity to performance-based cues when learning
            // vertically 
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].vp[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].vp[1] 
                    );

            mean_vp += z;
            ss_vp += z * z;
            
            // sensitivity to confirmity-based cues when learning
            // vertically
            z = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].vc[0] 
                    +
                    Pop[patch_i].breeders[breeder_i].vc[1] 
                    );

            mean_vc += z;
            ss_vc += z * z;
        } // end for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
    } // end for (int patch_i = 0; patch_i < NPatches; ++patch_i)

    mean_phen_ad /= NPatches * NBreeder;
    double var_phen_ad = ss_phen_ad / (NPatches * NBreeder) - mean_phen_ad * mean_phen_ad;
    
    mean_phen_juv /= NPatches * NBreeder;
    double var_phen_juv = ss_phen_juv / (NPatches * NBreeder) - mean_phen_juv * mean_phen_juv;
    
    mean_phen_prestige_vert /= NPatches * NBreeder;

    double var_phen_prestige_vert = ss_phen_prestige_vert / (NPatches * NBreeder) 
        - mean_phen_prestige_vert * mean_phen_prestige_vert;
   
    mean_phen_prestige_horiz /= NPatches * NBreeder;

    double var_phen_prestige_horiz = ss_phen_prestige_horiz / (NPatches * NBreeder) 
        - mean_phen_prestige_horiz * mean_phen_prestige_horiz;

    mean_agen /= NPatches * NBreeder;
    double var_agen = ss_agen / (NPatches * NBreeder) - mean_agen * mean_agen;
   
    mean_amat /= NPatches * NBreeder;
    double var_amat = NPatches * NBreeder - mean_amat * mean_amat;
    
    mean_ajuv /= NPatches * NBreeder;
    double var_ajuv = ss_ajuv / (NPatches * NBreeder) -  mean_ajuv * mean_ajuv;
    
    mean_asoc_vert /= NPatches * NBreeder;
    double var_asoc_vert = ss_asoc_vert / (NPatches * NBreeder) - mean_asoc_vert * mean_asoc_vert;
    
    mean_asoc_horiz /= NPatches * NBreeder;
    double var_asoc_horiz = ss_asoc_horiz / (NPatches * NBreeder) - mean_asoc_horiz * mean_asoc_horiz;
    
    mean_bmat_phen /= NPatches * NBreeder;
    double var_bmat_phen = ss_bmat_phen / (NPatches * NBreeder) - mean_bmat_phen * mean_bmat_phen;
    
    mean_bmat_envt /= NPatches * NBreeder;
    double var_bmat_envt = ss_bmat_envt / (NPatches * NBreeder) - mean_bmat_envt * mean_bmat_envt;
    
    mean_hp /= NPatches * NBreeder;
    double var_hp = ss_hp / (NPatches * NBreeder) - mean_hp * mean_hp;
    
    mean_hc /= NPatches * NBreeder;
    double var_hc = ss_hc / (NPatches * NBreeder) - mean_hc * mean_hc;
    
    mean_vp /= NPatches * NBreeder;
    double var_vp = ss_vp / (NPatches * NBreeder) - mean_vp * mean_vp;
    
    mean_vc /= NPatches * NBreeder;
    double var_vc = ss_vc / (NPatches * NBreeder) - mean_vc * mean_vc;

    mean_g /= NPatches * NBreeder;
    double var_g = ss_g / (NPatches * NBreeder) - mean_g * mean_g;



    // Variance components at the logit scale

    double var_component_gen = ss2_gen_component / (NPatches * NBreeder) - 
        pow(ss1_gen_component / (NPatches * NBreeder),2);
    
    double var_component_ajuv  = ss2_ajuv_component / (NPatches * NBreeder) - 
        pow(ss1_ajuv_component / (NPatches * NBreeder),2);

    double var_component_amat  = ss2_amat_component / (NPatches * NBreeder) - 
        pow(ss1_amat_component / (NPatches * NBreeder),2);

    double var_component_amat_envt  = ss2_amat_envt_component / (NPatches * NBreeder) - 
        pow(ss1_amat_envt_component / (NPatches * NBreeder),2);
    
    double var_component_amat_phen  = ss2_amat_phen_component / (NPatches * NBreeder) - 
        pow(ss1_amat_phen_component / (NPatches * NBreeder),2);
    
    double cov_amat_ajuv = ss2_cov_amat_ajuv / (NPatches * NBreeder) - 
        ss1_amat_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);
    
    double cov_amat_envt_ajuv = ss2_cov_amat_envt_ajuv / (NPatches * NBreeder) - 
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);
    
    double cov_amat_phen_ajuv = ss2_cov_amat_phen_ajuv / (NPatches * NBreeder) - 
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);

    // variance of vertical social learning: total, conformism and prestige
    double var_component_asoc_vert = ss2_asoc_vert_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_vert_component / (NPatches * NBreeder),2);

    double var_component_asoc_vert_p = ss2_asoc_vert_p_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_vert_p_component / (NPatches * NBreeder),2);

    double var_component_asoc_vert_c = ss2_asoc_vert_c_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_vert_c_component / (NPatches * NBreeder),2);

    // variance of horizontal social learning: total, conformism and prestige
    double var_component_asoc_horiz = ss2_asoc_horiz_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_horiz_component / (NPatches * NBreeder),2);

    double var_component_asoc_horiz_p = ss2_asoc_horiz_p_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_horiz_p_component / (NPatches * NBreeder),2);

    double var_component_asoc_horiz_c = ss2_asoc_horiz_c_component / (NPatches * NBreeder) - 
        pow(ss1_asoc_horiz_c_component / (NPatches * NBreeder),2);



    
    // covariance between maternal total cue and vertical social learning
    double cov_amat_asoc_vert = ss2_cov_amat_asoc_vert / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    double cov_amat_asoc_vert_c = ss2_cov_amat_asoc_vert_c / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);
    
    double cov_amat_asoc_vert_p = ss2_cov_amat_asoc_vert_p / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);


    
    // covariance between maternal total cue and vertical social learning
    double cov_amat_asoc_horiz = ss2_cov_amat_asoc_horiz / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    double cov_amat_asoc_horiz_c = ss2_cov_amat_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);
    
    double cov_amat_asoc_horiz_p = ss2_cov_amat_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_amat_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);



    // covariance between maternal environmental cue and vertical social learning
    double cov_amat_envt_asoc_vert = ss2_cov_amat_envt_asoc_vert / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    double cov_amat_envt_asoc_vert_c = ss2_cov_amat_envt_asoc_vert_c / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);
    
    double cov_amat_envt_asoc_vert_p = ss2_cov_amat_envt_asoc_vert_p / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);


    // covariance between maternal environmental cue and horizontal social learning
    double cov_amat_envt_asoc_horiz = ss2_cov_amat_envt_asoc_horiz / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    double cov_amat_envt_asoc_horiz_c = ss2_cov_amat_envt_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);
    
    double cov_amat_envt_asoc_horiz_p = ss2_cov_amat_envt_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_amat_envt_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);


    // covariance between maternal phenotypic cue and vertical social learning
    double cov_amat_phen_asoc_vert = ss2_cov_amat_phen_asoc_vert / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    double cov_amat_phen_asoc_vert_c = ss2_cov_amat_phen_asoc_vert_c / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);
    
    double cov_amat_phen_asoc_vert_p = ss2_cov_amat_phen_asoc_vert_p / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);


    // covariance between maternal phenotypic cue and horizontal social learning
    double cov_amat_phen_asoc_horiz = ss2_cov_amat_phen_asoc_horiz / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    double cov_amat_phen_asoc_horiz_c = ss2_cov_amat_phen_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);
    
    double cov_amat_phen_asoc_horiz_p = ss2_cov_amat_phen_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_amat_phen_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);


    // covariance between vertical social learning and juvenile cues
    double cov_ajuv_asoc_vert = ss2_cov_ajuv_asoc_vert / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_vert_p = ss2_cov_ajuv_asoc_vert_p / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_vert_c = ss2_cov_ajuv_asoc_vert_c / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);

    
    // covariance between horizontal social learning and juvenile cues
    double cov_ajuv_asoc_horiz = ss2_cov_ajuv_asoc_horiz / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_horiz_p = ss2_cov_ajuv_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);
    
    double cov_ajuv_asoc_horiz_c = ss2_cov_ajuv_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_ajuv_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);


    // covariance between vertical social learning and genetic cues
    double cov_agen_asoc_vert = ss2_cov_agen_asoc_vert / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_vert_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_vert_p = ss2_cov_agen_asoc_vert_p / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_vert_p_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_vert_c = ss2_cov_agen_asoc_vert_c / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_vert_c_component / (NPatches * NBreeder);

    // covariance between horizontal social learning and genetic cues
    double cov_agen_asoc_horiz = ss2_cov_agen_asoc_horiz / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_horiz_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_horiz_p = ss2_cov_agen_asoc_horiz_p / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_horiz_p_component / (NPatches * NBreeder);
    
    double cov_agen_asoc_horiz_c = ss2_cov_agen_asoc_horiz_c / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_asoc_horiz_c_component / (NPatches * NBreeder);

    // covariance between genetic cue and juvenile cue 
    double cov_agen_ajuv = ss2_cov_agen_ajuv / (NPatches * NBreeder) -
        ss1_gen_component / (NPatches * NBreeder) * ss1_ajuv_component / (NPatches * NBreeder);

    freq_high /= NPatches;

    DataFile << generation << ";"
        << mean_phen_ad << ";"
        << mean_phen_juv << ";"
        << mean_phen_prestige_vert << ";"
        << mean_phen_prestige_horiz << ";"
        << mean_agen << ";"
        << mean_ajuv << ";"
        << mean_amat << ";"
        << mean_asoc_vert << ";"
        << mean_asoc_horiz << ";"
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
        << var_agen << ";"
        << var_ajuv << ";"
        << var_amat << ";"
        << var_asoc_vert << ";"
        << var_asoc_horiz << ";"
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
        << var_survival[1] << ";" 

        << var_component_gen << ";"
        << var_component_ajuv << ";"
        << var_component_amat << ";"
        << var_component_amat_envt << ";"
        << var_component_amat_phen << ";"

        << cov_amat_ajuv << ";"
        << cov_amat_envt_ajuv << ";"
        << cov_amat_phen_ajuv << ";"

        << var_component_asoc_vert << ";"
        << var_component_asoc_vert_p << ";"
        << var_component_asoc_vert_c << ";"

        << var_component_asoc_horiz << ";"
        << var_component_asoc_horiz_p << ";"
        << var_component_asoc_horiz_c << ";"

        << cov_amat_asoc_vert << ";"
        << cov_amat_asoc_vert_c << ";"
        << cov_amat_asoc_vert_p << ";"

        << cov_amat_asoc_horiz << ";"
        << cov_amat_asoc_horiz_c << ";"
        << cov_amat_asoc_horiz_p << ";"

        << cov_amat_envt_asoc_vert << ";"
        << cov_amat_envt_asoc_vert_c << ";"
        << cov_amat_envt_asoc_vert_p << ";"
        
        << cov_amat_envt_asoc_horiz << ";"
        << cov_amat_envt_asoc_horiz_c << ";"
        << cov_amat_envt_asoc_horiz_p << ";"

        << cov_amat_phen_asoc_vert << ";"
        << cov_amat_phen_asoc_vert_c << ";"
        << cov_amat_phen_asoc_vert_p << ";"

        << cov_amat_phen_asoc_horiz << ";"
        << cov_amat_phen_asoc_horiz_c << ";"
        << cov_amat_phen_asoc_horiz_p << ";"


        << cov_ajuv_asoc_vert << ";"
        << cov_ajuv_asoc_vert_c << ";"
        << cov_ajuv_asoc_vert_p << ";"

        << cov_ajuv_asoc_horiz << ";"
        << cov_ajuv_asoc_horiz_c << ";"
        << cov_ajuv_asoc_horiz_p << ";"

        << cov_agen_asoc_vert << ";"
        << cov_agen_asoc_vert_c << ";"
        << cov_agen_asoc_vert_p << ";"
        
        << cov_agen_asoc_horiz << ";"
        << cov_agen_asoc_horiz_c << ";"
        << cov_agen_asoc_horiz_p << ";"

        << cov_agen_ajuv << ";"
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
                
                // vertical social cue weighting
                Pop[patch_i].breeders[breeder_i].asoc_vert[allele_i] = init_asoc_vert;
                
                // horizontal social cue weighting
                Pop[patch_i].breeders[breeder_i].asoc_horiz[allele_i] = init_asoc_horiz;
                
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

            double cue_sum = -init_agen * init_g
                -init_ajuv * cue_ad_envt_high;

            // as this is generation t=0, forget about 
            // most cues when determining phenotype
            Pop[patch_i].breeders[breeder_i].phen_ad = 1.0 /
                (1.0 + exp(cue_sum));
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

            double sgnU = U < 0.0 ? -1.0 : U > 0.0 ? 1.0 : 0.0;

            // effect size of Laplace
            double x = -sdmu/sqrt(2.0) * sgnU * log(1.0 - 2 * fabs(U));

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
        ,bool const offspring_envt_high
        ,double const phen_prestige_vert
        ,double const xconformist_vert
        )
{
    // set up a bernoulli distribution that returns 0s or 1s
    // at equal probability to sample alleles from the first or
    // second set of chromosomes of diploid individuals
    bernoulli_distribution allele_sample(0.5);

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
    
    
    // inheritance of vertical social cue sensitivity values 
    offspring.asoc_vert[0] = mutation(
            mother.asoc_vert[allele_sample(rng_r)],
            mu_asoc_vert,
            sdmu_a);

    clamp(offspring.asoc_vert[0], amin, amax);

    offspring.asoc_vert[1] = mutation(
            father.asoc_vert[allele_sample(rng_r)],
            mu_asoc_vert,
            sdmu_a);

    clamp(offspring.asoc_vert[1], amin, amax);

    double asoc_vert_phen = 0.5 * (offspring.asoc_vert[0] + offspring.asoc_vert[1]);
    
    
    // inheritance of horizontal social cue sensitivity values 
    offspring.asoc_horiz[0] = mutation(
            mother.asoc_horiz[allele_sample(rng_r)],
            mu_asoc_horiz,
            sdmu_a);

    clamp(offspring.asoc_horiz[0], amin, amax);

    offspring.asoc_horiz[1] = mutation(
            father.asoc_horiz[allele_sample(rng_r)],
            mu_asoc_horiz,
            sdmu_a);

    clamp(offspring.asoc_horiz[1], amin, amax);


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

    // adult cue will be received after potential envt'al change
    //
    // has the mother observed a high cue or a low one?
    double dmat_weighting = mother.cue_ad_envt_high ? -1.0 : 1.0;

    // store the maternal phenotype for stats purposes
    offspring.phen_mat = mother.phen_ad;
    offspring.maternal_cue = mother.cue_ad_envt_high;

    // express sensitivity to maternal phenotype
    double b_phen = 0.5 * (offspring.bmat_phen[0] + offspring.bmat_phen[1]);
    assert(b_phen >= bmin);
    assert(b_phen <= bmax);

    // express sensitivity to maternal environment 
    double b_envt = 0.5 * (offspring.bmat_envt[0] + offspring.bmat_envt[1]);
    assert(b_envt >= bmin);
    assert(b_phen <= bmax);

    // generate maternal cue
    double xoff = 1.0 /
        (1.0 + exp(
                   -b_phen * (mother.phen_ad - 0.5) 
                   + 
                   dmat_weighting * b_envt));

    // generate maternal cue again 
    // but then holding the maternal environment constant
    // this is for stats purposes
    double xoff_phen_only = 1.0 /
        (1.0 + exp(
                   -b_phen * (mother.phen_ad - 0.5)));

    // generate maternal cue again 
    // but then holding the maternal phenotype constant
    // this is for stats purposes
    double xoff_envt_only = 1.0 /
        (1.0 + exp(dmat_weighting * b_envt));

    // noise in the maternal cue
    normal_distribution<> maternal_noise(0.0, sdmat);

    double mnoise = maternal_noise(rng_r);

    // calculate final value of maternal cue
    offspring.xmat = xoff + mnoise;
    offspring.xmat_phen_only = xoff_phen_only + mnoise;
    offspring.xmat_envt_only = xoff_envt_only + mnoise;
    clamp(offspring.xmat, 0.0, 1.0);

    // social learning
    offspring.xconformist_vert = xconformist_vert;
    offspring.phen_prestige_vert = phen_prestige_vert;


    // generate vertical socially learnt cue
    double xsoc_vert = 1.0 / (1.0 + exp(
                - vp_phen * (offspring.phen_prestige_vert - 0.5)
                - vc_phen * offspring.xconformist_vert));

    // also calculate vertical social cues for prestige only
    double xsoc_vert_c = 1.0 / (1.0 + exp(
                - vc_phen * offspring.xconformist_vert));

    double xsoc_vert_p = 1.0 / (1.0 + exp(
                - vp_phen * (offspring.phen_prestige_vert - 0.5)));

    normal_distribution<> social_noise(0.0, sdsoc_vert);

    double socnoise = social_noise(rng_r);

    offspring.xsoc_vert = xsoc_vert + socnoise;
    offspring.xsoc_vert_c = xsoc_vert_c + socnoise;
    offspring.xsoc_vert_p = xsoc_vert_p + socnoise;

    clamp(offspring.xsoc_vert, 0.0, 1.0);



    // expressing a juvenile phenotype
    offspring.phen_juv = 1.0 / 
        (1.0 + exp(
                   -amat_phen * offspring.xmat +
                   -agen_phen * sum_genes +
                   -ajuv_phen * offspring.cue_juv_envt_high
                   -asoc_vert_phen * offspring.xsoc_vert
                   ));
    // 
    offspring.phen_ad = NAN;
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
        assert(Pop[patch_i].n_breeders == NBreeder);
        
        // keep track of the number of breeders in a
        // high envt
        if (Pop[patch_i].envt_high)
        {
            n_high_patches += NBreeder;
        }

        // breeders endure survival selection
        for (int breeder_i = 0; breeder_i < Pop[patch_i].n_breeders; ++breeder_i)
        {
            // check whether adult phenotypes indeed exist
            assert(abs(::isnormal(Pop[patch_i].breeders[breeder_i].phen_ad)) > 0);


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

                --breeder_i;
                --Pop[patch_i].n_breeders;
            }
        } // end for (int breeder_i

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
    uniform_int_distribution<> random_local_breeder(
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

        if (phen > 0.5)
        {
            ++nhi;
        }
        else if (phen < 0.5)
        {
            ++nlo;
        }
    }

    // give conformist cue which is -1, 0 or 1
    xconformist = nhi > nlo ? 1 : nlo == nhi ? 0 : -1;
} // end void social_learning
 

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
    uniform_int_distribution<> patch_sampler(
            0,
            NPatches - 1);

    for (int patch_i = 0; patch_i < NPatches; ++patch_i)
    {
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            Individual Kid;

            // offspring is born in local patch hence sample from local parents
            if (uniform(rng_r) < 1.0 - m 
                    &&
                    Pop[patch_i].n_breeders > 0)
            {
                // set up a random number generator to 
                // sample from the remaining breeders
                uniform_int_distribution<> random_local_breeder(
                        0,
                        Pop[patch_i].n_breeders - 1);

                // prepare cues in local patch from social learning
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
        
                uniform_int_distribution<> random_remote_breeder(
                0,
                Pop[random_remote_patch].n_breeders - 1);
               
                // prepare cues in remote patch (where parents are breeding)
                // from social learning
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
            
            }

            Pop[patch_i].breeders_t1[breeder_i] = Kid;
                
            assert((int)Pop[patch_i].breeders_t1[breeder_i].g[0].size() == nloci_g);

        } // for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
    } // end for (int patch_i = 0

    bool cue_ad_envt_high;

    // auxiliary variables for horizontal social learning
    double hp, hc, asoc_horiz, xsoc_horiz;

    // random number for errors in horizontal social learning
    normal_distribution<> social_noise(0.0, sdsoc_horiz);

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
        
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            Pop[patch_i].breeders[breeder_i] = 
                Pop[patch_i].breeders_t1[breeder_i];
    
            assert((int)Pop[patch_i].breeders[breeder_i].g[0].size() == nloci_g);
            
            Pop[patch_i].n_breeders = NBreeder;
        
            // give breeder an environmental cue as adult
            Pop[patch_i].breeders[breeder_i].cue_ad_envt_high = 
                cue_ad_envt_high;
        } // all juveniles are now assigned a breeding position

        // horizontal social learning 
        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            // ad_phen should be NaN as it is yet to be set
            assert(abs(::isnan(Pop[patch_i].breeders[breeder_i].phen_ad)) > 0);

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

            asoc_horiz = 0.5 * (
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[0]
                    + 
                    Pop[patch_i].breeders[breeder_i].asoc_horiz[1]);

            Pop[patch_i].breeders[breeder_i].xconformist_horiz = xconformist;
            Pop[patch_i].breeders[breeder_i].phen_prestige_horiz = prestige_phen;

            // generate the horizontally learnt cue and add error
            xsoc_horiz = 1.0 /
                (1.0 + exp(
                           -hp * (prestige_phen - 0.5)
                           -hc * xconformist)
                ) + social_noise(rng_r);

            clamp(xsoc_horiz, 0.0, 1.0);

            // assign the horizontal social learnt 
            // conformism bias cue to the breeding individual
            Pop[patch_i].breeders[breeder_i].xsoc_horiz = 
                xsoc_horiz;

            // add social learning to the adult phenotype
            Pop[patch_i].breeders[breeder_i].phen_ad = 1.0 / 
                (1.0 + exp(
                            log(1.0 / Pop[patch_i].breeders[breeder_i].phen_juv - 1.0) 
                            -asoc_horiz * xsoc_horiz));
        }

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

    // auxiliary variable to store current generation
    int generation;

    for (generation = 0; generation < number_generations; ++generation)
    {
        // survival of adult breeders followed by reproduction
        survive();

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
