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

// random number generation
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// various functions, such as unique filename creation
#include "auxiliary.h"

#define DEBUG

// standard namespace
using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_r; // gnu scientific rng 

// parameters & variables:

// number of individuals in population
const int NPatches = 400;
const int NBreeder = 100;

// number of generations
int number_generations = 50000;

// environmental switch rate
//
// parameters below will be changed on the command line
double p = 0.0;

// whether survival selection has a sigmoidal or a quadratic function
bool survival_sigmoidal = true;

int nloci_g = 3;

// mutation rates
double mu_theta = 0.0;
double mu_phi = 0.0;
double sdmu_theta = 0.0;
double sdmu_phi = 0.0;

int seed = 0;

// write out the data every nth generation
int data_every_generation = 10;

struct Individual 
{
    // the actual phenotype of an individual
    // affecting survival
    double u;

    // loci for the genetic cue (unlinked)
    double g[nloci_g][2];

    // evolving strategy locus for the genetic cue
    double agen[2];

    // evolving strategy locus for the maternal cue
    double amat[2];
    
    // evolving strategy locus for the binary envt'al cue
    double ajuv[2];
};

struct Patch
{
    Individual breeders[NBreeder];
    bool state_high;  // high-state patch yes/no
};

Patch Pop[NPatch];


// survival function
double survival(double const u, bool const state_high)
{
    if (sigmoidal_survival)
    {
        double scalar = state_high ? 3.5 : -2.5;

        return(1.0 / (1.0 + exp(scalar + 6.0 * u)));
    }

    return(state_high ? 
            1.0 - 0.8 * (1.0 - u) * (1.0 - u)
            :
            1.0 - 0.8 * u * u
            );
}


// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{
    init_theta_a = atof(argv[1]);
    init_theta_b = atof(argv[2]);
    init_phi_a = atof(argv[3]);
    init_phi_b = atof(argv[4]);
    pmort = atof(argv[5]);
    pgood_init = atof(argv[6]);
    decay_good = atof(argv[7]);
    rgood = atof(argv[8]);
    rbad = atof(argv[9]);
    arrival_resource_decay = atof(argv[10]);
    resource_reproduce_threshold = atof(argv[11]);
    mu_theta = atof(argv[12]);
    mu_phi = atof(argv[13]);
    sdmu_theta = atof(argv[14]);
    sdmu_phi = atof(argv[15]);
    max_migration_cost = atof(argv[16]);
    min_migration_cost = atof(argv[17]);
    migration_cost_decay = atof(argv[18]);
    migration_cost_nonlinear_decay = atof(argv[19]);
    migration_cost_power = atof(argv[20]);
    tmax = atoi(argv[21]);

    // set the random seed
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng_r = gsl_rng_alloc(T);
    gsl_rng_set(rng_r, seed);
}

// write down all parameters in the file
void write_parameters(ofstream &DataFile)
{
    DataFile << endl << endl
            << "init_theta_a;" << init_theta_a << endl
            << "init_theta_b;" << init_theta_b << endl
            << "init_phi_a;" << init_phi_a << endl
            << "init_phi_b;" << init_phi_b << endl
            << "pmort;" << pmort << endl
            << "pgood_init;" << pgood_init << endl
            << "decay_good;" << decay_good << endl
            << "rgood;" << rgood << endl
            << "rbad;" << rbad << endl
            << "arrival_resource_decay;" << arrival_resource_decay << endl
            << "resource_reproduce_threshold;" << resource_reproduce_threshold << endl
            << "mu_theta;" << mu_theta << endl
            << "mu_phi;" << mu_phi << endl
            << "sdmu_theta;" << sdmu_theta << endl
            << "sdmu_phi;" << sdmu_phi << endl
            << "tmax;" << tmax << endl
            << "N;" << N << endl
            << "max_migration_cost;" << max_migration_cost << endl
            << "min_migration_cost;" << min_migration_cost << endl
            << "migration_cost_decay;" << migration_cost_decay << endl
            << "migration_cost_nonlinear_decay;" << migration_cost_nonlinear_decay << endl
            << "migration_cost_power;" << migration_cost_power << endl
            << "seed;" << seed << endl;
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
    // loop through all individuals in the wintering ground
    // and give them values
    for (int patch_i = 0; patch_i < NPatch; ++patch_i)
    {
        Pop[patch_i].high_state = gsl_rng_uniform(rng_r) < p;

        for (int breeder_i = 0; breeder_i < NBreeder; ++breeder_i)
        {
            for (int g_loc_i = 0; g_loc_i < nloci_g; ++g_loc_i)
            {
                for (int allele_i = 0; allele_i < 2; ++allele_i)
                {
                    Pop[patch_i].breeders[breeder_i].g[g_loc_i][allele_i] = init_g;
                }
            }
        }
    }

    NWinter = N;
}

// individuals in both summer and winter populations
// die at a certain rate. If this function is called in winter
// the summer pool will be empty so no individuals die there
// If this function is called in summer, however, both the summer
// ground individuals die, as well as the individuals who have
// stayed at the wintering ground
void mortality()
{
    for (int i = 0; i < NWinter;++i)
    {
        // individual dies; replace with end of the stack individual
        if (gsl_rng_uniform(rng_r) < pmort)
        {
            WinterPop[i] = WinterPop[NWinter - 1];
            --NWinter;
            --i;
        }
    }

    for (int i = 0; i < NSummer;++i)
    {
        // individual dies; replace with end of the stack individual
        if (gsl_rng_uniform(rng_r) < pmort)
        {
            SummerPop[i] = SummerPop[NSummer - 1];
            --NSummer;
            --i;
        }
    }
}

// remove individuals from the staging pool and put them
// back in the original population
void clear_staging_pool()
{
    // put individuals from staging pool (which haven't migrated) 
    // back in the original population
    for (int i = 0; i < NStaging; ++i)
    {
        WinterPop[NWinter++] = StagingPool[i];

    }

    // just double check that NWinter does not exceed max population size
    assert(NWinter <= N);

    NStaging = 0;
}

// the dynamics of the population at the wintering ground
void winter_dynamics(int t)
{
    // individuals forage
    // individuals accumulate resources
    // individuals make dispersal decisions

    // probability of encountering a good resource
    double pgood = pgood_init - decay_good * t;

    // set lower boundary to the probability
    if (pgood <= 0)
    {
        pgood = 0;
    }

    // foraging of individuals who are just at the wintering site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < NWinter; ++i)
    {
        if (gsl_rng_uniform(rng_r) < pgood) // good resource chosen
        {
            WinterPop[i].resources += rgood;
        }
        else
        {
            WinterPop[i].resources += rbad;
        }
    
    } // ok resource dynamic done


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < NStaging; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
        if (gsl_rng_uniform(rng_r) < pgood) // good resource chosen
        {
            StagingPool[i].resources += rgood;
        }
        else
        {
            StagingPool[i].resources += rbad;
        }
    }

    assert(NWinter <= N);
    assert(NWinter >= 0);

    double psignal = 0.0;

    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < NWinter; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1])
            + 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * WinterPop[i].resources;

        // does individual want to signal to others to be ready for departure?
        if (gsl_rng_uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[NStaging] = WinterPop[i];
            ++NStaging; // increment the number of individuals in the staging pool

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // delete this individual from the winter population
            WinterPop[i] = WinterPop[NWinter - 1];

            // decrement the number of individuals in the winter population
            --NWinter;
            --i;
        }
    }

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int NSummer_old = NSummer;

    int NFlock = 0;

    double pdisperse = 0.0;

    int NStaging_start = NStaging;

    // actual dispersal
    for (int i = 0; i < NStaging; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * NStaging_start;

        // yes individual goes
        if (gsl_rng_uniform(rng_r) < pdisperse)
        {
            SummerPop[NSummer] = StagingPool[i];
            ++NSummer;
            
            assert(NSummer <= N);


            // delete this individual from the staging population
            StagingPool[i] = StagingPool[NStaging - 1];

            // decrement the number of individuals in the staging population
            --NStaging;
            --i;

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        }
    }

    double total_migration_cost;

    // update resource levels for all new individuals that have just
    // been added to the summer pool dependent on their flock size
    for (int i = NSummer_old; i < NSummer; ++i)
    {
        total_migration_cost = max_migration_cost - migration_cost_decay * NFlock - migration_cost_nonlinear_decay * pow(NFlock,migration_cost_power);

        if (total_migration_cost < min_migration_cost)
        {
            total_migration_cost = min_migration_cost;
        }

        assert(total_migration_cost >= 0.0);
        assert(total_migration_cost <= 1.0);

        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources * total_migration_cost;

        // and reduce it by time of arrival
        // TODO think more about this function
        SummerPop[i].resources -= arrival_resource_decay * t;

        if (SummerPop[i].resources < 0)
        {
            SummerPop[i].resources = 0;
        }
    }

    // add current dispersal flock size to stats
    mean_flock_size_winter += NFlock;
    mean_staging_size_winter += NStaging_start;
} // end winter_dynamics

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

// create a new offspring
void create_offspring(Individual &mother, Individual &father, Individual &offspring)
{
    offspring.resources = 0.0;

    // inherit theta loci

    // each parental allele has probability 0.5 to make it into offspring
    offspring.theta_a[0] = mutation(mother.theta_a[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);

    offspring.theta_a[1] = mutation(father.theta_a[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);

    

    
    offspring.theta_b[0] = mutation(mother.theta_b[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);
    offspring.theta_b[1] = mutation(father.theta_b[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);

    // inherit phi loci
    offspring.phi_a[0] = mutation(mother.phi_a[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);
    offspring.phi_a[1] = mutation(father.phi_a[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);
    
    offspring.phi_b[0] = mutation(mother.phi_b[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);
    offspring.phi_b[1] = mutation(father.phi_b[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);

    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        // put boundaries on the elevation between 0 and 1
        // to help with the interpretation of the evolved values of the slope
        if (offspring.theta_a[allele_i] < 0.0)
        {
            offspring.theta_a[allele_i] = 0.0;
        }
        else if (offspring.theta_a[allele_i] > 1.0)
        {
            offspring.theta_a[allele_i] = 1.0;
        }
        
        if (offspring.phi_a[allele_i] < 0.0)
        {
            offspring.phi_a[allele_i] = 0.0;
        }
        else if (offspring.phi_a[allele_i] > 1.0)
        {
            offspring.phi_a[allele_i] = 1.0;
        }
    }   
}



// in summery they reproduce dependent on 
// resources and arrival time
void summer_reproduction(ofstream &DataFile)
{
    // auxiliary variables storing current mom and dad
    Individual mother, father;

    // auxiliary variable specifing the rounded amount of a mother's
    // resources
    int resource_integer;

    /// auxilary variable specifying the id of the randomly sampled
    // fater
    int father_id;

    // see if population is extinct
    if (NSummer == 1)
    {
        // quit if extinct 
        write_parameters(DataFile);
        exit(1);
    }

    // use a flexible array for the kids
    vector<Individual> Kids;

    // mating dynamic. Presumes that there an even 
    // number of individuals
    // so we just discard the last individual
    for (int i = 0; i < NSummer; ++i)
    {
        // get the mother
        mother = SummerPop[i];

        // if mom does not meet minimum standards
        // no reproduction through female function
        if (mother.resources < resource_reproduce_threshold)
        {
            continue;
        }

        // now randomly select a father
        do {
            // sample integer uniformly between 0 and NSummer
            // (not including NSummer itself)
            father_id = gsl_rng_uniform_int(rng_r, NSummer);
        }
        while (father_id == i);

        father = SummerPop[father_id];

        // translate maternal resources to numbers of offspring
        //
        // first round to lowest integer
        resource_integer = floor(mother.resources);

        // TODO (slightly digressing): can we come up with an analytical 
        // description of this rounding process of w into integer values?
        if (gsl_rng_uniform(rng_r) < mother.resources - resource_integer)
        {
            // make an additional offspring
            ++resource_integer;
        }
        
        // for each parent create the offspring
        for (int kid_i = 0; kid_i < resource_integer; ++kid_i)
        {
            Individual kid;

            create_offspring(mother, father, kid);

            // add kid to the stack
            Kids.push_back(kid);
        }
    }

    NKids = Kids.size();

    // number of dead individuals is the max population
    // minus the current individuals in the summer population
    // minus the current individuals who stayed at the wintering ground
    int Ndead = N - NSummer - NWinter;

    int random_kid = 0;

    // recruit new individuals to the summer pool
    for (int i = 0; i < Ndead; ++i)
    {
        // no kids left to recruit
        if (Kids.size() == 0)
        {
            break;
        }

        random_kid = gsl_rng_uniform_int(rng_r, Kids.size());

        // add random kid to population
        SummerPop[NSummer++] = Kids[random_kid];

        //  delete random kid as it has been sampled
        Kids[random_kid] = Kids[Kids.size()];

        Kids.pop_back();
    }
}

// gaining resources at breeding ground
// & fly back
void summer_dynamics(int t)
{

    // probability of encountering a good resource
    double pgood = pgood_init - decay_good * t;

    // foraging of individuals who are just at the wintering site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < NSummer; ++i)
    {
        if (gsl_rng_uniform(rng_r) < pgood) // good resource chosen
        {
            SummerPop[i].resources += rgood;
        }
        else
        {
            SummerPop[i].resources += rbad;
        }
    
    } // ok resource dynamic done


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < NStaging; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
        if (gsl_rng_uniform(rng_r) < pgood) // good resource chosen
        {
            StagingPool[i].resources += rgood;
        }
        else
        {
            StagingPool[i].resources += rbad;
        }
    }

    assert(NSummer<= N);
    assert(NSummer>= 0);

    double psignal = 0.0;
    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < NSummer; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1])
            + 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]) * SummerPop[i].resources;

        // does individual want to signal to others to be ready for departure?
        if (gsl_rng_uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[NStaging] = SummerPop[i];
            ++NStaging; // increment the number of individuals in the staging pool

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // delete this individual from the winter population
            SummerPop[i] = SummerPop[NSummer - 1];

            // decrement the number of individuals in the winter population
            --NSummer;
            --i;
        }
    }

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int NWinter_old = NWinter;

    int NFlock = 0;

    double pdisperse = 0.0;

    int NStaging_start = NStaging;

    // actual dispersal
    for (int i = 0; i < NStaging; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * NStaging_start;

        // yes individual goes
        if (gsl_rng_uniform(rng_r) < pdisperse)
        {
            WinterPop[NWinter] = StagingPool[i];
            ++NWinter;
            
            assert(NWinter <= N);


            // delete this individual from the staging population
            StagingPool[i] = StagingPool[NStaging - 1];

            // decrement the number of individuals in the staging population
            --NStaging;
            --i;

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        }
    }

    // update resource levels for all new individuals that have just
    // been added to the pool dependent on their flock size
    for (int i = NWinter_old; i < NWinter; ++i)
    {
        // TODO: think about relationship between flock size and
        // resource reduction
        WinterPop[i].resources -= 1.0 / NFlock;

        // and reduce it by time of arrival
        // TODO think more about this function
        WinterPop[i].resources -= arrival_resource_decay * t;
    }
    
    // add current dispersal flock size to stats
    mean_flock_size_summer += NFlock;
    mean_staging_size_summer += NStaging_start;
}


// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    string filename = "sim_migration";
    create_filename(filename);
    ofstream DataFile(filename.c_str());  // output file 

    init_arguments(argc, argv);

    write_data_headers(DataFile);

    init_population();

    for (int generation = 0; generation < number_generations; ++generation)
    {
        mean_flock_size_winter = 0.0;
        mean_staging_size_winter = 0.0;

        // time during winter (i.e., days)
        // during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            winter_dynamics(t);
        }

        // now take averages over all timesteps that individuals can join groups
        mean_flock_size_winter /= tmax;
        mean_staging_size_winter /= tmax;
        
        // all individuals that wanted to migrate have migrated now
        // all remainers are going to stay at wintering ground
        clear_staging_pool();

        // let individuals die with a certain probability 
        mortality();


        // have individuals reproduce after they migrated to the summer spot
        summer_reproduction(DataFile);
        
        // set flock size stats to 0 before summer dynamics starts
        mean_flock_size_summer = 0.0;
        mean_staging_size_summer = 0.0;

        // time during summer (i.e., days)
        // during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            summer_dynamics(t);
        }
        
        // now take averages over all timesteps that individuals can join groups
        mean_flock_size_summer /= tmax;
        mean_staging_size_summer /= tmax;


        // all individuals who remain at the summer ground die
        NSummer = 0;
        NStaging = 0;

        // let individuals die with a certain probability 
        mortality();

        if (generation % skip == 0)
        {
            write_stats(DataFile, generation, 2);
        }
    }

    write_parameters(DataFile);
}
