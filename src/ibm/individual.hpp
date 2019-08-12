#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <vector>

class Individual
{
    public:
        
        // the actual phenotype of an individual
        // affecting survival
        double ad_phen;

        // diploid loci for the genetic cue (unlinked)
        std::vector < double > g[2];

        // evolving strategy locus for the genetic cue
        double agen[2];

        // evolving strategy locus for the maternal cue
        double amat[2];
        
        // evolving strategy locus for the binary envt'al cue
        double ajuv[2];

        // evolving strategy to weigh maternal phenotypic cues
        double bmat_phen[2];

        // evolving strategy to weigh maternal enviromental cues
        double bmat_envt[2];

        // juvenile cue
        bool cue_ad_envt_high;

        // adult cue
        bool cue_juv_envt_high;

        Individual();

        Individual(Individual const &other);

        void operator=(Individual const &other);
};

#endif
