#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <vector>

class Individual
{
    public:
        
        // the actual phenotype of an individual
        // affecting survival
        double phen_ad;
        
        // the maternal phenotype (stats purposes)
        double phen_mat;

        // phenotypic cue from prestige bias
        double phen_prestige;

        // diploid loci for the genetic cue (unlinked)
        std::vector < double > g[2];

        // maternal cue provided to the offspring (stats purposes)
        // including noise
        double xmat;
        
        // socially learnt cue
        double xsoc;
        
        // socially learnt conformism-based cue 
        double xconformist;

        // evolving strategy locus for the genetic cue
        double agen[2];

        // evolving strategy locus for the maternal cue
        double amat[2];
        
        // evolving strategy locus for the binary envt'al cue
        double ajuv[2];

        // evolving strategy locus for the social learning strategy
        double asoc[2];

        // evolving strategy to weigh maternal enviromental cues
        double bmat_envt[2];

        // evolving strategy to weigh maternal phenotypic cues
        double bmat_phen[2];

        // evolving strategy to weigh cues from prestige-based social learning
        double dp[2];
        
        // evolving strategy to weigh cues from conformism-based social learning
        double dc[2];

        // adult cue
        bool cue_ad_envt_high;

        // juvenile cue
        bool cue_juv_envt_high;

        Individual();

        Individual(Individual const &other);

        void operator=(Individual const &other);
};

#endif
