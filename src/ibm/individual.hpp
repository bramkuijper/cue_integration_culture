#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

#include <vector>

class Individual
{
    public:
        
        // an individual's juvenile phenotype
        // which is the adult phenotype prior to 
        // horizontal social learning
        double phen_juv;
        double phen_juv_logistic;
        
        // the actual phenotype of an individual
        // affecting survival
        double phen_ad;
        double phen_ad_logistic;
        
        // the maternal phenotype (stats purposes)
        double phen_mat;
        
        // the maternal phenotype (stats purposes)
        double phen_mat_error;
        
        // phenotypic cue from horizontally learnt prestige bias
        double phen_prestige_horiz;
        double phen_prestige_horiz_error; // same but with error

        // phenotypic cue from vertically learnt prestige bias
        double phen_prestige_vert;
        double phen_prestige_vert_error; // same but with error

        // diploid loci for the genetic cue (unlinked)
        std::vector < double > g[2];

        // horizontal socially learnt conformism-based cue 
        double xconformist_horiz;
        
        // horizontal socially learnt conformism-based cue including noise
        double xconformist_horiz_error;

        // vertical socially learnt conformism-based cue 
        double xconformist_vert;
        
        // vertical socially learnt conformism-based cue including noise
        double xconformist_vert_error;

        // basic intercept locus against which 
        // other cues are interpreted
        double aintercept[2];

        // evolving strategy locus for the genetic cue
        double agen[2];

        // evolving strategy locus for the binary envt'al cue
        double ajuv[2];

        // evolving strategy to weigh maternal enviromental cues
        double bmat_envt[2];

        // evolving strategy to weigh maternal phenotypic cues
        double bmat_phen[2];

        // evolving strategy to weigh cues from 
        // vertical
        // prestige-based social learning
        double vp[2];
        
        // evolving strategy to weigh cues from 
        // vertical
        // conformism-based social learning
        double vc[2];
        
        
        // evolving strategy to weigh cues from 
        // horizontal
        // prestige-based social learning
        double hp[2];
        
        // evolving strategy to weigh cues from 
        // horizontal
        // conformism-based social learning
        double hc[2];
        
        // juvenile cue
        bool cue_juv_envt_high;
        
        // adult cue
        bool cue_ad_envt_high;

        bool envt_high_selection;
        bool envt_high_previous;

        // the maternal environmental cue (stats purposes)
        bool maternal_envt_cue;
        bool maternal_envt_cue_error;


        // default constructor
        Individual();

        // copy constructor
        Individual(Individual const &other);

        // assignment operator
        void operator=(Individual const &other);
};

#endif
