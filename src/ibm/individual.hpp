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
        
        // the actual phenotype of an individual
        // affecting survival
        double phen_ad;
        
        // the maternal phenotype (stats purposes)
        double phen_mat;

        
        // phenotypic cue from horizontally learnt prestige bias
        double phen_prestige_horiz;

        // phenotypic cue from vertically learnt prestige bias
        double phen_prestige_vert;

        // diploid loci for the genetic cue (unlinked)
        std::vector < double > g[2];

        // maternal cue provided to the offspring (stats purposes)
        // including noise
        double xmat;
        
        // horizontal socially learnt cue
        double xsoc_horiz;
        
        // vertical socially learnt cue
        double xsoc_vert;
        
        // horizontal socially learnt conformism-based cue 
        double xconformist_horiz;
        
        // vertical socially learnt conformism-based cue 
        double xconformist_vert;

        // evolving strategy locus for the genetic cue
        double agen[2];

        // evolving strategy locus for the maternal cue
        double amat[2];
        
        // evolving strategy locus for the binary envt'al cue
        double ajuv[2];

        // evolving strategy locus for the vertical social learning strategy
        double asoc_vert[2];
        
        // evolving strategy locus for the vertical social learning strategy
        double asoc_horiz[2];

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
        
        // the maternal cue (stats purposes)
        bool maternal_cue;

        // adult cue
        bool cue_ad_envt_high;

        // juvenile cue
        bool cue_juv_envt_high;

        // default constructor
        Individual();

        // copy constructor
        Individual(Individual const &other);

        // assignment operator
        void operator=(Individual const &other);
};

#endif
