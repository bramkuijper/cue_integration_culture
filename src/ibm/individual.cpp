#include "individual.hpp"

Individual::Individual():
    ad_phen{0.0},
    ad_mat{0.0},
    xmat{0.0},
    agen{0.0,0.0},
    amat{0.0,0.0},
    ajuv{0.0,0.0},
    bmat_envt{0.0,0.0},
    bmat_phen{0.0,0.0},
    cue_ad_envt_high{false},
    cue_juv_envt_high{false}
{
}

Individual::Individual(Individual const &other):
    ad_phen{other.ad_phen},
    ad_mat{other.ad_mat},
    g{other.g[0],other.g[1]},
    xmat{other.xmat},
    agen{other.agen[0],other.agen[1]},
    amat{other.amat[0],other.amat[1]},
    ajuv{other.ajuv[0],other.ajuv[1]},
    bmat_envt{other.bmat_envt[0],other.bmat_envt[1]},
    bmat_phen{other.bmat_phen[0],other.bmat_phen[1]},
    cue_ad_envt_high{other.cue_ad_envt_high},
    cue_juv_envt_high{other.cue_juv_envt_high}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    ad_phen = other.ad_phen;
    ad_mat = other.ad_mat;
    xmat = other.xmat;

    agen[0] = other.agen[0];
    agen[1] = other.agen[1];

    amat[0] = other.amat[0];
    amat[1] = other.amat[1];

    ajuv[0] = other.ajuv[0];
    ajuv[1] = other.ajuv[1];
    
    bmat_phen[0] = other.bmat_phen[0];
    bmat_phen[1] = other.bmat_phen[1];
    
    bmat_envt[0] = other.bmat_envt[0];
    bmat_envt[1] = other.bmat_envt[1];

    cue_ad_envt_high = other.cue_ad_envt_high;
    cue_juv_envt_high = other.cue_juv_envt_high;

    g[0] = other.g[0];
    g[1] = other.g[1];

} // end void Individual::operator=()
