#include "individual.hpp"

Individual::Individual():
    phen_ad{0.0},
    phen_mat{0.0},
    phen_prestige{0.0},
    xmat{0.0},
    xsoc{0.0},
    xconformist{0.0},
    agen{0.0,0.0},
    amat{0.0,0.0},
    ajuv{0.0,0.0},
    asoc{0.0,0.0},
    bmat_envt{0.0,0.0},
    bmat_phen{0.0,0.0},
    dp{0.0,0.0},
    dc{0.0,0.0},
    cue_ad_envt_high{false},
    cue_juv_envt_high{false}
{
}

Individual::Individual(Individual const &other):
    phen_ad{other.phen_ad},
    phen_mat{other.phen_mat},
    phen_prestige{other.phen_prestige},
    g{other.g[0],other.g[1]},
    xmat{other.xmat},
    xsoc{other.xsoc},
    xconformist{other.xconformist},
    agen{other.agen[0],other.agen[1]},
    amat{other.amat[0],other.amat[1]},
    ajuv{other.ajuv[0],other.ajuv[1]},
    asoc{other.asoc[0],other.asoc[1]},
    bmat_envt{other.bmat_envt[0],other.bmat_envt[1]},
    bmat_phen{other.bmat_phen[0],other.bmat_phen[1]},
    dp{other.dp[0],other.dp[1]},
    dc{other.dc[0],other.dc[1]},
    cue_ad_envt_high{other.cue_ad_envt_high},
    cue_juv_envt_high{other.cue_juv_envt_high}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    phen_ad = other.phen_ad;
    phen_mat = other.phen_mat;
    phen_prestige = other.phen_prestige;
    xmat = other.xmat;
    xconformist = other.xconformist;
    xsoc = other.xsoc;

    agen[0] = other.agen[0];
    agen[1] = other.agen[1];

    amat[0] = other.amat[0];
    amat[1] = other.amat[1];

    ajuv[0] = other.ajuv[0];
    ajuv[1] = other.ajuv[1];
    
    asoc[0] = other.asoc[0];
    asoc[1] = other.asoc[1];
    
    bmat_phen[0] = other.bmat_phen[0];
    bmat_phen[1] = other.bmat_phen[1];
    
    bmat_envt[0] = other.bmat_envt[0];
    bmat_envt[1] = other.bmat_envt[1];

    dp[0] = other.dp[0];
    dp[1] = other.dp[1];
    
    dc[0] = other.dc[0];
    dc[1] = other.dc[1];

    cue_ad_envt_high = other.cue_ad_envt_high;
    cue_juv_envt_high = other.cue_juv_envt_high;

    g[0] = other.g[0];
    g[1] = other.g[1];

} // end void Individual::operator=()
