#include "individual.hpp"

Individual::Individual():
    phen_juv{0.0},
    phen_ad{0.0},
    phen_mat{0.0},
    phen_mat_error{0.0},
    phen_prestige_horiz{0.0},
    phen_prestige_horiz_error{0.0},
    phen_prestige_vert{0.0},
    phen_prestige_vert_error{0.0},
    xconformist_horiz{0.0},
    xconformist_horiz_error{0.0},
    xconformist_vert{0.0},
    xconformist_vert_error{0.0},
    aintercept{0.0,0.0},
    agen{0.0,0.0},
    ajuv{0.0,0.0},
    bmat_envt{0.0,0.0},
    bmat_phen{0.0,0.0},
    vp{0.0,0.0},
    vc{0.0,0.0},
    hp{0.0,0.0},
    hc{0.0,0.0},
    cue_juv_envt_high{false},
    cue_ad_envt_high{false},
    maternal_envt_cue{false},
    maternal_envt_cue_error{false}
{
}

Individual::Individual(Individual const &other):
    phen_juv{other.phen_juv},
    phen_ad{other.phen_ad},
    phen_mat{other.phen_mat},
    phen_mat_error{other.phen_mat_error},
    phen_prestige_horiz{other.phen_prestige_horiz},
    phen_prestige_horiz_error{other.phen_prestige_horiz_error},
    phen_prestige_vert{other.phen_prestige_vert},
    phen_prestige_vert_error{other.phen_prestige_vert_error},
    g{other.g[0],other.g[1]},
    xconformist_horiz{other.xconformist_horiz},
    xconformist_horiz_error{other.xconformist_horiz_error},
    xconformist_vert{other.xconformist_vert},
    xconformist_vert_error{other.xconformist_vert_error},
    aintercept{other.aintercept[0],other.aintercept[1]},
    agen{other.agen[0],other.agen[1]},
    ajuv{other.ajuv[0],other.ajuv[1]},
    bmat_envt{other.bmat_envt[0],other.bmat_envt[1]},
    bmat_phen{other.bmat_phen[0],other.bmat_phen[1]},
    vp{other.vp[0],other.vp[1]},
    vc{other.vc[0],other.vc[1]},
    hp{other.hp[0],other.hp[1]},
    hc{other.hc[0],other.hc[1]},
    cue_juv_envt_high{other.cue_juv_envt_high},
    cue_ad_envt_high{other.cue_ad_envt_high},
    maternal_envt_cue{other.maternal_envt_cue},
    maternal_envt_cue_error{other.maternal_envt_cue_error}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    phen_juv = other.phen_juv;
    phen_ad = other.phen_ad;
    phen_mat = other.phen_mat;
    phen_mat_error = other.phen_mat_error;

    phen_prestige_horiz = other.phen_prestige_horiz;
    phen_prestige_horiz_error = other.phen_prestige_horiz_error;
    phen_prestige_vert = other.phen_prestige_vert;
    phen_prestige_vert_error = other.phen_prestige_vert_error;

    xconformist_horiz = other.xconformist_horiz;
    xconformist_horiz_error = other.xconformist_horiz_error;
    xconformist_vert = other.xconformist_vert;
    xconformist_vert_error = other.xconformist_vert_error;

    agen[0] = other.agen[0];
    agen[1] = other.agen[1];
    
    ajuv[0] = other.ajuv[0];
    ajuv[1] = other.ajuv[1];
    
    aintercept[0] = other.aintercept[0];
    aintercept[1] = other.aintercept[1];
    
    bmat_phen[0] = other.bmat_phen[0];
    bmat_phen[1] = other.bmat_phen[1];
    
    bmat_envt[0] = other.bmat_envt[0];
    bmat_envt[1] = other.bmat_envt[1];

    hp[0] = other.hp[0];
    hp[1] = other.hp[1];
    
    hc[0] = other.hc[0];
    hc[1] = other.hc[1];
    
    vp[0] = other.vp[0];
    vp[1] = other.vp[1];
    
    vc[0] = other.vc[0];
    vc[1] = other.vc[1];

    cue_juv_envt_high = other.cue_juv_envt_high;
    cue_ad_envt_high = other.cue_ad_envt_high;

    g[0] = other.g[0];
    g[1] = other.g[1];
    
    maternal_envt_cue = other.maternal_envt_cue;
    maternal_envt_cue_error = other.maternal_envt_cue_error;
} // end void Individual::operator=()
