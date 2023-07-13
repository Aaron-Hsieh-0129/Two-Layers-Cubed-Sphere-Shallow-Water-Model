#include "iteration.hpp"

#if defined(SecondOrderSpace)
void Iteration::ph_pt_2(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_h = model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j];
                dy_for_h = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1];

                psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                              ((model.sqrtG[i+1][j] * model.csswm[p].h[i+1][j] * (model.gUpper[i+1][j][0] * model.csswm[p].u[i+1][j] + model.gUpper[i+1][j][1] * model.csswm[p].v[i+1][j])) - 
                               (model.sqrtG[i-1][j] * model.csswm[p].h[i-1][j] * (model.gUpper[i-1][j][0] * model.csswm[p].u[i-1][j] + model.gUpper[i-1][j][1] * model.csswm[p].v[i-1][j])));

                psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                              ((model.sqrtG[i][j+1] * model.csswm[p].h[i][j+1] * (model.gUpper[i][j+1][2] * model.csswm[p].u[i][j+1] + model.gUpper[i][j+1][3] * model.csswm[p].v[i][j+1])) - 
                               (model.sqrtG[i][j-1] * model.csswm[p].h[i][j-1] * (model.gUpper[i][j-1][2] * model.csswm[p].u[i][j-1] + model.gUpper[i][j-1][3] * model.csswm[p].v[i][j-1])));
            
                model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                
                #ifdef DIFFUSION
                    model.csswm[p].hp[i][j] += D2T * KX * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx_for_h, 2) + 
                                               D2T * KY * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy_for_h, 2);
                #endif
            }
        }
    }
    return;
}

void Iteration::pu_pt_2(CSSWM &model) {
    double dx_for_u = 0, dy_for_u = 0;
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;
    double f;
    #ifdef Mountain
        double pgHs_px = 0.;
    #endif

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_u = model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j];
                dy_for_u = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1];

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #else
                    f = 0;
                #endif
                
                pgH_px = GRAVITY / dx_for_u * (model.csswm[p].h[i+1][j] - model.csswm[p].h[i-1][j]);

                pU2_px = 0.5 / dx_for_u * 
                        (model.gUpper[i+1][j][0] * pow(model.csswm[p].u[i+1][j], 2) - 
                         model.gUpper[i-1][j][0] * pow(model.csswm[p].u[i-1][j], 2));

                pUV_px = 1 / dx_for_u * 
                        (model.gUpper[i+1][j][1] * model.csswm[p].u[i+1][j] * model.csswm[p].v[i+1][j] - 
                         model.gUpper[i-1][j][1] * model.csswm[p].u[i-1][j] * model.csswm[p].v[i-1][j]);

                pV2_px = 0.5 / dx_for_u * 
                        (model.gUpper[i+1][j][3] * pow(model.csswm[p].v[i+1][j], 2) -
                         model.gUpper[i-1][j][3] * pow(model.csswm[p].v[i-1][j], 2));

                rotationU = (((model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_u) - 
                             ((model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_u) + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][2] * model.csswm[p].u[i][j] + model.gUpper[i][j][3] * model.csswm[p].v[i][j]);
            

                #ifdef Mountain
                    pgHs_px = GRAVITY / dx_for_u * (model.csswm[p].hs[i+1][j] - model.csswm[p].hs[i-1][j]);
                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pgHs_px - pU2_px - pUV_px - pV2_px + rotationU);
                #else
                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);
                #endif

                #ifdef DIFFUSION
                    model.csswm[p].up[i][j] += D2T * KX * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx_for_u, 2) + 
                                               D2T * KY * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy_for_u, 2);
                #endif
            }
        }
    }
    return;
}

void Iteration::pv_pt_2(CSSWM &model) {
    double dx_for_v = 0, dy_for_v = 0;
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;
    double f;
    #ifdef Mountain
        double pgHs_py = 0.;
    #endif

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                dx_for_v = model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j];
                dy_for_v = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1];

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #else
                    f = 0;
                #endif

                pgH_py = GRAVITY / dy_for_v * (model.csswm[p].h[i][j+1] - model.csswm[p].h[i][j-1]);

                pU2_py = 0.5 / dy_for_v * 
                        (model.gUpper[i][j+1][0] * pow(model.csswm[p].u[i][j+1], 2) - 
                         model.gUpper[i][j-1][0] * pow(model.csswm[p].u[i][j-1], 2));

                pUV_py = 1 / dy_for_v * 
                        (model.gUpper[i][j+1][1] * model.csswm[p].u[i][j+1] * model.csswm[p].v[i][j+1] - 
                         model.gUpper[i][j-1][1] * model.csswm[p].u[i][j-1] * model.csswm[p].v[i][j-1]);

                pV2_py = 0.5 / dy_for_v * 
                        (model.gUpper[i][j+1][3] * pow(model.csswm[p].v[i][j+1], 2) -
                         model.gUpper[i][j-1][3] * pow(model.csswm[p].v[i][j-1], 2));

                rotationV = (((model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_v) - 
                             ((model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_v) + model.sqrtG[i][j] * f) * 
                            (model.gUpper[i][j][0] * model.csswm[p].u[i][j] + model.gUpper[i][j][1] * model.csswm[p].v[i][j]);

                
                #ifdef Mountain
                    pgHs_py = GRAVITY / dy_for_v * (model.csswm[p].hs[i][j+1] - model.csswm[p].hs[i][j-1]);;
                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pgHs_py - pU2_py - pUV_py - pV2_py - rotationV);
                #else
                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                #endif

                #ifdef DIFFUSION
                    model.csswm[p].vp[i][j] += D2T * KX * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx_for_v, 2) + 
                                               D2T * KY * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy_for_v, 2);
                #endif
            }
        }
    }
    return;
}
#elif defined(FourthOrderSpace)
void Iteration::ph_pt_4(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx_for_h = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy_for_h = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);
                
                for (int k = 0; k < NZ; k++) {
                    psqrtGHU_px = (1. / (model.sqrtG[i][j] * 12. * dx_for_h)) * 
                                  (-1.*(model.sqrtG[i+2][j] * model.csswm[p].h[i+2][j][k] * (model.gUpper[i+2][j][0] * model.csswm[p].u[i+2][j][k] + model.gUpper[i+2][j][1] * model.csswm[p].v[i+2][j][k]))
                                   +8.*(model.sqrtG[i+1][j] * model.csswm[p].h[i+1][j][k] * (model.gUpper[i+1][j][0] * model.csswm[p].u[i+1][j][k] + model.gUpper[i+1][j][1] * model.csswm[p].v[i+1][j][k]))
                                   -8.*(model.sqrtG[i-1][j] * model.csswm[p].h[i-1][j][k] * (model.gUpper[i-1][j][0] * model.csswm[p].u[i-1][j][k] + model.gUpper[i-1][j][1] * model.csswm[p].v[i-1][j][k]))
                                   +1.*(model.sqrtG[i-2][j] * model.csswm[p].h[i-2][j][k] * (model.gUpper[i-2][j][0] * model.csswm[p].u[i-2][j][k] + model.gUpper[i-2][j][1] * model.csswm[p].v[i-2][j][k])));

                    psqrtGHU_py = (1. / (model.sqrtG[i][j] * 12. * dy_for_h)) * 
                                  (-1.*(model.sqrtG[i][j+2] * model.csswm[p].h[i][j+2][k] * (model.gUpper[i][j+2][2] * model.csswm[p].u[i][j+2][k] + model.gUpper[i][j+2][3] * model.csswm[p].v[i][j+2][k])) 
                                   +8.*(model.sqrtG[i][j+1] * model.csswm[p].h[i][j+1][k] * (model.gUpper[i][j+1][2] * model.csswm[p].u[i][j+1][k] + model.gUpper[i][j+1][3] * model.csswm[p].v[i][j+1][k]))
                                   -8.*(model.sqrtG[i][j-1] * model.csswm[p].h[i][j-1][k] * (model.gUpper[i][j-1][2] * model.csswm[p].u[i][j-1][k] + model.gUpper[i][j-1][3] * model.csswm[p].v[i][j-1][k]))
                                   +1.*(model.sqrtG[i][j-2] * model.csswm[p].h[i][j-2][k] * (model.gUpper[i][j-2][2] * model.csswm[p].u[i][j-2][k] + model.gUpper[i][j-2][3] * model.csswm[p].v[i][j-2][k])));
                
                    model.csswm[p].hp[i][j][k] = model.csswm[p].hm[i][j][k] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                    
                    #if defined(DIFFUSION)
                        model.csswm[p].hp[i][j][k] += D2T * KX * (model.csswm[p].hm[i+1][j][k] - 2. * model.csswm[p].hm[i][j][k] + model.csswm[p].hm[i-1][j][k]) / pow(dx_for_h, 2) + 
                                                      D2T * KY * (model.csswm[p].hm[i][j+1][k] - 2. * model.csswm[p].hm[i][j][k] + model.csswm[p].hm[i][j-1][k]) / pow(dy_for_h, 2);
                    #elif defined(DIFFUSION4)
                        model.csswm[p].hp[i][j][k] += D2T * KX * (-model.csswm[p].hm[i+2][j][k] + 16.*model.csswm[p].hm[i+1][j][k] - 30.*model.csswm[p].hm[i][j][k] + 16.*model.csswm[p].hm[i-1][j][k] - model.csswm[p].hm[i-2][j][k]) / (12.*pow(dx_for_h, 2)) + 
                                                      D2T * KY * (-model.csswm[p].hm[i][j+2][k] + 16.*model.csswm[p].hm[i][j+1][k] - 30.*model.csswm[p].hm[i][j][k] + 16.*model.csswm[p].hm[i][j-1][k] - model.csswm[p].hm[i][j-2][k]) / (12.*pow(dy_for_h, 2));
                    #elif defined(DIFFUSIONLARS)
                        model.csswm[p].hp[i][j][k] += D2T * KX * (model.csswm[p].h[i+1][j][k] - model.csswm[p].hp[i][j][k] - model.csswm[p].hm[i][j][k] + model.csswm[p].h[i-1][j][k]) / pow(dx_for_h, 2) + 
                                                      D2T * KY * (model.csswm[p].h[i][j+1][k] - model.csswm[p].hp[i][j][k] - model.csswm[p].hm[i][j][k] + model.csswm[p].h[i][j-1][k]) / pow(dy_for_h, 2);
                    #endif
                }
            }
        }
    }
    return;
}

void Iteration::pu_pt_4(CSSWM &model) {
    double dx_for_u = 0, dy_for_u = 0;
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;
    double f;
    #ifdef Mountain
        double pgHs_px = 0.;
    #endif
    double pgH_px_all = 0;

    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx_for_u = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy_for_u = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #else
                    f = 0;
                #endif
                pgH_px_all = 0;

                for (int k = 0; k < NZ; k++) {
                    pgH_px_all += GRAVITY / (12.*dx_for_u) * (-1*model.csswm[p].h[i+2][j][k] + 8*model.csswm[p].h[i+1][j][k] - 8*model.csswm[p].h[i-1][j][k] + 1*model.csswm[p].h[i-2][j][k]);
                }
                #if defined (Mountain)
                    pgH_px_all += GRAVITY / (12.*dx_for_u) * (-1.*model.csswm[p].hs[i+2][j] + 8.*model.csswm[p].hs[i+1][j] - 8.*model.csswm[p].hs[i-1][j] + 1.*model.csswm[p].hs[i-2][j]);
                #endif

                for (int k = 0; k < NZ; k++) {
                    if (k == 0) {
                        #if defined (Mountain)
                            pgH_px = RHO2 / RHO1 * 
                                     (pgH_px_all + (GRAVITY * (RHO1-RHO2)/RHO2) / (12.*dx_for_u) * (-1.*model.csswm[p].h[i+2][j][1] + 8.*model.csswm[p].h[i+1][j][1] - 8.*model.csswm[p].h[i-1][j][1] + 1.*model.csswm[p].h[i-2][j][1] + 
                                                                                                    -1.*model.csswm[p].hs[i+2][j] + 8.*model.csswm[p].hs[i+1][j] - 8.*model.csswm[p].hs[i-1][j] + 1.*model.csswm[p].hs[i-2][j]));
                        #else
                            pgH_px = RHO2 / RHO1 * 
                                     (pgH_px_all + (GRAVITY * (RHO1-RHO2)/RHO2) / (12.*dx_for_u) * (-1.*model.csswm[p].h[i+2][j][1] + 8.*model.csswm[p].h[i+1][j][1] - 8.*model.csswm[p].h[i-1][j][1] + 1.*model.csswm[p].h[i-2][j][1]));
                        #endif
                        
                    }
                    else if (k == 1) {
                        pgH_px = pgH_px_all;
                    }

                    pU2_px = 0.5 / (12.*dx_for_u) * 
                            (-1.*(model.gUpper[i+2][j][0] * pow(model.csswm[p].u[i+2][j][k], 2))
                             +8.*(model.gUpper[i+1][j][0] * pow(model.csswm[p].u[i+1][j][k], 2))
                             -8.*(model.gUpper[i-1][j][0] * pow(model.csswm[p].u[i-1][j][k], 2))
                             +1.*(model.gUpper[i-2][j][0] * pow(model.csswm[p].u[i-2][j][k], 2)));

                    pUV_px = 1./ (12.*dx_for_u) * 
                            (-1.*(model.gUpper[i+2][j][1] * model.csswm[p].u[i+2][j][k] * model.csswm[p].v[i+2][j][k])
                             +8.*(model.gUpper[i+1][j][1] * model.csswm[p].u[i+1][j][k] * model.csswm[p].v[i+1][j][k])
                             -8.*(model.gUpper[i-1][j][1] * model.csswm[p].u[i-1][j][k] * model.csswm[p].v[i-1][j][k])
                             +1.*(model.gUpper[i-2][j][1] * model.csswm[p].u[i-2][j][k] * model.csswm[p].v[i-2][j][k]));

                    pV2_px = 0.5 / (12.*dx_for_u) * 
                            (-1.*(model.gUpper[i+2][j][3] * pow(model.csswm[p].v[i+2][j][k], 2))
                             +8.*(model.gUpper[i+1][j][3] * pow(model.csswm[p].v[i+1][j][k], 2))
                             -8.*(model.gUpper[i-1][j][3] * pow(model.csswm[p].v[i-1][j][k], 2))
                             +1.*(model.gUpper[i-2][j][3] * pow(model.csswm[p].v[i-2][j][k], 2)));

                    rotationU = (((-1.*model.csswm[p].v[i+2][j][k] + 8.*model.csswm[p].v[i+1][j][k] - 8.*model.csswm[p].v[i-1][j][k] + 1.*model.csswm[p].v[i-2][j][k]) / (12.*dx_for_u)) - 
                                 ((-1.*model.csswm[p].u[i][j+2][k] + 8.*model.csswm[p].u[i][j+1][k] - 8.*model.csswm[p].u[i][j-1][k] + 1.*model.csswm[p].u[i][j-2][k]) / (12.*dy_for_u)) + model.sqrtG[i][j] * f) 
                                    * (model.gUpper[i][j][2] * model.csswm[p].u[i][j][k] + model.gUpper[i][j][3] * model.csswm[p].v[i][j][k]);
                

                    model.csswm[p].up[i][j][k] = model.csswm[p].um[i][j][k] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);

                    #if defined(DIFFUSION)
                        model.csswm[p].up[i][j][k] += D2T * KX * (model.csswm[p].um[i+1][j][k] - 2. * model.csswm[p].um[i][j][k] + model.csswm[p].um[i-1][j][k]) / pow(dx_for_u, 2) + 
                                                      D2T * KY * (model.csswm[p].um[i][j+1][k] - 2. * model.csswm[p].um[i][j][k] + model.csswm[p].um[i][j-1][k]) / pow(dy_for_u, 2);
                    #elif defined(DIFFUSION4)
                        model.csswm[p].up[i][j][k] += D2T * KX * (-model.csswm[p].um[i+2][j][k] + 16.*model.csswm[p].um[i+1][j][k] - 30.*model.csswm[p].um[i][j][k] + 16.*model.csswm[p].um[i-1][j][k] - model.csswm[p].um[i-2][j][k]) / (12.*pow(dx_for_u, 2)) + 
                                                      D2T * KY * (-model.csswm[p].um[i][j+2][k] + 16.*model.csswm[p].um[i][j+1][k] - 30.*model.csswm[p].um[i][j][k] + 16.*model.csswm[p].um[i][j-1][k] - model.csswm[p].um[i][j-2][k]) / (12.*pow(dy_for_u, 2));
                    #elif defined(DIFFUSIONLARS)
                        model.csswm[p].up[i][j][k] += D2T * KX * (model.csswm[p].u[i+1][j][k] - model.csswm[p].up[i][j][k] - model.csswm[p].um[i][j][k] + model.csswm[p].u[i-1][j][k]) / pow(dx_for_u, 2) + 
                                                      D2T * KY * (model.csswm[p].u[i][j+1][k] - model.csswm[p].up[i][j][k] - model.csswm[p].um[i][j][k] + model.csswm[p].u[i][j-1][k]) / pow(dy_for_u, 2);
                    #endif
                }
            }
        }
    }
    return;
}

void Iteration::pv_pt_4(CSSWM &model) {
    double dx_for_v = 0, dy_for_v = 0;
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;
    double f;
    #ifdef Mountain
        double pgHs_py = 0.;
    #endif
    double pgH_py_all = 0;

    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                dx_for_v = 0.5 * (model.csswm[p].x[i+1][j] - model.csswm[p].x[i-1][j]);
                dy_for_v = 0.5 * (model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j-1]);

                #if defined(SteadyGeostrophy) || defined(Mountain)
                    f = 2 * OMEGA * (-cos(model.csswm[p].lon[i][j]) * cos(model.csswm[p].lat[i][j]) * sin(ALPHA0) + sin(model.csswm[p].lat[i][j]) * cos(ALPHA0));
                #elif defined(Barotropic) || defined(RossbyHaurwitz)
                    f = 2 * OMEGA * sin(model.csswm[p].lat[i][j]);
                #else
                    f = 0;
                #endif
                pgH_py_all = 0;

                for (int k = 0; k < NZ; k++) {
                    pgH_py_all += GRAVITY / (12.*dy_for_v) * (-1.*model.csswm[p].h[i][j+2][k] + 8.*model.csswm[p].h[i][j+1][k] - 8.*model.csswm[p].h[i][j-1][k] + 1.*model.csswm[p].h[i][j-2][k]);
                }
                #if defined (Mountain)
                    pgH_py_all += GRAVITY / (12.*dy_for_v) * (-1.*model.csswm[p].hs[i][j+2] + 8.*model.csswm[p].hs[i][j+1] - 8.*model.csswm[p].hs[i][j-1] + 1.*model.csswm[p].hs[i][j-2]);
                #endif

                for (int k = 0; k < NZ; k++) {
                    if (k == 0) {
                        #if defined (Mountain)
                            pgH_py = RHO2 / RHO1 * 
                                     (pgH_py_all + (GRAVITY * (RHO1-RHO2)/RHO2) / (12.*dy_for_v) * (-1.*model.csswm[p].h[i][j+2][1] + 8.*model.csswm[p].h[i][j+1][1] - 8.*model.csswm[p].h[i][j-1][1] + 1.*model.csswm[p].h[i][j-2][1] + 
                                                                                                    -1.*model.csswm[p].hs[i][j+2] + 8.*model.csswm[p].hs[i][j+1] - 8.*model.csswm[p].hs[i][j-1] + 1.*model.csswm[p].hs[i][j-2]));
                        #else
                            pgH_py = RHO2 / RHO1 * 
                                     (pgH_py_all + (GRAVITY * (RHO1-RHO2)/RHO2) / (12.*dy_for_v) * (-1.*model.csswm[p].h[i][j+2][1] + 8.*model.csswm[p].h[i][j+1][1] - 8.*model.csswm[p].h[i][j-1][1] + 1.*model.csswm[p].h[i][j-2][1]));
                        #endif
                        
                    }
                    else if (k == 1) {
                        pgH_py = pgH_py_all;
                    }

                    pU2_py = 0.5 /(12.*dy_for_v) * 
                            (-1.*(model.gUpper[i][j+2][0] * pow(model.csswm[p].u[i][j+2][k], 2))
                             +8.*(model.gUpper[i][j+1][0] * pow(model.csswm[p].u[i][j+1][k], 2))
                             -8.*(model.gUpper[i][j-1][0] * pow(model.csswm[p].u[i][j-1][k], 2))
                             +1.*(model.gUpper[i][j-2][0] * pow(model.csswm[p].u[i][j-2][k], 2)));

                    pUV_py = 1. / (12.*dy_for_v) * 
                            (-1.*(model.gUpper[i][j+2][1] * model.csswm[p].u[i][j+2][k] * model.csswm[p].v[i][j+2][k])
                             +8.*(model.gUpper[i][j+1][1] * model.csswm[p].u[i][j+1][k] * model.csswm[p].v[i][j+1][k])
                             -8.*(model.gUpper[i][j-1][1] * model.csswm[p].u[i][j-1][k] * model.csswm[p].v[i][j-1][k])
                             +1.*(model.gUpper[i][j-2][1] * model.csswm[p].u[i][j-2][k] * model.csswm[p].v[i][j-2][k]));

                    pV2_py = 0.5 / (12.*dy_for_v) * 
                            (-1.*(model.gUpper[i][j+2][3] * pow(model.csswm[p].v[i][j+2][k], 2))
                             +8.*(model.gUpper[i][j+1][3] * pow(model.csswm[p].v[i][j+1][k], 2))
                             -8.*(model.gUpper[i][j-1][3] * pow(model.csswm[p].v[i][j-1][k], 2))
                             +1.*(model.gUpper[i][j-2][3] * pow(model.csswm[p].v[i][j-2][k], 2)));

                    rotationV = (((-1.*model.csswm[p].v[i+2][j][k] + 8.*model.csswm[p].v[i+1][j][k] - 8.*model.csswm[p].v[i-1][j][k] + 1.*model.csswm[p].v[i-2][j][k]) / (12.*dx_for_v)) - 
                                 ((-1.*model.csswm[p].u[i][j+2][k] + 8.*model.csswm[p].u[i][j+1][k] - 8.*model.csswm[p].u[i][j-1][k] + 1.*model.csswm[p].u[i][j-2][k]) / (12.*dy_for_v)) + model.sqrtG[i][j] * f) 
                                    * (model.gUpper[i][j][0] * model.csswm[p].u[i][j][k] + model.gUpper[i][j][1] * model.csswm[p].v[i][j][k]);

                    model.csswm[p].vp[i][j][k] = model.csswm[p].vm[i][j][k] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);

                    #if defined(DIFFUSION)
                        model.csswm[p].vp[i][j][k] += D2T * KX * (model.csswm[p].vm[i+1][j][k] - 2. * model.csswm[p].vm[i][j][k] + model.csswm[p].vm[i-1][j][k]) / pow(dx_for_v, 2) + 
                                                      D2T * KY * (model.csswm[p].vm[i][j+1][k] - 2. * model.csswm[p].vm[i][j][k] + model.csswm[p].vm[i][j-1][k]) / pow(dy_for_v, 2);
                    #elif defined(DIFFUSION4)
                        model.csswm[p].vp[i][j][k] += D2T * KX * (-model.csswm[p].vm[i+2][j][k] + 16.*model.csswm[p].vm[i+1][j][k] - 30.*model.csswm[p].vm[i][j][k] + 16.*model.csswm[p].vm[i-1][j][k] - model.csswm[p].vm[i-2][j][k]) / (12.*pow(dx_for_v, 2)) + 
                                                      D2T * KY * (-model.csswm[p].vm[i][j+2][k] + 16.*model.csswm[p].vm[i][j+1][k] - 30.*model.csswm[p].vm[i][j][k] + 16.*model.csswm[p].vm[i][j-1][k] - model.csswm[p].vm[i][j-2][k]) / (12.*pow(dy_for_v, 2));
                    #elif defined(DIFFUSIONLARS)
                        model.csswm[p].vp[i][j][k] += D2T * KX * (model.csswm[p].v[i+1][j][k] - model.csswm[p].vp[i][j][k] - model.csswm[p].vm[i][j][k] + model.csswm[p].v[i-1][j][k]) / pow(dx_for_v, 2) + 
                                                      D2T * KY * (model.csswm[p].v[i][j+1][k] - model.csswm[p].vp[i][j][k] - model.csswm[p].vm[i][j][k] + model.csswm[p].v[i][j-1][k]) / pow(dy_for_v, 2);
                    #endif
                }
            }
        }
    }
    return;
}

#endif

void Iteration::leap_frog(CSSWM &model) {
    Outputs::create_all_directory();
    #ifdef NCOUTPUT
        Outputs::grid_nc(model);
    #endif
    #ifdef TXTOUTPUT
        Outputs::grid(model);
    #endif
    int n = 0;
    // double timenow = 0.;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;

    while (n < nmax) {
        std::cout << n << std::endl;

        if (n % OUTPUTINTERVAL == 0) {
            #ifdef TXTOUTPUT
                Outputs::h(n, model);
                Outputs::u(n, model);
                Outputs::v(n, model);
            #endif

            #ifdef NCOUTPUT
                Outputs::huv_nc(n, model);
            #endif
        }

        n++;
        // timenow = n * DT;

        // Integrate
        #if defined(SecondOrderSpace)
            ph_pt_2(model);
            #ifndef Advection
                pu_pt_2(model);
                pv_pt_2(model);
            #endif

            // Boundary exchange and interpolation
            model.BP_h(model);
            #ifndef Advection
                // model.BP_wind_convert(model);
                // model.BP_wind_interpolation(model);
                model.BP_wind_interpolation2(model);
            #endif
        #elif defined(FourthOrderSpace)
            ph_pt_4(model);
            #ifndef Advection
                pu_pt_4(model);
                pv_pt_4(model);
            #endif

            // Boundary exchange and interpolation
            model.BP_h(model);
            #ifndef Advection
                // model.BP_wind_convert(model);
                // model.BP_wind_interpolation(model);
                model.BP_wind_interpolation2(model);
            #endif
        #endif

        

        // Time filter
        #ifdef TIMEFILTER
            for (int p = 0; p < 6; p++) {
                for (int i = 0; i < NX; i++) {
                    for (int j = 0; j < NY; j++) {
                        for (int k = 0; k < NZ; k++) {
                            model.csswm[p].h[i][j][k] += TIMETS * (model.csswm[p].hp[i][j][k] - 2 * model.csswm[p].h[i][j][k] + model.csswm[p].hm[i][j][k]);
                            model.csswm[p].u[i][j][k] += TIMETS * (model.csswm[p].up[i][j][k] - 2 * model.csswm[p].u[i][j][k] + model.csswm[p].um[i][j][k]);
                            model.csswm[p].v[i][j][k] += TIMETS * (model.csswm[p].vp[i][j][k] - 2 * model.csswm[p].v[i][j][k] + model.csswm[p].vm[i][j][k]);
                        }
                    }
                }
            }
        #endif

        // next step
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    for (int k = 0; k < NZ; k++) {
                        model.csswm[p].hm[i][j][k] = model.csswm[p].h[i][j][k];
                        model.csswm[p].h[i][j][k] = model.csswm[p].hp[i][j][k];

                        model.csswm[p].um[i][j][k] = model.csswm[p].u[i][j][k];
                        model.csswm[p].u[i][j][k] = model.csswm[p].up[i][j][k];

                        model.csswm[p].vm[i][j][k] = model.csswm[p].v[i][j][k];
                        model.csswm[p].v[i][j][k] = model.csswm[p].vp[i][j][k];
                    }
                }
            }
        }
    }
    
    return;
}