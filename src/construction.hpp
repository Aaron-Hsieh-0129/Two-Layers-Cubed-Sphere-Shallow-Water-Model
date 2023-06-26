#include <cmath>
#include <iostream>
#include <iomanip>
#include "define.hpp"

class CSSWM {
public:
    class patch {
    public: 
        // constructer
        patch();

        double hp[NX][NY], h[NX][NY], hm[NX][NY];
        double up[NX][NY], u[NX][NY], um[NX][NY];
        double vp[NX][NY], v[NX][NY], vm[NX][NY];

        #ifdef Mountain
            double hs[NX][NY];
        #endif

        double lon[NX][NY], lat[NX][NY];
        double lon_original[NX][NY];

        double x[NX][NY], y[NX][NY];

        double A[NX][NY][4], IA[NX][NY][4];

        #if defined(SecondOrderSpace)
            double IP1_L[NX][4], IP1_R[NX][4], IP1_U[NX][4], IP1_D[NX][4]; 
        #elif defined(FourthOrderSpace)
            double IP_ouTTer_L[NX][4], IP_ouTTer_R[NX][4], IP_ouTTer_U[NX][4], IP_ouTTer_D[NX][4]; 
            double IP_ouTer_L[NX][4], IP_ouTer_R[NX][4], IP_ouTer_U[NX][4], IP_ouTer_D[NX][4]; 
        #endif
    };

    CSSWM();
    patch csswm[6];
    double sqrtG[NX][NY], gamma[NX][NY], gLower[NX][NY][4], gUpper[NX][NY][4];
    double alpha2D[NX][NY], beta2D[NX][NY];
    
    #if defined(SecondOrderSpace)
        int match[24][8];
        int checkIP[NX][2];
    #elif defined(FourthOrderSpace)
        int match_ouTTer[24][8];
        int match_ouTer[24][8];
        int checkIP_ouTTer[NX][2];
        int checkIP_ouTer[NX][2];
    #endif

    // ***********************************************************************************
    // In construction.cpp
    void get_gUpper(double ans[4], double alpha, double beta);
    void get_gLower(double ans[4], double alpha, double beta);
    void get_A(double ans[4], int p, double alpha, double beta);
    void get_IA(double ans[4], int p, double alpha, double beta);
    // ***********************************************************************************

    // ***********************************************************************************
    // In transform.cpp
    double Cube2Sphere_U(CSSWM &, int, int, int);
    double Cube2Sphere_V(CSSWM &, int, int, int);
    double Sphere2Cube_U(CSSWM &, int, int, int);
    double Sphere2Cube_V(CSSWM &, int, int, int);
    double Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_V(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_BV2AU(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_BU2AV(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2);
    double Cube2Cube_U_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v);
    double Cube2Cube_V_2(double gLower[4], double IA[4], double A[4], double gUpper[4], double u, double v);
    void matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]);
    void Cube2Cube_matrix();
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_h.cpp
    void BP_h(CSSWM &);
    #if defined(Mountain)
        void BP_hs(CSSWM &);
    #endif
    double interpolate(double A1, double A2, double V1, double V2, double B);
    // ***********************************************************************************

    // ***********************************************************************************
    // In bp_wind.cpp
    void BP_wind_convert(CSSWM &);
    void BP_wind_interpolation(CSSWM &);
    void BP_wind_interpolation2(CSSWM &);
    // ***********************************************************************************



private:
    void Construct_gamma_sqrtG_GUpper(double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double sqrtG[NX][NY], double gUpper[NX][NY][4], double gLower[NX][NY][4]);
    void Construct_p0123_lonlat_xy_AIA(int p, double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
    void Construct_p4_lonlat_xy_AIA(int p, double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
    void Construct_p5_lonlat_xy_AIA(int p, double alpha2D[NX][NY], double beta2D[NX][NY], double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
};