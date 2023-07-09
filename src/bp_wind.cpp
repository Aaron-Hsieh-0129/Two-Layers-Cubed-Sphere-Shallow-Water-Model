#include "construction.hpp"

void CSSWM::BP_wind_convert(CSSWM &model) {
    for (int k = 0; k < NZ; k++) {
        for (int idx = 0; idx < NX; idx++) {
            #if defined(SecondOrderSpace)
                // patch0 (Left, Right, Up, Down)
                model.csswm[0].up[0][idx][k] = model.Cube2Cube_U(model, 0, 3, 0, idx, NX-2, idx);
                model.csswm[0].vp[0][idx][k] = model.Cube2Cube_V(model, 0, 3, 0, idx, NX-2, idx);

                model.csswm[0].up[NX-1][idx][k] = model.Cube2Cube_U(model, 0, 1, NX-1, idx, 1, idx);
                model.csswm[0].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 0, 1, NX-1, idx, 1, idx);

                model.csswm[0].up[idx][NY-1][k] = model.Cube2Cube_U(model, 0, 4, idx, NY-1, idx, 1);
                model.csswm[0].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 0, 4, idx, NY-1, idx, 1);

                model.csswm[0].up[idx][0][k] = model.Cube2Cube_U(model, 0, 5, idx, 0, idx, NY-2);
                model.csswm[0].vp[idx][0][k] = model.Cube2Cube_V(model, 0, 5, idx, 0, idx, NY-2);

                // patch1 (Left, Right, Up, Down)
                model.csswm[1].up[0][idx][k] = model.Cube2Cube_U(model, 1, 0, 0, idx, NX-2, idx);
                model.csswm[1].vp[0][idx][k] = model.Cube2Cube_V(model, 1, 0, 0, idx, NX-2, idx);

                model.csswm[1].up[NX-1][idx][k] = model.Cube2Cube_U(model, 1, 2, NX-1, idx, 1, idx);
                model.csswm[1].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 1, 2, NX-1, idx, 1, idx);

                model.csswm[1].up[idx][NY-1][k] = model.Cube2Cube_BV2AU(model, 1, 4, idx, NY-1, NX-2, idx);
                model.csswm[1].vp[idx][NY-1][k] = model.Cube2Cube_BU2AV(model, 1, 4, idx, NY-1, NX-2, idx);

                model.csswm[1].up[idx][0][k] = model.Cube2Cube_BV2AU(model, 1, 5, idx, 0, NX-2, NY-1-idx);
                model.csswm[1].vp[idx][0][k] = model.Cube2Cube_BU2AV(model, 1, 5, idx, 0, NX-2, NY-1-idx);

                // patch2 (Left, Right, Up, Down)
                model.csswm[2].up[0][idx][k] = model.Cube2Cube_U(model, 2, 1, 0, idx, NX-2, idx);
                model.csswm[2].vp[0][idx][k] = model.Cube2Cube_V(model, 2, 1, 0, idx, NX-2, idx);

                model.csswm[2].up[NX-1][idx][k] = model.Cube2Cube_U(model, 2, 3, NX-1, idx, 1, idx);
                model.csswm[2].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 2, 3, NX-1, idx, 1, idx);

                model.csswm[2].up[idx][NY-1][k] = model.Cube2Cube_U(model, 2, 4, idx, NY-1, NX-1-idx, NY-2);
                model.csswm[2].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 2, 4, idx, NY-1, NX-1-idx, NY-2);

                model.csswm[2].up[idx][0][k] = model.Cube2Cube_U(model, 2, 5, idx, 0, NX-1-idx, 1);
                model.csswm[2].vp[idx][0][k] = model.Cube2Cube_V(model, 2, 5, idx, 0, NX-1-idx, 1);

                // patch3 (Left, Right, Up, Down)
                model.csswm[3].up[0][idx][k] = model.Cube2Cube_U(model, 3, 2, 0, idx, NX-2, idx);
                model.csswm[3].vp[0][idx][k] = model.Cube2Cube_V(model, 3, 2, 0, idx, NX-2, idx);

                model.csswm[3].up[NX-1][idx][k] = model.Cube2Cube_U(model, 3, 0, NX-1, idx, 1, idx);
                model.csswm[3].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 3, 0, NX-1, idx, 1, idx);

                model.csswm[3].up[idx][NY-1][k] = model.Cube2Cube_BV2AU(model, 3, 4, idx, NY-1, 1, NY-1-idx);
                model.csswm[3].vp[idx][NY-1][k] = model.Cube2Cube_BU2AV(model, 3, 4, idx, NY-1, 1, NY-1-idx);

                model.csswm[3].up[idx][0][k] = model.Cube2Cube_BV2AU(model, 3, 5, idx, 0, 1, idx);
                model.csswm[3].vp[idx][0][k] = model.Cube2Cube_BU2AV(model, 3, 5, idx, 0, 1, idx);

                // patch4 (Left, Right, Up, Down)
                model.csswm[4].up[0][idx][k] = model.Cube2Cube_BV2AU(model, 4, 3, 0, idx, NX-1-idx, NY-2);
                model.csswm[4].vp[0][idx][k] = model.Cube2Cube_BU2AV(model, 4, 3, 0, idx, NX-1-idx, NY-2);

                model.csswm[4].up[NX-1][idx][k] = model.Cube2Cube_BV2AU(model, 4, 1, NX-1, idx, idx, NY-2);
                model.csswm[4].vp[NX-1][idx][k] = model.Cube2Cube_BU2AV(model, 4, 1, NX-1, idx, idx, NY-2);

                model.csswm[4].up[idx][NY-1][k] = model.Cube2Cube_U(model, 4, 2, idx, NY-1, NX-1-idx, NY-2);
                model.csswm[4].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 4, 2, idx, NY-1, NX-1-idx, NY-2);

                model.csswm[4].up[idx][0][k] = model.Cube2Cube_U(model, 4, 0, idx, 0, idx, NY-2);
                model.csswm[4].vp[idx][0][k] = model.Cube2Cube_V(model, 4, 0, idx, 0, idx, NY-2);

                // patch5 (Left, Right, Up, Down)
                model.csswm[5].up[0][idx][k] = model.Cube2Cube_BV2AU(model, 5, 3, 0, idx, idx, 1);
                model.csswm[5].vp[0][idx][k] = model.Cube2Cube_BU2AV(model, 5, 3, 0, idx, idx, 1);

                model.csswm[5].up[NX-1][idx][k] = model.Cube2Cube_BV2AU(model, 5, 1, NX-1, idx, NX-1-idx, 1);
                model.csswm[5].vp[NX-1][idx][k] = model.Cube2Cube_BU2AV(model, 5, 1, NX-1, idx, NX-1-idx, 1);

                model.csswm[5].up[idx][NY-1][k] = model.Cube2Cube_U(model, 5, 0, idx, NY-1, idx, 1);
                model.csswm[5].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 5, 0, idx, NY-1, idx, 1);

                model.csswm[5].up[idx][0][k] = model.Cube2Cube_U(model, 5, 2, idx, 0, NX-1-idx, 1);
                model.csswm[5].vp[idx][0][k] = model.Cube2Cube_V(model, 5, 2, idx, 0, NX-1-idx, 1);
            #elif defined(FourthOrderSpace)
                // patch0 (Left, Right, Up, Down)
                model.csswm[0].up[0][idx][k] = model.Cube2Cube_U(model, 0, 3, 0, idx, NX-4, idx, k);
                model.csswm[0].vp[0][idx][k] = model.Cube2Cube_V(model, 0, 3, 0, idx, NX-4, idx, k);

                model.csswm[0].up[NX-1][idx][k] = model.Cube2Cube_U(model, 0, 1, NX-1, idx, 3, idx, k);
                model.csswm[0].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 0, 1, NX-1, idx, 3, idx, k);

                model.csswm[0].up[idx][NY-1][k] = model.Cube2Cube_U(model, 0, 4, idx, NY-1, idx, 3, k);
                model.csswm[0].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 0, 4, idx, NY-1, idx, 3, k);

                model.csswm[0].up[idx][0][k] = model.Cube2Cube_U(model, 0, 5, idx, 0, idx, NY-4, k);
                model.csswm[0].vp[idx][0][k] = model.Cube2Cube_V(model, 0, 5, idx, 0, idx, NY-4, k);

                model.csswm[0].up[1][idx][k] = model.Cube2Cube_U(model, 0, 3, 0, idx, NX-3, idx, k);
                model.csswm[0].vp[1][idx][k] = model.Cube2Cube_V(model, 0, 3, 0, idx, NX-3, idx, k);

                model.csswm[0].up[NX-2][idx][k] = model.Cube2Cube_U(model, 0, 1, NX-1, idx, 2, idx, k);
                model.csswm[0].vp[NX-2][idx][k] = model.Cube2Cube_V(model, 0, 1, NX-1, idx, 2, idx, k);

                model.csswm[0].up[idx][NY-2][k] = model.Cube2Cube_U(model, 0, 4, idx, NY-1, idx, 2, k);
                model.csswm[0].vp[idx][NY-2][k] = model.Cube2Cube_V(model, 0, 4, idx, NY-1, idx, 2, k);

                model.csswm[0].up[idx][1][k] = model.Cube2Cube_U(model, 0, 5, idx, 0, idx, NY-3, k);
                model.csswm[0].vp[idx][1][k] = model.Cube2Cube_V(model, 0, 5, idx, 0, idx, NY-3, k);


                // patch1 (Left, Right, Up, Down)
                model.csswm[1].up[0][idx][k] = model.Cube2Cube_U(model, 1, 0, 0, idx, NX-4, idx, k);
                model.csswm[1].vp[0][idx][k] = model.Cube2Cube_V(model, 1, 0, 0, idx, NX-4, idx, k);

                model.csswm[1].up[NX-1][idx][k] = model.Cube2Cube_U(model, 1, 2, NX-1, idx, 3, idx, k);
                model.csswm[1].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 1, 2, NX-1, idx, 3, idx, k);

                model.csswm[1].up[idx][NY-1][k] = model.Cube2Cube_BV2AU(model, 1, 4, idx, NY-1, NX-4, idx, k);
                model.csswm[1].vp[idx][NY-1][k] = model.Cube2Cube_BU2AV(model, 1, 4, idx, NY-1, NX-4, idx, k);

                model.csswm[1].up[idx][0][k] = model.Cube2Cube_BV2AU(model, 1, 5, idx, 0, NX-4, NY-1-idx, k);
                model.csswm[1].vp[idx][0][k] = model.Cube2Cube_BU2AV(model, 1, 5, idx, 0, NX-4, NY-1-idx, k);

                model.csswm[1].up[1][idx][k] = model.Cube2Cube_U(model, 1, 0, 0, idx, NX-3, idx, k);
                model.csswm[1].vp[1][idx][k] = model.Cube2Cube_V(model, 1, 0, 0, idx, NX-3, idx, k);

                model.csswm[1].up[NX-2][idx][k] = model.Cube2Cube_U(model, 1, 2, NX-1, idx, 2, idx, k);
                model.csswm[1].vp[NX-2][idx][k] = model.Cube2Cube_V(model, 1, 2, NX-1, idx, 2, idx, k);

                model.csswm[1].up[idx][NY-2][k] = model.Cube2Cube_BV2AU(model, 1, 4, idx, NY-1, NX-3, idx, k);
                model.csswm[1].vp[idx][NY-2][k] = model.Cube2Cube_BU2AV(model, 1, 4, idx, NY-1, NX-3, idx, k);

                model.csswm[1].up[idx][1][k] = model.Cube2Cube_BV2AU(model, 1, 5, idx, 0, NX-3, NY-1-idx, k);
                model.csswm[1].vp[idx][1][k] = model.Cube2Cube_BU2AV(model, 1, 5, idx, 0, NX-3, NY-1-idx, k);

                // patch2 (Left, Right, Up, Down)
                model.csswm[2].up[0][idx][k] = model.Cube2Cube_U(model, 2, 1, 0, idx, NX-4, idx, k);
                model.csswm[2].vp[0][idx][k] = model.Cube2Cube_V(model, 2, 1, 0, idx, NX-4, idx, k);

                model.csswm[2].up[NX-1][idx][k] = model.Cube2Cube_U(model, 2, 3, NX-1, idx, 3, idx, k);
                model.csswm[2].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 2, 3, NX-1, idx, 3, idx, k);

                model.csswm[2].up[idx][NY-1][k] = model.Cube2Cube_U(model, 2, 4, idx, NY-1, NX-1-idx, NY-4, k);
                model.csswm[2].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 2, 4, idx, NY-1, NX-1-idx, NY-4, k);

                model.csswm[2].up[idx][0][k] = model.Cube2Cube_U(model, 2, 5, idx, 0, NX-1-idx, 3, k);
                model.csswm[2].vp[idx][0][k] = model.Cube2Cube_V(model, 2, 5, idx, 0, NX-1-idx, 3, k);

                model.csswm[2].up[1][idx][k] = model.Cube2Cube_U(model, 2, 1, 0, idx, NX-3, idx, k);
                model.csswm[2].vp[1][idx][k] = model.Cube2Cube_V(model, 2, 1, 0, idx, NX-3, idx, k);

                model.csswm[2].up[NX-2][idx][k] = model.Cube2Cube_U(model, 2, 3, NX-1, idx, 2, idx, k);
                model.csswm[2].vp[NX-2][idx][k] = model.Cube2Cube_V(model, 2, 3, NX-1, idx, 2, idx, k);

                model.csswm[2].up[idx][NY-2][k] = model.Cube2Cube_U(model, 2, 4, idx, NY-1, NX-1-idx, NY-3, k);
                model.csswm[2].vp[idx][NY-2][k] = model.Cube2Cube_V(model, 2, 4, idx, NY-1, NX-1-idx, NY-3, k);

                model.csswm[2].up[idx][1][k] = model.Cube2Cube_U(model, 2, 5, idx, 0, NX-1-idx, 2, k);
                model.csswm[2].vp[idx][1][k] = model.Cube2Cube_V(model, 2, 5, idx, 0, NX-1-idx, 2, k);

                // patch3 (Left, Right, Up, Down)
                model.csswm[3].up[0][idx][k] = model.Cube2Cube_U(model, 3, 2, 0, idx, NX-4, idx, k);
                model.csswm[3].vp[0][idx][k] = model.Cube2Cube_V(model, 3, 2, 0, idx, NX-4, idx, k);

                model.csswm[3].up[NX-1][idx][k] = model.Cube2Cube_U(model, 3, 0, NX-1, idx, 3, idx, k);
                model.csswm[3].vp[NX-1][idx][k] = model.Cube2Cube_V(model, 3, 0, NX-1, idx, 3, idx, k);

                model.csswm[3].up[idx][NY-1][k] = model.Cube2Cube_BV2AU(model, 3, 4, idx, NY-1, 3, NY-1-idx, k);
                model.csswm[3].vp[idx][NY-1][k] = model.Cube2Cube_BU2AV(model, 3, 4, idx, NY-1, 3, NY-1-idx, k);

                model.csswm[3].up[idx][0][k] = model.Cube2Cube_BV2AU(model, 3, 5, idx, 0, 3, idx, k);
                model.csswm[3].vp[idx][0][k] = model.Cube2Cube_BU2AV(model, 3, 5, idx, 0, 3, idx, k);

                model.csswm[3].up[1][idx][k] = model.Cube2Cube_U(model, 3, 2, 0, idx, NX-3, idx, k);
                model.csswm[3].vp[1][idx][k] = model.Cube2Cube_V(model, 3, 2, 0, idx, NX-3, idx, k);

                model.csswm[3].up[NX-2][idx][k] = model.Cube2Cube_U(model, 3, 0, NX-1, idx, 2, idx, k);
                model.csswm[3].vp[NX-2][idx][k] = model.Cube2Cube_V(model, 3, 0, NX-1, idx, 2, idx, k);

                model.csswm[3].up[idx][NY-2][k] = model.Cube2Cube_BV2AU(model, 3, 4, idx, NY-1, 2, NY-1-idx, k);
                model.csswm[3].vp[idx][NY-2][k] = model.Cube2Cube_BU2AV(model, 3, 4, idx, NY-1, 2, NY-1-idx, k);

                model.csswm[3].up[idx][1][k] = model.Cube2Cube_BV2AU(model, 3, 5, idx, 0, 2, idx, k);
                model.csswm[3].vp[idx][1][k] = model.Cube2Cube_BU2AV(model, 3, 5, idx, 0, 2, idx, k);

                // patch4 (Left, Right, Up, Down)
                model.csswm[4].up[0][idx][k] = model.Cube2Cube_BV2AU(model, 4, 3, 0, idx, NX-1-idx, NY-4, k);
                model.csswm[4].vp[0][idx][k] = model.Cube2Cube_BU2AV(model, 4, 3, 0, idx, NX-1-idx, NY-4, k);

                model.csswm[4].up[NX-1][idx][k] = model.Cube2Cube_BV2AU(model, 4, 1, NX-1, idx, idx, NY-4, k);
                model.csswm[4].vp[NX-1][idx][k] = model.Cube2Cube_BU2AV(model, 4, 1, NX-1, idx, idx, NY-4, k);

                model.csswm[4].up[idx][NY-1][k] = model.Cube2Cube_U(model, 4, 2, idx, NY-1, NX-1-idx, NY-4, k);
                model.csswm[4].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 4, 2, idx, NY-1, NX-1-idx, NY-4, k);

                model.csswm[4].up[idx][0][k] = model.Cube2Cube_U(model, 4, 0, idx, 0, idx, NY-4, k);
                model.csswm[4].vp[idx][0][k] = model.Cube2Cube_V(model, 4, 0, idx, 0, idx, NY-4, k);

                model.csswm[4].up[1][idx][k] = model.Cube2Cube_BV2AU(model, 4, 3, 0, idx, NX-1-idx, NY-3, k);
                model.csswm[4].vp[1][idx][k] = model.Cube2Cube_BU2AV(model, 4, 3, 0, idx, NX-1-idx, NY-3, k);

                model.csswm[4].up[NX-2][idx][k] = model.Cube2Cube_BV2AU(model, 4, 1, NX-1, idx, idx, NY-3, k);
                model.csswm[4].vp[NX-2][idx][k] = model.Cube2Cube_BU2AV(model, 4, 1, NX-1, idx, idx, NY-3, k);

                model.csswm[4].up[idx][NY-2][k] = model.Cube2Cube_U(model, 4, 2, idx, NY-1, NX-1-idx, NY-3, k);
                model.csswm[4].vp[idx][NY-2][k] = model.Cube2Cube_V(model, 4, 2, idx, NY-1, NX-1-idx, NY-3, k);

                model.csswm[4].up[idx][1][k] = model.Cube2Cube_U(model, 4, 0, idx, 0, idx, NY-3, k);
                model.csswm[4].vp[idx][1][k] = model.Cube2Cube_V(model, 4, 0, idx, 0, idx, NY-3, k);

                // patch5 (Left, Right, Up, Down)
                model.csswm[5].up[0][idx][k] = model.Cube2Cube_BV2AU(model, 5, 3, 0, idx, idx, 3, k);
                model.csswm[5].vp[0][idx][k] = model.Cube2Cube_BU2AV(model, 5, 3, 0, idx, idx, 3, k);

                model.csswm[5].up[NX-1][idx][k] = model.Cube2Cube_BV2AU(model, 5, 1, NX-1, idx, NX-1-idx, 3, k);
                model.csswm[5].vp[NX-1][idx][k] = model.Cube2Cube_BU2AV(model, 5, 1, NX-1, idx, NX-1-idx, 3, k);

                model.csswm[5].up[idx][NY-1][k] = model.Cube2Cube_U(model, 5, 0, idx, NY-1, idx, 3, k);
                model.csswm[5].vp[idx][NY-1][k] = model.Cube2Cube_V(model, 5, 0, idx, NY-1, idx, 3, k);

                model.csswm[5].up[idx][0][k] = model.Cube2Cube_U(model, 5, 2, idx, 0, NX-1-idx, 3, k);
                model.csswm[5].vp[idx][0][k] = model.Cube2Cube_V(model, 5, 2, idx, 0, NX-1-idx, 3, k);

                model.csswm[5].up[1][idx][k] = model.Cube2Cube_BV2AU(model, 5, 3, 0, idx, idx, 2, k);
                model.csswm[5].vp[1][idx][k] = model.Cube2Cube_BU2AV(model, 5, 3, 0, idx, idx, 2, k);

                model.csswm[5].up[NX-2][idx][k] = model.Cube2Cube_BV2AU(model, 5, 1, NX-1, idx, NX-1-idx, 2, k);
                model.csswm[5].vp[NX-2][idx][k] = model.Cube2Cube_BU2AV(model, 5, 1, NX-1, idx, NX-1-idx, 2, k);

                model.csswm[5].up[idx][NY-2][k] = model.Cube2Cube_U(model, 5, 0, idx, NY-1, idx, 2, k);
                model.csswm[5].vp[idx][NY-2][k] = model.Cube2Cube_V(model, 5, 0, idx, NY-1, idx, 2, k);

                model.csswm[5].up[idx][1][k] = model.Cube2Cube_U(model, 5, 2, idx, 0, NX-1-idx, 2, k);
                model.csswm[5].vp[idx][1][k] = model.Cube2Cube_V(model, 5, 2, idx, 0, NX-1-idx, 2, k);
            #endif
        }
    }
}


void CSSWM::BP_wind_interpolation(CSSWM &model) {

    double B, A1, A2, V1, V2, V3, V4;

    int p1, p2, i1, j1, i2, j2, reversed, lonlat;
    double uIP, vIP;
    double alpha, beta;
    double alpha_B, beta_B;
    double gLower[4], IA[4], A[4], gUpper[4];
    for (int k = 0; k < NZ; k++) {
        for (int pp = 0; pp < 24; pp++) {
            #if defined(SecondOrderSpace)
                p1 = model.match[pp][0], p2 = model.match[pp][1], i1 = model.match[pp][2], j1 = model.match[pp][3], i2 = model.match[pp][4], j2 = model.match[pp][5], reversed = model.match[pp][6], lonlat = model.match[pp][7];
            #elif defined(FourthOrderSpace)
                p1 = model.match_ouTTer[pp][0], p2 = model.match_ouTTer[pp][1], i1 = model.match_ouTTer[pp][2], j1 = model.match_ouTTer[pp][3], i2 = model.match_ouTTer[pp][4], j2 = model.match_ouTTer[pp][5], reversed = model.match_ouTTer[pp][6], lonlat = model.match_ouTTer[pp][7];
            #endif
            for (int idx = 0; idx < NX; idx++) {
                int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                #if defined(SecondOrderSpace)
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP[NX-1-idx][0] : model.checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP[NY-1-idx][0] : model.checkIP[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP[NX-1-idx][1] : model.checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP[NY-1-idx][1] : model.checkIP[idx][1] : j2;
                #elif defined(FourthOrderSpace)
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTTer[NX-1-idx][0] : model.checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTTer[NY-1-idx][0] : model.checkIP_ouTTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTTer[NX-1-idx][1] : model.checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTTer[NY-1-idx][1] : model.checkIP_ouTTer[idx][1] : j2;
                #endif
                if (lonlat == 0) {
                    B = model.csswm[p1].lat[I1][J1];
                    A1 = model.csswm[p2].lat[I2_1][J2_1], A2 = model.csswm[p2].lat[I2_2][J2_2];
                    V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                    V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];

                    alpha = model.alpha2D[I1][J1];
                    beta = model.beta2D[I1][J1];

                    uIP = model.interpolate(A1, A2, V1, V2, B);
                    vIP = model.interpolate(A1, A2, V3, V4, B);

                    alpha_B = model.interpolate(A1, A2, model.alpha2D[I2_1][J2_1], model.alpha2D[I2_2][J2_2], B);
                    beta_B = model.interpolate(A1, A2, model.beta2D[I2_1][J2_1], model.beta2D[I2_2][J2_2], B);

                    model.get_gLower(gLower, alpha, beta);
                    model.get_IA(IA, p1, alpha, beta);
                    model.get_A(A, p2, alpha_B, beta_B);
                    model.get_gUpper(gUpper, alpha_B, beta_B);
                }
                else {
                    B = model.csswm[p1].lon[I1][J1];
                    A1 = model.csswm[p2].lon[I2_1][J2_1], A2 = model.csswm[p2].lon[I2_2][J2_2];
                    V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                    V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                    
                    alpha = model.alpha2D[I1][J1];
                    beta = model.beta2D[I1][J1];

                    uIP = model.interpolate(A1, A2, V1, V2, B);
                    vIP = model.interpolate(A1, A2, V3, V4, B);

                    alpha_B = model.interpolate(A1, A2, model.alpha2D[I2_1][J2_1], model.alpha2D[I2_2][J2_2], B);
                    beta_B = model.interpolate(A1, A2, model.beta2D[I2_1][J2_1], model.beta2D[I2_2][J2_2], B);

                    model.get_gLower(gLower, alpha, beta);
                    model.get_IA(IA, p1, alpha, beta);
                    model.get_A(A, p2, alpha_B, beta_B);
                    model.get_gUpper(gUpper, alpha_B, beta_B);
                }
                model.csswm[p1].up[I1][J1][k] = model.Cube2Cube_U_2(gLower, IA, A, gUpper, uIP, vIP);
                model.csswm[p1].vp[I1][J1][k] = model.Cube2Cube_V_2(gLower, IA, A, gUpper, uIP, vIP);

                // std::cout << "(" << I2_1 << ", " << J2_1 << ")" << ", (" << I2_2 << ", " << J2_2 << ")" << std::endl;
            }
        }

        #if defined(FourthOrderSpace)
            for (int pp = 0; pp < 24; pp++) {
                p1 = model.match_ouTer[pp][0], p2 = model.match_ouTer[pp][1], i1 = model.match_ouTer[pp][2], j1 = model.match_ouTer[pp][3], i2 = model.match_ouTer[pp][4], j2 = model.match_ouTer[pp][5], reversed = model.match_ouTer[pp][6], lonlat = model.match_ouTer[pp][7];

                for (int idx = 0; idx < NX; idx++) {
                    int I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                    int I2_1 = i2 == -1 ? reversed ? model.checkIP_ouTer[NX-1-idx][0] : model.checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? model.checkIP_ouTer[NY-1-idx][0] : model.checkIP_ouTer[idx][0] : j2;
                    int I2_2 = i2 == -1 ? reversed ? model.checkIP_ouTer[NX-1-idx][1] : model.checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? model.checkIP_ouTer[NY-1-idx][1] : model.checkIP_ouTer[idx][1] : j2;
                    if (lonlat == 0) {
                        B = model.csswm[p1].lat[I1][J1];
                        A1 = model.csswm[p2].lat[I2_1][J2_1], A2 = model.csswm[p2].lat[I2_2][J2_2];
                        V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                        V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];

                        alpha = model.alpha2D[I1][J1];
                        beta = model.beta2D[I1][J1];

                        uIP = model.interpolate(A1, A2, V1, V2, B);
                        vIP = model.interpolate(A1, A2, V3, V4, B);

                        alpha_B = model.interpolate(A1, A2, model.alpha2D[I2_1][J2_1], model.alpha2D[I2_2][J2_2], B);
                        beta_B = model.interpolate(A1, A2, model.beta2D[I2_1][J2_1], model.beta2D[I2_2][J2_2], B);

                        model.get_gLower(gLower, alpha, beta);
                        model.get_IA(IA, p1, alpha, beta);
                        model.get_A(A, p2, alpha_B, beta_B);
                        model.get_gUpper(gUpper, alpha_B, beta_B);
                    }
                    else {
                        B = model.csswm[p1].lon[I1][J1];
                        A1 = model.csswm[p2].lon[I2_1][J2_1], A2 = model.csswm[p2].lon[I2_2][J2_2];
                        V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                        V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];

                        if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                        if (A1 > B && B < A2) B += 2 * M_PI;
                        
                        alpha = model.alpha2D[I1][J1];
                        beta = model.beta2D[I1][J1];

                        uIP = model.interpolate(A1, A2, V1, V2, B);
                        vIP = model.interpolate(A1, A2, V3, V4, B);

                        alpha_B = model.interpolate(A1, A2, model.alpha2D[I2_1][J2_1], model.alpha2D[I2_2][J2_2], B);
                        beta_B = model.interpolate(A1, A2, model.beta2D[I2_1][J2_1], model.beta2D[I2_2][J2_2], B);

                        model.get_gLower(gLower, alpha, beta);
                        model.get_IA(IA, p1, alpha, beta);
                        model.get_A(A, p2, alpha_B, beta_B);
                        model.get_gUpper(gUpper, alpha_B, beta_B);
                    }
                    model.csswm[p1].up[I1][J1][k] = model.Cube2Cube_U_2(gLower, IA, A, gUpper, uIP, vIP);
                    model.csswm[p1].vp[I1][J1][k] = model.Cube2Cube_V_2(gLower, IA, A, gUpper, uIP, vIP);
                }
            }
        #endif
    }
}


void CSSWM::BP_wind_interpolation2(CSSWM &model) {
    double B, A1, A2, V1, V2, V3, V4;

    int p1, p2, i1, j1, i2, j2, reversed, lonlat;
    double uIP, vIP;
    int I1, I2_1, I2_2, J1, J2_1, J2_2;
    for (int k = 0; k < NZ; k++) {
        for (int pp = 0; pp < 24; pp++) {
            #if defined(SecondOrderSpace)
                p1 = match[pp][0], p2 = match[pp][1], i1 = match[pp][2], j1 = match[pp][3], i2 = match[pp][4], j2 = match[pp][5], reversed = match[pp][6], lonlat = match[pp][7];
            #elif defined(FourthOrderSpace)
                p1 = match_ouTTer[pp][0], p2 = match_ouTTer[pp][1], i1 = match_ouTTer[pp][2], j1 = match_ouTTer[pp][3], i2 = match_ouTTer[pp][4], j2 = match_ouTTer[pp][5], reversed = match_ouTTer[pp][6], lonlat = match_ouTTer[pp][7];
            #endif
            for (int idx = 0; idx < NX; idx++) {
                I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                #if defined(SecondOrderSpace)
                    I2_1 = i2 == -1 ? reversed ? checkIP[NX-1-idx][0] : checkIP[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP[NY-1-idx][0] : checkIP[idx][0] : j2;
                    I2_2 = i2 == -1 ? reversed ? checkIP[NX-1-idx][1] : checkIP[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP[NY-1-idx][1] : checkIP[idx][1] : j2;
                #elif defined(FourthOrderSpace)
                    I2_1 = i2 == -1 ? reversed ? checkIP_ouTTer[NX-1-idx][0] : checkIP_ouTTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP_ouTTer[NY-1-idx][0] : checkIP_ouTTer[idx][0] : j2;
                    I2_2 = i2 == -1 ? reversed ? checkIP_ouTTer[NX-1-idx][1] : checkIP_ouTTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP_ouTTer[NY-1-idx][1] : checkIP_ouTTer[idx][1] : j2;
                #endif
                if (lonlat == 0) {
                    B = model.csswm[p1].lat[I1][J1];
                    A1 = model.csswm[p2].lat[I2_1][J2_1], A2 = model.csswm[p2].lat[I2_2][J2_2];
                    V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                    V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];
                }
                else {
                    B = model.csswm[p1].lon[I1][J1];
                    A1 = model.csswm[p2].lon[I2_1][J2_1], A2 = model.csswm[p2].lon[I2_2][J2_2];
                    V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                    V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];

                    if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                    if (A1 > B && B < A2) B += 2 * M_PI;
                }

                uIP = model.interpolate(A1, A2, V1, V2, B);
                vIP = model.interpolate(A1, A2, V3, V4, B);

                if (i1 == -1) {
                    if (j1 == 0) {
                        #if defined(SecondOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP1_D[I1][0] * uIP + model.csswm[p1].IP1_D[I1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP1_D[I1][2] * uIP + model.csswm[p1].IP1_D[I1][3] * vIP;
                        #elif defined(FourthOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTTer_D[I1][0] * uIP + model.csswm[p1].IP_ouTTer_D[I1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTTer_D[I1][2] * uIP + model.csswm[p1].IP_ouTTer_D[I1][3] * vIP;
                        #endif
                    }
                    else if (j1 == NY-1) {
                        #if defined(SecondOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP1_U[I1][0] * uIP + model.csswm[p1].IP1_U[I1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP1_U[I1][2] * uIP + model.csswm[p1].IP1_U[I1][3] * vIP;
                        #elif defined(FourthOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTTer_U[I1][0] * uIP + model.csswm[p1].IP_ouTTer_U[I1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTTer_U[I1][2] * uIP + model.csswm[p1].IP_ouTTer_U[I1][3] * vIP;
                        #endif
                    }
                }
                else if (j1 == -1) {
                    if (i1 == 0) {
                        #if defined(SecondOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP1_L[J1][0] * uIP + model.csswm[p1].IP1_L[J1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP1_L[J1][2] * uIP + model.csswm[p1].IP1_L[J1][3] * vIP;
                        #elif defined(FourthOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTTer_L[J1][0] * uIP + model.csswm[p1].IP_ouTTer_L[J1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTTer_L[J1][2] * uIP + model.csswm[p1].IP_ouTTer_L[J1][3] * vIP;
                        #endif
                    }
                    else if (i1 == NX-1) {
                        #if defined(SecondOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP1_R[J1][0] * uIP + model.csswm[p1].IP1_R[J1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP1_R[J1][2] * uIP + model.csswm[p1].IP1_R[J1][3] * vIP;
                        #elif defined(FourthOrderSpace)
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTTer_R[J1][0] * uIP + model.csswm[p1].IP_ouTTer_R[J1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTTer_R[J1][2] * uIP + model.csswm[p1].IP_ouTTer_R[J1][3] * vIP;
                        #endif
                    }
                }
            }
        }

        #if defined(FourthOrderSpace)
            for (int pp = 0; pp < 24; pp++) {
                p1 = match_ouTer[pp][0], p2 = match_ouTer[pp][1], i1 = match_ouTer[pp][2], j1 = match_ouTer[pp][3], i2 = match_ouTer[pp][4], j2 = match_ouTer[pp][5], reversed = match_ouTer[pp][6], lonlat = match_ouTer[pp][7];
                for (int idx = 0; idx < NX; idx++) {
                    if (lonlat == 0) {
                        I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                        I2_1 = i2 == -1 ? reversed ? checkIP_ouTer[NX-1-idx][0] : checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP_ouTer[NY-1-idx][0] : checkIP_ouTer[idx][0] : j2;
                        I2_2 = i2 == -1 ? reversed ? checkIP_ouTer[NX-1-idx][1] : checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP_ouTer[NY-1-idx][1] : checkIP_ouTer[idx][1] : j2;

                        B = model.csswm[p1].lat[I1][J1];
                        A1 = model.csswm[p2].lat[I2_1][J2_1], A2 = model.csswm[p2].lat[I2_2][J2_2];
                        V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                        V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];
                    }
                    else {
                        I1 = i1 == -1 ? idx : i1, J1 = j1 == -1 ? idx : j1;
                        I2_1 = i2 == -1 ? reversed ? checkIP_ouTer[NX-1-idx][0] : checkIP_ouTer[idx][0] : i2, J2_1 = j2 == -1 ? reversed ? checkIP_ouTer[NY-1-idx][0] : checkIP_ouTer[idx][0] : j2;
                        I2_2 = i2 == -1 ? reversed ? checkIP_ouTer[NX-1-idx][1] : checkIP_ouTer[idx][1] : i2, J2_2 = j2 == -1 ? reversed ? checkIP_ouTer[NY-1-idx][1] : checkIP_ouTer[idx][1] : j2;

                        B = model.csswm[p1].lon[I1][J1];
                        A1 = model.csswm[p2].lon[I2_1][J2_1], A2 = model.csswm[p2].lon[I2_2][J2_2];
                        V1 = model.csswm[p2].up[I2_1][J2_1][k], V2 = model.csswm[p2].up[I2_2][J2_2][k];
                        V3 = model.csswm[p2].vp[I2_1][J2_1][k], V4 = model.csswm[p2].vp[I2_2][J2_2][k];

                        if (A1 > A2 && (p1 == 0 || p2 == 0))  A2 += 2 * M_PI;
                        if (A1 > B && B < A2) B += 2 * M_PI;
                    }

                    uIP = model.interpolate(A1, A2, V1, V2, B);
                    vIP = model.interpolate(A1, A2, V3, V4, B);

                    if (i1 == -1) {
                        if (j1 == 1) {
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTer_D[I1][0] * uIP + model.csswm[p1].IP_ouTer_D[I1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTer_D[I1][2] * uIP + model.csswm[p1].IP_ouTer_D[I1][3] * vIP;
                        }
                        else if (j1 == NY-2) {
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTer_U[I1][0] * uIP + model.csswm[p1].IP_ouTer_U[I1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTer_U[I1][2] * uIP + model.csswm[p1].IP_ouTer_U[I1][3] * vIP;
                        }
                    }
                    else if (j1 == -1) {
                        if (i1 == 1) {
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTer_L[J1][0] * uIP + model.csswm[p1].IP_ouTer_L[J1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTer_L[J1][2] * uIP + model.csswm[p1].IP_ouTer_L[J1][3] * vIP;
                        }
                        else if (i1 == NX-2) {
                            model.csswm[p1].up[I1][J1][k] = model.csswm[p1].IP_ouTer_R[J1][0] * uIP + model.csswm[p1].IP_ouTer_R[J1][1] * vIP;
                            model.csswm[p1].vp[I1][J1][k] = model.csswm[p1].IP_ouTer_R[J1][2] * uIP + model.csswm[p1].IP_ouTer_R[J1][3] * vIP;
                        }
                    }
                }
            }
        #endif
    }
}
