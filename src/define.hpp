#ifndef __DEFINE_H__
#define __DEFINE_H__

#define FILLVALUE (-99999999999999);

#define RADIUS (6371220.)
#define GRAVITY (9.80616)
#define OMEGA (7.292E-5)

#define DX (2)
#define DY (2)
#define DT (180.)
#define TIMEEND (86400 * 3 * 2)
#define OUTPUTPATH "/data/Aaron/CSSWM2Layers/Advection/45/0/outputs/"
#define OUTPUTINTERVAL (1)
// #define SecondOrderSpace
#define FourthOrderSpace
#define NCOUTPUT
// #define TXTOUTPUT

#if defined(SecondOrderSpace)
    #define NX ((int) (90/DX + 2))
    #define NY ((int) (90/DY + 2))
#elif defined(FourthOrderSpace) 
    #define NX ((int) (90/DX + 4))
    #define NY ((int) (90/DY + 4))
#endif

#define D2T (2. * DT)

// Jung
#define ALPHA0 (0)
// #define Advection
// #define GravityWave
// #define SteadyGeostrophy
#define Barotropic
// #define Mountain
// #define RossbyHaurwitz

// #define TrueSol
#define DIFFUSION
#define TIMEFILTER
#define KX (200000.)
#define KY (200000.)
#define TIMETS (0.06)

#endif