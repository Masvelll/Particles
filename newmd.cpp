#include <iostream>

#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * (x) * (x))
#define DO_MOL for (n=0; n < nMol; n++)

#define VAdd(v1, v2, v3) (v1).x = (v2).x + (v3).x, (v1).y = (v2).y + (v3).y
#define VSub(v1, v2, v3) (v1).x = (v2).x - (v3).x, (v1).y = (v2).y - (v3).y
#define VMul(v1, v2, v3) (v1).x = (v2).x * (v3).x, (v1).y = (v2).y * (v3).y
#define VDiv(v1, v2, v3) (v1).x = (v2).x / (v3).x, (v1).y = (v2).y / (v3).y

#define VDot(v1, v2) ((v1).x * (v2).x + (v1).y * (v2).y)
#define VSAdd(v1, v2, s3, v3) ((v1).x = (v2).x + (s3) * (v3).x), ((v1).y = (v2).y + (s3) * (v3).y)
#define VSet(v, sx, sy) (v).x = sx, (v).y = sy

#define VSetAll(v, s) VSet (v, s, s)
#define VZero(v) VSetAll (v, 0)
#define VVSAdd(v1, s2, v2) VSAdd (v1, v2, s2, v2)
#define VLenSq(v) VDot (v, v)

#define VScale(v, s) (v).x *= s, (v).y *= s
#define VVAdd(v1, v2) VAdd(v1, v1, v2)


#define VWrap(v, t)                             \
    if (v.t >= 0.5 region.t) v.t -= region.t    \
    else if (v.t < -0.5 * region.t) v.t += region.t

#define VWrapAll(v) \ 
    {VWrap (v, x);  \
     VWrap (v, y);}
typedef double real;

typedef struct{
    real x, y;
} VecR;

typedef struct{
    int64_t x, y;
} VecI;

typedef struct{
    VecR r, rv, ra;
} Mol;

typedef struct {
    real val, sum, sum2;
} Prop;

Mol *mol;
VecR region, vSum;
VecI initUcell;
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum;
int16_t moreCycles, nMol, stepAvg, stepCount, stepEquilm, stepLimit;

void InitCoords (){
    VecR c, gap;
    int16_t n, nx, ny;

    VDiv (gap, region, initUcell);
    n = 0;
    for (ny = 0; ny < initUcell.y; ny ++){
        for (nx = 0; nx < initUcell.x; nx ++) {
            VSet (c, nx + 0.5, ny + 0.5);
            VMul (c, c, gap);
            VVSAdd (c, -0.5, region);
            mol[n].r = c;
            ++ n;
        }
    }
}

void InitVels () {
    int16_t n;

    VZero (vSum);
    DO_MOL {
        //VRand (&mol[n].rv);
        VScale (mol[n].rv, velMag);
        VVAdd (vSum, mol[n].rv);
    }
    DO_MOL VVSAdd (mol[n].rv, -1. / nMol, vSum);
}

void InitAccels () {
    int16_t n;

    DO_MOL VZero (mol[n].ra);
}

void SetupJob(){
    AllocArrays();
    stepCount = 0;
    InitCoords();
    InitVels();
    InitAccels();
    AccumProps(0);
}

void ComputeForces(){
    VecR dr;
    real fcVal, rr, rrCut, rri, rri3;
    int16_t j1, j2, n;

    rrCut = Sqr(rCut);
    DO_MOL VZero (mol[n].ra);
    uSum = 0.;
    for (j1 = 0; j1 <nMol -1; j1 ++){
        for (j2 = j1 + 1; j2 < nMol; j2 ++) {
            VSub (dr, mol[j1].r, mol[j2].r);
            VWrapAll (dr);
            rr = VLenSq (dr);
            if (rr < rrCut) {
                rri = 1. / rr;
                rri3 = Cube(rri);
                fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                VVSAdd (mol[j1].ra, fcVal, dr);
                VVSAdd (mol[j2].ra, - fcVal, dr);
                uSum += 4. * rri3 * (rri3 - 1.) + 1.;
                virSum += fcVal * rr;
            }
        }
    }

}

void LeapfrogStep (int64_t part){
    int n;

    if (part == 1) {
        DO_MOL {
            VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
            VVSAdd (mol[n].r, deltaT, mol[n].rv);
        }
    } else {
        DO_MOL VVSAdd (mol[n].rv, 0.5 * deltaT, mol[n].ra);
    }
}

void ApplyBoundaryCond() {
    int64_t n;

    DO_MOL VWrapAll (mol[n].r);
}

void SingleStep(){
    ++stepCount;
    timeNow = stepCount * deltaT;
    LeapfrogStep(1);
    ApplyBoundaryCond();
    ComputeForces();
    LeapfrogStep(2);
    EvalProps();
    AccumProps(1);
    if (stepCount % stepAvg == 0) {
        AccumProps(2);
        PrintSummary(stdout);
        AccumProps(0);
    }
};


int main() {

    

    GetNameList();
    PrintNameList(stdout);
    SetParams();
    SetupJob();
    moreCycles = 1;
    while (moreCycles) {
        SingleStep();
        if (stepCount >= stepLimit) moreCycles = 0;
    }
}