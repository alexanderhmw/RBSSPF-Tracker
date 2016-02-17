#ifndef TRACKER_CUH
#define TRACKER_CUH

#include<cuda.h>
#include<cuda_runtime.h>
#include<thrust/random/linear_congruential_engine.h>
#include<thrust/random/uniform_real_distribution.h>
#include<thrust/random/normal_distribution.h>
#include<thrust/generate.h>

#include<chrono>
#include<vector>
#include<list>
#include<random>
#include<time.h>
#include<algorithm>

#define PI 3.14159265359
#define DEG2RAD(ang) (ang*PI/180)

#define MINSIGMA 1e-2
#define UNCERTAINTHRESH 0.4
#define MAXSIGMA 1e6

#define RQPN 256
#define SPN 4
#define MAXPN (SPN*RQPN)
#define MAXBEAMNUM 2048

#define CUDAFREE(pointer) if(pointer!=NULL){cudaFree(pointer);pointer=NULL;}

#define THREAD_1D 256
#define THREAD_2D 16
#define GetKernelDim_1D(numBlocks, threadsPerBlock, dim) int numBlocks=(dim+THREAD_1D-1)/THREAD_1D; int threadsPerBlock=THREAD_1D;
#define GetKernelDim_2D(numBlocks, threadsPerBlock, xdim, ydim) dim3 numBlocks(int((xdim+THREAD_2D-1)/THREAD_2D), int((ydim+THREAD_2D-1)/THREAD_2D)); dim3 threadsPerBlock(THREAD_2D, THREAD_2D);
#define GetThreadID_1D(id) int id=blockDim.x*blockIdx.x+threadIdx.x;
#define GetThreadID_2D(xid,yid) int xid=blockDim.x*blockIdx.x+threadIdx.x;int yid=blockDim.y*blockIdx.y+threadIdx.y;

#define NEARESTRING 3.35
#define MINBEAM 2
#define MAXBEAM 100

//==========================

#define SIGMA 0.01

#define COST0 1
#define WEIGHT0 -2
#define COST1 2
#define WEIGHT1 -8
#define COST2 0
#define WEIGHT2 0
#define COST3 1.6
#define WEIGHT3 -5.12

#define MARGIN0 0.2
#define MARGIN1 0.1
#define MARGIN2 0.1

//==========================

#define SSPFFLAG 1
#define REJECTFLAG 0

//==========================
#define MOTIONMIN {DEG2RAD(-60),-10,-0.5,DEG2RAD(-90)}
#define MOTIONMAX {DEG2RAD(60),30,0.5,DEG2RAD(90)}
#define MOTIONPREC {DEG2RAD(1),1,0.001,DEG2RAD(1)}
#define INITMOTIONOFFSET {DEG2RAD(60),20,0.5,DEG2RAD(90)}
#define UPDATEMOTIONOFFSET_SSPF {DEG2RAD(30),10,0.15,DEG2RAD(30)}
#define UPDATEMOTIONOFFSET_PF {DEG2RAD(10),3,0.05,DEG2RAD(10)}

#define GEOMETRYMIN {DEG2RAD(-30),0,0,0,0}
#define GEOMETRYMAX {DEG2RAD(30),3,3,5,5}
#define GEOMETRYPREC {DEG2RAD(1),0.1,0.1,0.1,0.1}
#define INITGEOMETRYOFFSET {DEG2RAD(30),1.5,1.5,2.5,2.5}
#define UPDATEGEOMETRYOFFSET {DEG2RAD(1),1.5,1.5,2.5,2.5}

struct GeometrySampleParam
{
    double theta;
    double wl,wr,lf,lb;
};

struct MotionSampleParam
{
    double a,v,k,omega;
};

struct TrackerSampleControl
{
    int id;
    bool pfflag;

    MotionSampleParam motionmin;
    MotionSampleParam motionmax;
    MotionSampleParam motionprec;
    MotionSampleParam motionoffset;
    MotionSampleParam motionzoom;

    double motioniteration;
    double motionanneal;
    double motionannealratio;

    GeometrySampleParam geometrymin;
    GeometrySampleParam geometrymax;
    GeometrySampleParam geometryprec;
    GeometrySampleParam geometryoffset;
    GeometrySampleParam geometryzoom;

    double geometryiteration;
    double geometryanneal;
    double geometryannealratio;
};

//====================================================

enum TrackerStatus
{
    StatusInitGeometry,
    StatusInitMotion,
    StatusUpdateTracker_SSPF,
    StatusUpdateTracker_PF,
    StatusFinishTracker
};

struct TrackerState
{
    double x,y,theta;
    double wl,wr,lf,lb;
    double a,v,k,omega;
};

struct TrackerGeometry
{
    double cx[4],cy[4];
    double cn[4],sa[4];
    double dx[4],dy[4];
};

struct TrackerResult //CPU
{
    int id;
    TrackerStatus status;
    TrackerState mean;
    TrackerState sigma;
};

struct TrackerParticle //GPU-core
{
    double weight;
    int count;
    TrackerState state;
    TrackerGeometry geometry;
};

//====================================================

struct LaserScan //CPU
{
    double timestamp;
    double x,y,theta;
    int beamnum;
    double length[MAXBEAMNUM];
};

struct EgoMotion //GPU-cache
{
    bool validflag=0;
    double x,y,theta;
    double timestamp;
    double dx=0,dy=0,dtheta=0;
    double dt=0;
};

extern "C" void cuda_openTracker();
extern "C" void cuda_closeTracker();
extern "C" void cuda_SetLaserScan(LaserScan & scan);
extern "C" void cuda_UpdateTracker(int trackernum, TrackerResult *trackers);

#endif // TRACKER_CUH

