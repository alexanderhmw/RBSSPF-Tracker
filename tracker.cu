#include"tracker.cuh"

using namespace std::chrono;

LaserScan h_scan;
EgoMotion h_egomotion;

int h_seed[MAXPN];
thrust::minstd_rand * d_rng=NULL;

//==============================================================================

__global__
void kernelSetRandomSeed(int * seed, thrust::minstd_rand * rng, int tmppnum)
{
    GetThreadID_1D(id);
    if(id>=tmppnum) return;
    rng[id]=thrust::minstd_rand(seed[id]);
    return;
}

extern "C" void cudaopenTracker()
{
    cudacloseTracker();
    //==============================
    //initialize rand seed
    int * d_seed;
    cudaMalloc(&d_seed,sizeof(int)*MAXPN);
    thrust::generate(h_seed,h_seed+MAXPN,rand);
    cudaMemcpy(d_seed,h_seed,sizeof(int)*MAXPN,cudaMemcpyHostToDevice);
    cudaMalloc(&d_rng,sizeof(thrust::minstd_rand)*MAXPN);
    GetKernelDim_1D(blocks,threads,MAXPN);
    kernelSetRandomSeed<<<blocks,threads>>>(d_seed,d_rng,MAXPN);
    CUDAFREE(d_seed);
}

extern "C" void cudacloseTracker()
{
    CUDAFREE(d_rng);
}

//==============================================================================

extern "C" void cudaSetLaserScan(LaserScan & scan)
{
    h_scan=scan;
    if(h_egomotion.validflag)
    {
        double tmpdx=h_egomotion.x-h_scan.x;
        double tmpdy=h_egomotion.y-h_scan.y;
        double c=cos(h_scan.theta);
        double s=sin(h_scan.theta);
        h_egomotion.dx=c*tmpdx+s*tmpdy;
        h_egomotion.dy=-s*tmpdx+c*tmpdy;
        h_egomotion.dtheta=h_egomotion.theta-h_scan.theta;
        h_egomotion.dt=h_scan.timestamp-h_egomotion.timestamp;
    }
    h_egomotion.x=h_scan.x;
    h_egomotion.y=h_scan.y;
    h_egomotion.theta=h_scan.theta;
    h_egomotion.timestamp=h_scan.timestamp;
    h_egomotion.validflag=1;
}

//==============================================================================

__global__
void kernelInitParticle(TrackerState value, TrackerParticle * particles, int pnum)
{
    GetThreadID_1D(pid);
    if(pid>=pnum) return;
    particles[pid].state=value;
}

//==============================================================================

__host__ __device__
void deviceEgoMotion(TrackerParticle & particle, EgoMotion & egomotion)
{
    double c=cos(egomotion.dtheta);
    double s=sin(egomotion.dtheta);
    double tmpx=c*particle.state.x-s*particle.state.y+egomotion.dx;
    double tmpy=s*particle.state.x+c*particle.state.y+egomotion.dy;
    particle.state.x=tmpx;
    particle.state.y=tmpy;
    particle.state.theta+=egomotion.dtheta;
}

__host__ __device__
void deviceAckermannModel(TrackerParticle & particle, EgoMotion & egomotion)
{
    double c=cos(particle.state.theta);
    double s=sin(particle.state.theta);

    if(particle.state.k==0)
    {
        particle.state.x=particle.state.x+c*particle.state.v*egomotion.dt;
        particle.state.y=particle.state.y+s*particle.state.v*egomotion.dt;
        particle.state.a=0;
    }
    else
    {
        double c0=cos(particle.state.theta+particle.state.a);
        double s0=sin(particle.state.theta+particle.state.a);

        particle.state.omega=particle.state.v*particle.state.k;
        double dtheta=particle.state.omega*egomotion.dt;
        particle.state.theta+=dtheta;

        double c1=cos(particle.state.theta+particle.state.a);
        double s1=sin(particle.state.theta+particle.state.a);
        double R=1/particle.state.k;

        particle.state.x=particle.state.x+R*(-s0+s1);
        particle.state.y=particle.state.y+R*(c0-c1);
    }
}

__host__ __device__
void deviceBuildModel(TrackerParticle & particle)
{
    double c=cos(particle.state.theta);
    double s=sin(paritcle.state.theta);

    particle.geometry.cx[0]=c*particle.state.lf-s*particle.state.wl+particle.state.x; particle.geometry.cy[0]=s*particle.state.lf+c*particle.state.wl+particle.state.y;
    particle.geometry.cx[1]=c*particle.state.lf+s*particle.state.wr+particle.state.x; particle.geometry.cy[1]=s*particle.state.lf-c*particle.state.wr+particle.state.y;
    particle.geometry.cx[2]=-c*particle.state.lb+s*particle.state.wr+particle.state.x; particle.geometry.cy[2]=-s*particle.state.lb-c*particle.state.wr+particle.state.y;
    particle.geometry.cx[3]=-c*particle.state.lb-s*particle.state.wl+particle.state.x; particle.geometry.cy[3]=-s*particle.state.lb+c*particle.state.wl+particle.state.y;

    double width=particle.state.wl+particle.state.wr;
    double length=particle.state.lf+particle.state.lb;
    particle.geometry.dx[0]=(particle.geometry.cx[1]-particle.geometry.cx[0])/width; particle.geometry.dy[0]=(particle.geometry.cy[1]-particle.geometry.cy[0])/width;
    particle.geometry.dx[1]=(particle.geometry.cx[2]-particle.geometry.cx[1]/length); particle.geometry.dy[1]=(particle.geometry.cy[2]-particle.geometry.cy[1])/length;
    particle.geometry.dx[2]=(particle.geometry.cx[3]-particle.geometry.cx[2])/width; particle.geometry.dy[2]=(particle.geometry.cy[3]-particle.geometry.cy[2])/width;
    particle.geometry.dx[3]=(particle.geometry.cx[0]-particle.geometry.cx[3])/length; particle.geometry.dy[3]=(particle.geometry.cy[0]-particle.geometry.cy[3])/length;

    for(int i=0;i<4;i++)
    {
        particle.geometry.cn[i]=sqrt(particle.geometry.cx[i]*particle.geometry.cx[i]+particle.geometry.cy[i]*particle.geometry.cy[i]);
        particle.geometry.sa[i]=(particle.geometry.cx[i]*particle.geometry.dy[i]-particle.geometry.cy[i]*particle.geometry.dx[i])/particle.geometry.cn[i];
    }
}

//==============================================================================

__global__
void kernelMotionUpSample(TrackerParticle * particles, TrackerSampleControl * controls, TrackerParticle * tmpparticles, int tmppnum, thrust::minstd_rand * rng, EgoMotion egomotion)
{
    GetThreadID_1D(tmppid);
    if(tmppid>=tmppnum) return;
    int pid=int(tmppid/SPN);
    int vid=int(pid/RQPN);
    int rid=tmppid%MAXPN;

    TrackerSampleControl control=controls[vid];
    if(control.motioniteration<0) return;

    TrackerParticle particle;
    particle.weight=0;
    particle.count=0;
    particle.state=particles[pid].state;

    if(control.pfflag)
    {
        particle.state.v=thrust::random::normal_distribution<double>(particle.state.v,control.motionoffset.v)(rng[rid]);
        particle.state.v=particle.state.v>control.motionmin.v?particle.state.v:control.motionmin.v;
        particle.state.v=particle.state.v<control.motionmax.v?particle.state.v:control.motionmax.v;

        particle.state.omega=thrust::random::normal_distribution<double>(particle.state.omega,control.motionoffset.omega)(rng[rid]);
        particle.state.omega=particle.state.omega>control.motionmin.omega?particle.state.omega:control.motionmin.omega;
        particle.state.omega=particle.state.omega<control.motionmax.omega?particle.state.omega:control.motionmax.omega;

        particle.state.a=thrust::random::normal_distribution<double>(particle.state.a,control.motionoffset.a)(rng[rid]);
        particle.state.a=particle.state.a>control.motionmin.a?particle.state.a:control.motionmin.a;
        particle.state.a=particle.state.a<control.motionmax.a?particle.state.a:control.motionmax.a;
    }
    else
    {
        double vmin=particle.state.v-control.motionoffset.v; vmin=vmin>control.motionmin.v?vmin:control.motionmin.v;
        double vmax=particle.state.v+control.motionoffset.v; vmax=vmax<control.motionmax.v?vmax:control.motionmax.v;
        particle.state.v=thrust::random::uniform_real_distribution<double>(vmin,vmax)(rng[rid]);

        double omegamin=particle.state.omega-control.motionoffset.omega; omegamin=omegamin>control.motionmin.omega?omegamin:control.motionmin.omega;
        double omegamax=particle.state.omega+control.motionoffset.omega; omegamax=omegamax<control.motionmax.omega?omegamax:control.motionmax.omega;
        particle.state.omega=thrust::random::uniform_real_distribution<double>(omegamin,omegamax)(rng[rid]);

        double amin=particle.state.a-control.motionoffset.a; amin=amin>control.motionmin.a?amin:control.motionmin.a;
        double amax=particle.state.a+control.motionoffset.a; amax=amax<control.motionmax.a?amax:control.motionmax.a;
        particle.state.a=thrust::random::uniform_real_distribution<double>(amin,amax)(rng[rid]);
    }

    if(particle.state.v==0)
    {
        particle.state.k=0;
        particle.state.omega=0;
        particle.state.a=0;
    }
    else
    {
        particle.state.k=particle.state.omega/particle.state.v;
        particle.state.k=particle.state.k>control.motionmin.k?particle.state.k:control.motionmin.k;
        particle.state.k=particle.state.k<control.motionmax.k?particle.state.k:control.motionmax.k;
        particle.state.omega=particle.state.v*particle.state.k;
        deviceAckermannModel(particle,egomotion);
    }
    deviceEgoMotion(particle,egomotion);
    deviceBuildModel(particle);

    tmpparticles[tmppid]=particle;
}

__global__
void kernelGeometryUpSample(TrackerParticle * particles, TrackerSampleControl * controls, TrackerParticle * tmpparticles, int tmppnum, thrust::minstd_rand * rng)
{
    GetThreadID_1D(tmppid);
    if(tmppid>=tmppnum) return;
    int pid=int(tmppid/SPN);
    int vid=int(pid/RQPN);
    int rid=tmppid%MAXPN;

    TrackerSampleControl control=controls[vid];
    if(control.geometryiteration<0) return;

    TrackerParticle particle;
    particle.weight=0;
    particle.count=0;
    particle.state=particles[pid].state;

    if(control.geometryoffset.theta>control.geometryprec.theta)
    {
        double thetamin=particle.state.theta-control.geometryoffset.theta; thetamin=thetamin>control.geometrymin.theta?thetamin:control.geometrymin.theta;
        double thetamax=particle.state.theta+control.geometryoffset.theta; thetamax=thetamax<control.geometrymax.theta?thetamax:control.geometrymax.theta;
        particle.state.theta=thrust::random::uniform_real_distribution<double>(thetamin,thetamax)(rng[rid]);
    }

    double wlmin=particle.state.wl-control.geometryoffset.wl; wlmin=wlmin>control.geometrymin.wl?wlmin:control.geometrymin.wl;
    double wlmax=particle.state.wl+control.geometryoffset.wl; wlmax=wlmax<control.geometrymax.wl?wlmax:control.geometrymax.wl;
    particle.state.wl=thrust::random::uniform_real_distribution<double>(wlmin,wlmax)(rng[rid]);

    double wrmin=particle.state.wr-control.geometryoffset.wr; wrmin=wrmin>control.geometrymin.wr?wrmin:control.geometrymin.wr;
    double wrmax=particle.state.wr+control.geometryoffset.wr; wrmax=wrmax<control.geometrymax.wr?wrmax:control.geometrymax.wr;
    particle.state.wr=thrust::random::uniform_real_distribution<double>(wrmin,wrmax)(rng[rid]);

    double lfmin=particle.state.lf-control.geometryoffset.lf; lfmin=lfmin>control.geometrymin.lf?lfmin:control.geometrymin.lf;
    double lfmax=particle.state.lf+control.geometryoffset.lf; lfmax=lfmax<control.geometrymax.lf?lfmax:control.geometrymax.lf;
    particle.state.lf=thrust::random::uniform_real_distribution<double>(lfmin,lfmax)(rng[rid]);

    double lbmin=particle.state.lb-control.geometryoffset.lb; lbmin=lbmin>control.geometrymin.lb?lbmin:control.geometrymin.lb;
    double lbmax=particle.state.lb+control.geometryoffset.lb; lbmax=lbmax<control.geometrymax.lb?lbmax:control.geometrymax.lb;
    particle.state.lb=thrust::random::uniform_real_distribution<double>(lbmin,lbmax)(rng[rid]);

    deviceBuildModel(particle);

    tmpparticles[tmppid]=particle;
}

//==============================================================================

__global__
void kernelMeasureScan(TrackerParticle * tmpparticles, int tmppnum, TrackerSampleControl * controls, double * weight, double bear, double length, bool motionflag)
{
    GetThreadID_1D(tmppid);
    if(tmppid>=tmppnum) return;
    int vid=int(tmppid/MAXPN);

    double iteration=motionflag?controls[vid].motioniteration:controls[vid].geometryiteration;
    if(iteration<0) return;

    double anneal=motionflag?controls[vid].motionanneal:controls[vid].geometryanneal;
    TrackerParticle particle=tmpparticles[tmppid];

    double lx=cos(bear);
    double ly=sin(bear);
    for(int i=0;i<4;i++)
    {
        int j=(i+1)%4;
        double s0=particle.geometry.cx[i]*ly-lx*particle.geometry.cy[i];
        double s1=lx*particle.geometry.cy[j]-particle.geometry.cx[j]*ly;
        double cn=lx*particle.geometry.dx[j]+ly*particle.geometry.dy[j];

        bool flag=s0>0&&s1>0&&cn>0;

        if(flag)
        {
            double sa=particle.geometry.sa[i];
            double sb=lx*particle.geometry.dy[i]-particle.geometry.dx[i]*ly;
            double l=sa/sb*particle.geometry.cn[i];

            double l0=l-MARGIN0/cn;
            double l1=l-MARGIN1/cn;
            double l2=l;
            double l3=l+MARGIN2/cn;

            double tmplogweight;
            double delta,w1,w2;
            if(length<=l0)
            {
                delta=length-l0;
                w1=WEIGHT0-WEIGHT0;
                w2=WEIGHT1-WEIGHT0;
            }
            else if(lenght<=l1)
            {
                delta=length-l1;
                w1=WEIGHT1-WEIGHT0;
                w2=WEIGHT2-WEIGHT0;
            }
            else if(lenght<=l3)
            {
                delta=length-l2;
                w1=WEIGHT2-WEIGHT0;
                w2=2*w1;
                particle.count++;
            }
            else
            {
                delta=length-l3;
                w1=WEIGHT3-WEIGHT0;
                w2=WEIGHT2-WEIGHT0;
            }
            particle.weight+=(w1+(w2-w1)*exp(-delta*delta/SIGMA))/anneal;
        }
    }
    weight[tmppid]=particle.weight;
    tmpparticles[tmppid]=particle;
}

__global__
void kernelDownSample(TrackerParticle * particles, int pnum, int * sampleids, TrackerParticle * tmpparticles, bool poseupdateflag)
{
    GetThreadID_1D(pid);
    if(pid>=pnum) return;
    if(sampleids[pid]<0) return;

    if(poseupdateflag)
    {
        particles[pid]=tmpparticles[sampleids[pid]];
    }
    else
    {
        particles[pid].state.a=tmpparticles[sampleids[pid]].state.a;
        particles[pid].state.v=tmpparticles[sampleids[pid]].state.v;
        particles[pid].state.k=tmpparticles[sampleids[pid]].state.k;
        particles[pid].state.omega=tmpparticles[sampleids[pid]].state.omega;
    }
}

//==============================================================================

#define CALRATIO(ratio, vratio, maxratio, maxrange, minrange) \
    ratio=maxrange/minrange; vratio*=ratio; maxratio=ratio>maxratio?ratio:maxratio;
#define CALZOOM(zoom, maxrange, minrange, N) \
    zoom=log(maxrange/minrange)/log(2)/N;zoom=1/pow(2,zoom);

void hostCalculateMotionControl(MotionSampleParam & offset, MotionSampleParam & prec, double & iteration, double & anneal, double & annealratio, MotionSampleParam & zoom)
{
    double ratio=1,vratio=1,maxratio=1;

    CALRATIO(ratio,vratio,maxratio,offset.a,prec.a);
    CALRATIO(ratio,vratio,maxratio,offset.v,prec.v);
    CALRATIO(ratio,vratio,maxratio,offset.omega,prec.omega);

    iteration=log(vratio)/log(2);
    anneal=maxratio*maxratio;
    annealratio=pow(anneal,-1/iteration);

    CALZOOM(zoom.a,offset.a,prec.a,iteration);
    CALZOOM(zoom.v,offset.v,prec.v,iteration);
    CALZOOM(zoom.omega,offset.omega,prec.omega,iteration);
}

void hostSetupGeometryOffset(TrackerState & sigma, GeometrySampleParam & offset, GeometrySampleParam & prec)
{
    offset=UPDATEGEOMETRYOFFSET;
    if(sigma.wl<prec.wl) offset.wl=prec.wl;
    if(sigma.wr<prec.wr) offset.wr=prec.wr;
    if(sigma.lf<prec.lf) offset.lf=prec.lf;
    if(sigma.lb<prec.lb) offset.lb=prec.lb;
}

void hostCalculateGeometryControl(GeometrySampleParam & offset, GeometrySampleParam & prec, double & iteration, double & anneal, double & annealratio, GeometrySampleParam & zoom)
{
    double ratio=1,vratio=1,maxratio=1;

    CALRATIO(ratio,vratio,maxratio,offset.theta,prec.theta);
    CALRATIO(ratio,vratio,maxratio,offset.wl,prec.wl);
    CALRATIO(ratio,vratio,maxratio,offset.wr,prec.wr);
    CALRATIO(ratio,vratio,maxratio,offset.lf,prec.lf);
    CALRATIO(ratio,vratio,maxratio,offset.lb,prec.lb);

    iteration=log(vratio)/log(2);
    anneal=maxratio*maxratio;
    annealratio=pow(anneal,-1/iteration);

    CALZOOM(zoom.theta,offset.theta,prec.theta,iteration);
    CALZOOM(zoom.wl,offset.wl,prec.wl,iteration);
    CALZOOM(zoom.wr,offset.wr,prec.wr,iteration);
    CALZOOM(zoom.lf,offset.lf,prec.lf,iteration);
    CALZOOM(zoom.lb,offset.lb,prec.lb,iteration);
}

void hostInitializeControl(TrackerResult & tracker, TrackerSampleControl & control)
{
    control.motionmin=MOTIONMIN;
    control.motionmax=MOTIONMAX;
    control.motionprec=MOTIONPREC;

    control.geometrymin=GEOMETRYMIN;
    control.geometrymax=GEOMETRYMAX;
    control.geometryprec=GEOMETRYPREC;

    switch(tracker.status)
    {
    case StatusInitGeometry:
        {
            control.motioniteration=-1;
            control.geometryoffset=INITGEOMETRYOFFSET;
            hostCalculateGeometryControl(control.geometryoffset,control.geometryprec,control.geometryiteration,control.geometryanneal,control.geometryannealratio,control.geometryzoom);
        }
        break;
    case StatusInitMotion:
        {
            control.pfflag=0;
            control.motionoffset=INITMOTIONOFFSET;
            hostCalculateMotionControl(control.motionoffset,control.motionprec,control.motioniteration,control.motionanneal,control.motionannealratio,control.motionzoom);
            hostSetupGeometryOffset(tracker.sigma,control.geometryoffset,control.geometryprec);
            hostCalculateGeometryControl(control.geometryoffset,control.geometryprec,control.geometryiteration,control.geometryanneal,control.geometryannealratio,control.geometryzoom);
        }
        break;
    case StatusUpdateTracker_SSPF:
        {
            control.pfflag=0;
            control.motionoffset=UPDATEMOTIONOFFSET_SSPF;
            hostCalculateMotionControl(control.motionoffset,control.motionprec,control.motioniteration,control.motionanneal,control.motionannealratio,control.motionzoom);
            hostSetupGeometryOffset(tracker.sigma,control.geometryoffset,control.geometryprec);
            hostCalculateGeometryControl(control.geometryoffset,control.geometryprec,control.geometryiteration,control.geometryanneal,control.geometryannealratio,control.geometryzoom);
        }
        break;
    case StatusUpdateTracker_PF:
        {
            control.pfflag=1;
            control.motionoffset=UPDATEMOTIONOFFSET_PF;
            control.motioniteration=0;
            hostSetupGeometryOffset(tracker.sigma,control.geometryoffset,control.geometryprec);
            hostCalculateGeometryControl(control.geometryoffset,control.geometryprec,control.geometryiteration,control.geometryanneal,control.geometryannealratio,control.geometryzoom);
        }
    default:
        control.motioniteration=-1;
        control.geometryiteration=-1;
        break;
    }
    control.geometrymin.theta=tracker.mean.theta-control.geometryoffset.theta;
    control.geometrymax.theta=tracker.mean.theta+control.geometryoffset.theta;
}

//==============================================================================

void hostDownSampleID(int startid, double * weight, int * sampleids, double iteration)
{
    if(iteration<0)
    {
        for(int i=0;i<RQPN;i++)
        {
            sampleids[i]=-1;
        }
        return;
    }

    double maxlogweight=weight[0];
    double minlogweight=weight[0];
    for(int i=0;i<MAXPN;i++)
    {
        maxlogweight=maxlogweight>weight[i]?maxlogweight:weight[i];
        minlogweight=minlogweight<weight[i]?minlogweight:weight[i];
    }
    double maxscale=maxlogweight<30?1:30/maxlogweight;
    double minscale=minlogweight>-30?1:-30/minlogweight;
    double scale=maxscale<minscale?maxscale:minscale;

    weight[0]=exp(weight[0]*scale);
    for(int i=1;i<MAXPN;i++)
    {
        weight[i]=exp(weight[i]*scale);
        weight[i]+=weight[i-1];
    }

    double step=1.0/RQPN;
    int accuracy=1e6;
    double samplebase=(rand()%accuracy)*step/accuracy;
    double weightsum=weight[MAXPN-1];

    for(int i=0,j=0;i<RQPN;i++)
    {
        double sample=samplebase+i*step;
        while(j<MAXPN)
        {
            if(sample>weight[j]/weightsum)
            {
                j++;
                continue;
            }
            else
            {
                sampleids[i]=startid+j;
                break;
            }
        }
    }
}

//==============================================================================

void hostEstimateResult(TrackerParticle * particles, TrackerResult & result, bool geometryflag)
{
    if(geometryflag)
    {
        result.mean.wl=0;
        result.mean.wr=0;
        result.mean.lf=0;
        result.mean.lb=0;
    }
    else
    {
        result.mean.x=0;
        result.mean.y=0;
        result.mean.a=0;
        result.mean.v=0;
        result.mean.k=0;
        result.mean.omega=0;
    }
    result.mean.theta=0;

    TrackerState minstate=particles[0].state;
    TrackerState maxstate=particles[0].state;

    for(int i=0;i<RQPN;i++)
    {
        if(geometryflag)
        {
            result.mean.wl+=particles[i].state.wl;
            result.mean.wr+=particles[i].state.wr;
            result.mean.lf+=particles[i].state.lf;
            result.mean.lb+=particles[i].state.lb;

            minstate.wl=minstate.wl<particles[i].state.wl?minstate.wl:particles[i].state.wl;
            maxstate.wl=maxstate.wl>particles[i].state.wl?maxstate.wl:particles[i].state.wl;
            minstate.wr=minstate.wr<particles[i].state.wr?minstate.wr:particles[i].state.wr;
            maxstate.wr=maxstate.wr>particles[i].state.wr?maxstate.wr:particles[i].state.wr;
            minstate.lf=minstate.lf<particles[i].state.lf?minstate.lf:particles[i].state.lf;
            maxstate.lf=maxstate.lf>particles[i].state.lf?maxstate.lf:particles[i].state.lf;
            minstate.lb=minstate.lb<particles[i].state.lb?minstate.lb:particles[i].state.lb;
            maxstate.lb=maxstate.lb>particles[i].state.lb?maxstate.lb:particles[i].state.lb;
        }
        else
        {
            result.mean.x+=particles[i].state.x;
            result.mean.y+=particles[i].state.y;
            result.mean.a+=particles[i].state.a;
            result.mean.v+=particles[i].state.v;
            result.mean.k+=particles[i].state.k;
            result.mean.omega+=particles[i].state.omega;

            minstate.x=minstate.x<particles[i].state.x?minstate.x:particles[i].state.x;
            maxstate.x=maxstate.x>particles[i].state.x?maxstate.x:particles[i].state.x;
            minstate.y=minstate.y<particles[i].state.y?minstate.y:particles[i].state.y;
            maxstate.y=maxstate.y>particles[i].state.y?maxstate.y:particles[i].state.y;
            minstate.a=minstate.a<particles[i].state.a?minstate.a:particles[i].state.a;
            maxstate.a=maxstate.a>particles[i].state.a?maxstate.a:particles[i].state.a;
            minstate.v=minstate.v<particles[i].state.v?minstate.v:particles[i].state.v;
            maxstate.v=maxstate.v>particles[i].state.v?maxstate.v:particles[i].state.v;
            minstate.k=minstate.k<particles[i].state.k?minstate.k:particles[i].state.k;
            maxstate.k=maxstate.k>particles[i].state.k?maxstate.k:particles[i].state.k;
            minstate.omega=minstate.omega<particles[i].state.omega?minstate.omega:particles[i].state.omega;
            maxstate.omega=maxstate.omega>particles[i].state.omega?maxstate.omega:particles[i].state.omega;
        }
        result.mean.theta+=particles[i].state.theta;
        minstate.theta=minstate.theta<particles[i].state.theta?minstate.theta:particles[i].state.theta;
        maxstate.theta=maxstate.theta>particles[i].state.theta?maxstate.theta:particles[i].state.theta;
    }

    if(geometryflag)
    {
        result.mean.wl/=RQPN;
        result.mean.wr/=RQPN;
        result.mean.lf/=RQPN;
        result.mean.lb/=RQPN;

        result.sigma.wl=std::max(result.mean.wl-minstate.wl,maxstate.wl-result.mean.wl);
        result.sigma.wr=std::max(result.mean.wr-minstate.wr,maxstate.wr-result.mean.wr);
        result.sigma.lf=std::max(result.mean.lf-minstate.lf,maxstate.lf-result.mean.lf);
        result.sigma.lb=std::max(result.mean.lb-minstate.lb,maxstate.lb-result.mean.lb);

        result.sigma.wl=result.sigma.wl>MINSIGMA?result.sigma.wl:MINSIGMA;
        result.sigma.wl=result.sigma.wl<UNCERTAINTHRESH?result.sigma.wl:MAXSIGMA;
        result.sigma.wr=result.sigma.wr>MINSIGMA?result.sigma.wr:MINSIGMA;
        result.sigma.wr=result.sigma.wr<UNCERTAINTHRESH?result.sigma.wr:MAXSIGMA;
        result.sigma.lf=result.sigma.lf>MINSIGMA?result.sigma.lf:MINSIGMA;
        result.sigma.lf=result.sigma.lf<UNCERTAINTHRESH?result.sigma.lf:MAXSIGMA;
        result.sigma.lb=result.sigma.lb>MINSIGMA?result.sigma.lb:MINSIGMA;
        result.sigma.lb=result.sigma.lb<UNCERTAINTHRESH?result.sigma.lb:MAXSIGMA;
    }
    else
    {
        result.mean.x/=RQPN;
        result.mean.y/=RQPN;
        result.mean.a/=RQPN;
        result.mean.v/=RQPN;
        result.mean.k/=RQPN;
        result.mean.omega/=RQPN;

        result.sigma.x=std::max(result.mean.x-minstate.x,maxstate.x-result.mean.x);
        result.sigma.y=std::max(result.mean.y-minstate.y,maxstate.y-result.mean.y);
        result.sigma.a=std::max(result.mean.a-minstate.a,maxstate.a-result.mean.a);
        result.sigma.v=std::max(result.mean.v-minstate.v,maxstate.v-result.mean.v);
        result.sigma.k=std::max(result.mean.k-minstate.k,maxstate.k-result.mean.k);
        result.sigma.omega=std::max(result.mean.omega-minstate.omega,maxstate.omega-result.mean.omega);
    }
    result.mean.theta/=RQPN;
    result.sigma.theta=std::max(result.mean.theta-minstate.theta,maxstate.theta-result.mean.theta);
}

//==============================================================================

extern "C" void cudaUpdateTracker(int trackernum, TrackerResult * trackers)
{
    if(trackernum<=0) return;
    double density=2*PI/h_scan.beamnum;
    GetKernelDim_1D(maxblocks,maxthreads,MAXPN*trackernum);
    GetKernelDim_1D(rqblocks,rqthreads,RQPN*trackernum);

    //==============================
    //allocate memory of weight and sampleids
    double * h_weight=new double[MAXPN*trackernum];
    double * d_weight=NULL;
    cudaMalloc(&d_weight,sizeof(double)*MAXPN*trackernum);

    int * h_sampleids=new int[RQPN*trackernum];
    int * d_sampleids=NULL;
    cudaMalloc(&d_sampleids,sizeof(int)*RQPN*trackernum);

    //==============================
    //initialize sample control
    TrackerSampleControl * h_controls=new TrackerSampleControl[trackernum];
    TrackerSampleControl * d_controls=NULL;
    cudaMalloc(&d_controls,sizeof(TrackerSampleControl)*trackernum);

    double maxmotioniteration=-1;
    double maxgeometryiteration=-1;
    for(int i=0;i<trackernum;i++)
    {
        h_controls[i].id=i;
        hostInitializeControl(trackers[i],h_controls[i]);
        if(maxmotioniteration<h_controls[i].motioniteration) maxmotioniteration=h_controls[i].motioniteration;
    }

    //==============================
    //allocate memory of paritlces and tmpparticles
    TrackerParticle * h_particles=new TrackerParticle[RQPN*trackernum];
    TrackerParticle * d_particles;
    cudaMalloc(&d_particles,sizeof(TrackerParticle)*RQPN*trackernum);
    TrackerParticle * d_tmpparticles;
    cudaMalloc(&d_tmpparticles,sizeof(TrackerParticle)*MAXPN*trackernum);

    //==============================
    //initialize particles
    {
        GetKernelDim_1D(blocks,threads,RQPN);
        TrackerParticle * d_particles_iter=d_particles;
        for(int i=0;i<trackernum;i++)
        {
            kernelInitParticle<<<blocks,threads>>>(trackers[i].mean,d_particles_iter,RQPN);
            d_particles_iter=d_particles_iter+RQPN;
        }
    }

    //==============================
    //Motion estimate
    //SSPF-loop
    for(int i=1;i<=maxmotioniteration;i++)
    {
        //upsample
        cudaMemcpy(d_controls,h_controls,sizeof(TrackerSampleControl)*trackernum,cudaMemcpyHostToDevice);
        kernelMotionUpSample<<<maxblocks,maxthreads>>>(d_particles,d_controls,d_tmpparticles,MAXPN*trackernum,d_rng,h_egomotion);

        //measure
        for(int j=0;j<h_scan.beamnum;j++)
        {
            double bear=density*j-PI;
            kernelMeasureScan<<<maxblocks,maxthreads>>>(d_tmpparticles,MAXPN*trackernum,d_controls,d_weight,bear,h_scan.length[j],1);
        }

        //downsample
        cudaMemcpy(h_weight,d_weight,sizeof(double)*MAXPN*trackernum,cudaMemcpyDeviceToHost);
        for(int j=0;j<trackernum;j++)
        {
            hostDownSampleID(j*MAXPN,h_weight+j*MAXPN,h_sampleids+j*RQPN,h_controls[j].motioniteration);
        }
        cudaMemcpy(d_sampleids,h_sampleids,sizeof(int)*RQPN*trackernum,cudaMemcpyHostToDevice);
        kernelDownSample<<<rqblocks,rqthreads>>>(d_particles,RQPN*trackernum,d_sampleids,d_tmpparticles,0);

        //update controls
        for(int j=0;j<trackernum;j++)
        {
            if(h_controls[j].motioniteration>=1)
            {
                h_controls[j].motioniteration--;
                h_controls[j].motionoffset.a*=h_controls[j].motionzoom.a;
                h_controls[j].motionoffset.v*=h_controls[j].motionzoom.v;
                h_controls[j].motionoffset.omega*=h_controls[j].motionzoom.omega;
                h_controls[j].motionanneal*=h_controls[j].motionannealratio;
            }
        }
    }
    //Final estimate
    if(maxmotioniteration>=0)
    {
        //setup controls
        for(int j=0;j<trackernum;j++)
        {
            if(h_controls[j].motioniteration>=0)
            {
                h_controls[j].motioniteration=0;
                h_controls[j].motionoffset=MOTIONPREC;
                h_controls[j].motionanneal=1;
            }
        }

        //upsample
        cudaMemcpy(d_controls,h_controls,sizeof(TrackerSampleControl)*trackernum,cudaMemcpyHostToDevice);
        kernelMotionUpSample<<<maxblocks,maxthreads>>>(d_particles,d_controls,d_tmpparticles,MAXPN*trackernum,d_rng,h_egomotion);

        //measure
        for(int j=0;j<h_scan.beamnum;j++)
        {
            double bear=density*j-PI;
            kernelMeasureScan<<<maxblocks,maxthreads>>>(d_tmpparticles,MAXPN*trackernum,d_controls,d_weight,bear,h_scan.length[j],1);
        }

        //downsample
        cudaMemcpy(h_weight,d_weight,sizeof(double)*MAXPN*trackernum,cudaMemcpyDeviceToHost);
        for(int j=0;j<trackernum;j++)
        {
            hostDownSampleID(j*MAXPN,h_weight+j*MAXPN,h_sampleids+j*RQPN,h_controls[j].motioniteration);
        }
        cudaMemcpy(d_sampleids,h_sampleids,sizeof(int)*RQPN*trackernum,cudaMemcpyHostToDevice);
        kernelDownSample<<<rqblocks,rqthreads>>>(d_particles,RQPN*trackernum,d_sampleids,d_tmpparticles,1);

        //Motion result
        cudaMemcpy(h_particles,d_particles,sizeof(TrackerParticle)*RQPN*trackernum,cudaMemcpyDeviceToHost);
        for(int j=0;j<trackernum;j++)
        {
            if(h_controls[j].motioniteration>=0)
            {
                hostEstimateResult(h_particles,trackers[j],0);
                if(trackers[j].sigma.x<0.5&&trackers[j].sigma.y<0.5&&trackers[j].sigma.theta<DEG2RAD(10))
                {
                    trackers[j].status=StatusUpdateTracker_SSPF;
                }
                else
                {
                    trackers[j].status=StatusUpdateTracker_PF;
                    h_controls[j].geometryiteration=-1;
                }
            }
            if(maxgeometryiteration<h_controls[j].geometryiteration) maxgeometryiteration=h_controls[j].geometryiteration;
        }
    }

    //==============================
    //initialize particles
    {
        GetKernelDim_1D(blocks,threads,RQPN);
        TrackerParticle * d_particles_iter=d_particles;
        for(int i=0;i<trackernum;i++)
        {
            kernelInitParticle<<<blocks,threads>>>(trackers[i].mean,d_particles_iter,RQPN);
            d_particles_iter=d_particles_iter+RQPN;
        }
    }

    //==============================
    //Geometry estimate
    //SSPF-loop
    for(int i=1;i<=maxgeometryiteration;i++)
    {
        //upsample
        cudaMemcpy(d_controls,h_controls,sizeof(TrackerSampleControl)*trackernum,cudaMemcpyHostToDevice);
        kernelGeometryUpSample<<<maxblocks,maxthreads>>>(d_particles,d_controls,d_tmpparticles,MAXPN*trackernum,d_rng);

        //measure
        for(int j=0;j<h_scan.beamnum;j++)
        {
            double bear=density*j-PI;
            kernelMeasureScan<<<maxblocks,maxthreads>>>(d_tmpparticles,MAXPN*trackernum,d_controls,d_weight,bear,h_scan.length[j],0);
        }

        //downsample
        cudaMemcpy(h_weight,d_weight,sizeof(double)*MAXPN*trackernum,cudaMemcpyDeviceToHost);
        for(int j=0;j<trackernum;j++)
        {
            hostDownSampleID(j*MAXPN,h_weight+j*MAXPN,h_sampleids+j*RQPN,h_controls[j].geometryiteration);
        }
        cudaMemcpy(d_sampleids,h_sampleids,sizeof(int)*RQPN*trackernum,cudaMemcpyHostToDevice);
        kernelDownSample<<<rqblocks,rqthreads>>>(d_particles,RQPN*trackernum,d_sampleids,d_tmpparticles,1);

        //update controls
        for(int j=0;j<trackernum;j++)
        {
            if(h_controls[j].geometryiteration>=1)
            {
                h_controls[j].geometryiteration--;
                h_controls[j].geometryoffset.theta*=h_controls[j].geometryzoom.theta;
                h_controls[j].geometryoffset.wl*=h_controls[j].geometryzoom.wl;
                h_controls[j].geometryoffset.wr*=h_controls[j].geometryzoom.wr;
                h_controls[j].geometryoffset.lf*=h_controls[j].geometryzoom.lf;
                h_controls[j].geometryoffset.lb*=h_controls[j].geometryzoom.lb;
                h_controls[j].geometryanneal*=h_controls[j].geometryannealratio;
            }
        }
    }
    //Final estimate
    if(maxgeometryiteration>=0)
    {
        //setup controls
        for(int j=0;j<trackernum;j++)
        {
            if(h_controls[j].geometryiteration>=0)
            {
                h_controls[j].geometryiteration=0;
                h_controls[j].geometryoffset=GEOMETRYPREC;
                h_controls[j].geometryanneal=1;
            }
        }

        //upsample
        cudaMemcpy(d_controls,h_controls,sizeof(TrackerSampleControl)*trackernum,cudaMemcpyHostToDevice);
        kernelGeometryUpSample<<<maxblocks,maxthreads>>>(d_particles,d_controls,d_tmpparticles,MAXPN*trackernum,d_rng);

        //measure
        for(int j=0;j<h_scan.beamnum;j++)
        {
            double bear=density*j-PI;
            kernelMeasureScan<<<maxblocks,maxthreads>>>(d_tmpparticles,MAXPN*trackernum,d_controls,d_weight,bear,h_scan.length[j],0);
        }

        //downsample
        cudaMemcpy(h_weight,d_weight,sizeof(double)*MAXPN*trackernum,cudaMemcpyDeviceToHost);
        for(int j=0;j<trackernum;j++)
        {
            hostDownSampleID(j*MAXPN,h_weight+j*MAXPN,h_sampleids+j*RQPN,h_controls[j].geometryiteration);
        }
        cudaMemcpy(d_sampleids,h_sampleids,sizeof(int)*RQPN*trackernum,cudaMemcpyHostToDevice);
        kernelDownSample<<<rqblocks,rqthreads>>>(d_particles,RQPN*trackernum,d_sampleids,d_tmpparticles,1);

        //Geometry result
        cudaMemcpy(h_particles,d_particles,sizeof(TrackerParticle)*RQPN*trackernum,cudaMemcpyDeviceToHost);
        for(int j=0;j<trackernum;j++)
        {
            if(h_controls[j].geometryiteration>=0)
            {
                if(trackers[j].status==StatusInitGeometry)
                {
                    hostEstimateResult(h_particles,trackers[j],1);
                    trackers[j].status=StatusInitMotion;
                }
                else
                {
                    TrackerResult pretracker=trackers[j];
                    hostEstimateResult(h_particles,trackers[j],1);
                    trackers[j].sigma.theta=pretracker.sigma.theta;

                    trackers[j].mean.wl=(pretracker.mean.wl*pow(trackers[j].sigma.wl,2)+trackers[j].mean.wl*pow(pretracker.sigma.wl,2))/(pow(pretracker.sigma.wl,2)+pow(trackers[j].sigma.wl,2));
                    trackers[j].sigma.wl=(pow(pretracker.sigma.wl,2)*pow(trackers[j].sigma.wl,2))/(pow(pretracker.sigma.wl,2)+pow(trackers[j].sigma.wl,2));
                    trackers[j].sigma.wl=trackers[j].sigma.wl>MINSIGMA?trackers[j].sigma.wl:MINSIGMA;

                    trackers[j].mean.wr=(pretracker.mean.wr*pow(trackers[j].sigma.wr,2)+trackers[j].mean.wr*pow(pretracker.sigma.wr,2))/(pow(pretracker.sigma.wr,2)+pow(trackers[j].sigma.wr,2));
                    trackers[j].sigma.wr=(pow(pretracker.sigma.wr,2)*pow(trackers[j].sigma.wr,2))/(pow(pretracker.sigma.wr,2)+pow(trackers[j].sigma.wr,2));
                    trackers[j].sigma.wr=trackers[j].sigma.wr>MINSIGMA?trackers[j].sigma.wr:MINSIGMA;

                    trackers[j].mean.lf=(pretracker.mean.lf*pow(trackers[j].sigma.lf,2)+trackers[j].mean.lf*pow(pretracker.sigma.lf,2))/(pow(pretracker.sigma.lf,2)+pow(trackers[j].sigma.lf,2));
                    trackers[j].sigma.lf=(pow(pretracker.sigma.lf,2)*pow(trackers[j].sigma.lf,2))/(pow(pretracker.sigma.lf,2)+pow(trackers[j].sigma.lf,2));
                    trackers[j].sigma.lf=trackers[j].sigma.lf>MINSIGMA?trackers[j].sigma.lf:MINSIGMA;

                    trackers[j].mean.lb=(pretracker.mean.lb*pow(trackers[j].sigma.lb,2)+trackers[j].mean.lb*pow(pretracker.sigma.lb,2))/(pow(pretracker.sigma.lb,2)+pow(trackers[j].sigma.lb,2));
                    trackers[j].sigma.lb=(pow(pretracker.sigma.lb,2)*pow(trackers[j].sigma.lb,2))/(pow(pretracker.sigma.lb,2)+pow(trackers[j].sigma.lb,2));
                    trackers[j].sigma.lb=trackers[j].sigma.lb>MINSIGMA?trackers[j].sigma.lb:MINSIGMA;
                }
            }
        }
    }
    //==============================
    delete []h_weight;
    CUDAFREE(d_weight);
    delete []h_sampleids;
    CUDAFREE(d_sampleids);

    delete []h_controls;
    CUDAFREE(d_controls);

    delete []h_particles;
    CUDAFREE(d_particles);
    CUDAFREE(d_tmpparticles);
}
