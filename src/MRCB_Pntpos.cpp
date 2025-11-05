#include "MRCB.h"
/* constants */
#define NX_SPP      (3+NSYS)   /* # of estimated parameters */
#define ERR_ION     5.0            /* ionospheric delay std (m) */
#define ERR_BRDCI   0.5            /* broadcast iono model error factor */
#define MAXITR      10             /* max number of iteration for point pos */

/* pseudorange measurement error variance */
static double varerr(const prcopt_t *opt,double el,int sys)
{
    double fact=1.0,varr;
    varr = SQR(100.0 * opt->err[1] / sin(el));
    return SQR(fact) * varr;
}
/* psendorange with code bias correction */
static double prange(const obsd_t *obs,const nav_t *nav,double el,int iter,const prcopt_t *opt)
{
    double P1 = 0.0;
    int sys, prn;
    if (!(sys = satsys(obs->sat, &prn))) return 0.0;
    P1 = obs->P[0];
    return P1;
}
/* ionospheric correction */
static int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var, const prcopt_t* opt)
{
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel,opt);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    *ion = 0.0;
    *var = ionoopt == IONOOPT_OFF ? SQR(ERR_ION) : 0.0;
    return 1;
}
/* tropospheric correction */
static int tropcorr(gtime_t time, const nav_t* nav, const double* pos,
                    const double* azel, int tropopt, double* trp, double* var)
{
    /* Saastamoinen model */
    if (tropopt == TROPOPT_SAAS || tropopt == TROPOPT_EST || tropopt == TROPOPT_ESTG) {
        *trp = tropmodel(time, pos, azel, REL_HUMI);
        *var = SQR(ERR_SAAS / (sin(azel[1]) + 0.1));
        return 1;
    }
    /* no correction */
    *trp = 0.0;
    *var = 0.0;
    return 1;
}

/* pseudorange residuals */
static int rescode(int iter,const prcopt_t *opt,const obsd_t *obs,int n, const nav_t *nav,const double *rr,const double *dtr,
                   rscs_t *rscs,double *v,double *H,double *var,int *vsat,int *exc)
{
    double r,dion=0,vion=0,dtrp=0,vtrp=0.0,pos[3]={0},e[3]={0},P;
    int i,j,k,nv=0,sys,prn,sat,year=0;
    gtime_t time;
    ecef2pos(rr,pos);
    for (i=0;i<n&&i<MAXRCV*MAXOBS;++i) {
        /*initialize*/
		vsat[i]=0;
        dtrp = dion = 0.0;
        rscs[i].azel[0] = rscs[i].azel[1] = 0.0;
        sat = obs[i].sat;
        time = obs[i].time;
        /*check*/
		if(sat<=0) continue;
        if (!(sys = satsys(sat, &prn))) continue;
        if (satexclude(sat, rscs[i].var, rscs[i].svh, opt)) continue;
        if (opt->bds2 == 1) {
            if (sys == SYS_CMP && (prn > 18 || prn > 46 || prn == 31)) continue;
        }
        else if (opt->bds2 == 0) {
            if (sys == SYS_CMP && (prn <= 18 || prn > 46 || prn == 31)) continue;
        }
        if (exc[i]) continue;
        /*compute*/
    	k=sys2num(sys)-1;
    	/* geometric distance/azimuth/elevation angle */
    	if ((r=geodist(rscs[i].rs,rr,e))<=0.0||satazel(pos,e,rscs[i].azel)<opt->elmin) continue;
    	/* psudorange with code bias correction */
    	if ((P=prange(obs+i,nav,rscs[i].azel[1],iter,opt))==0.0) continue; 
        vsat[i] = 1;

        if (iter > 0) {
            /* ionospheric corrections   ionospheric delay (L1) (m) */
            if (!ionocorr(time, nav, sat, pos, rscs[i].azel, IONOOPT_BRDC, &dion, &vion, opt)) continue;
            /* tropospheric corrections */
            if (!tropcorr(time, nav, pos, rscs[i].azel, TROPOPT_SAAS, &dtrp, &vtrp)) {
                continue;
            }
        }
    	/* pseudorange residual */
    	v[nv]=P-(r+dtr[k]-CLIGHT*rscs[i].dts[0]+dion+dtrp);
    	/* design matrix */
    	for (j=0;j<NX_SPP;++j) H[nv*NX_SPP+j]=j<3?-e[j]:0.0;
    	H[nv*NX_SPP+k+3]=1.0;
        var[nv++]=varerr(opt,rscs[i].azel[1],sys);
    }
    return nv;
}
/* validate solution */
static int valsol(const rscs_t *rscs,const int *vsat,int n,const prcopt_t *opt,
                  const double *v,int nv,int nx,char *msg, sol_t* sol)
{
    double azels[MAXOBS*2],dop[4];
    int i,ns;
    /* large gdop check */
    for (i=ns=0;i<n;++i) {
        if (!vsat[i]) continue;
        azels[  ns*2]=rscs[i].azel[0];
        azels[1+ns*2]=rscs[i].azel[1];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    for (i = 0; i < 4; i++) {
        sol->dop[i] = dop[i];
    }
    if (dop[0]<=0.0||dop[0]>30.0) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position */
static int estpos(const obsd_t *obs,int n,const nav_t *nav,const prcopt_t *opt,rscs_t *rscs,sol_t *sol,int *vsat,char *msg)
{
    double sig,d,dx[NX_SPP],rr[3],dtr[NSYS]={0};
	double Q[NX_SPP*NX_SPP],v[MAXOBS],V[MAXOBS],H[MAXOBS*NX_SPP],H_[MAXOBS*NX_SPP],HH[MAXOBS*NX_SPP],var[MAXOBS];
	int i,j,k,nv,nx,ix[NX_SPP],exc[MAXRCV*MAXOBS]={0},info,stat;

    for (i = 0; i < 3; ++i) rr[i] =sol->rr[i];
	for (i=0;i<MAXITR;++i) {
		memset(v,0,sizeof(v)); memset(dx,0,sizeof(dx)); memset(var,0,sizeof(var));
		memset(H,0,sizeof(H)); memset(H_,0,sizeof(H_));
		memset(Q,0,sizeof(Q)); memset(ix,0,sizeof(ix));
		memset(V,0,sizeof(V)); memset(HH,0,sizeof(HH));
		/* pseudorange residuals */
		nv=rescode(i,opt,obs,n,nav,rr,dtr,rscs,v,H,var,vsat,exc);
		nx=0;
		for (j=0;j<NX_SPP;++j) {
			d=0.0;
			for (k=0;k<nv;++k) d+=fabs(H[k*NX_SPP+j]);
			if (d<=0) continue;
			ix[nx]=j;
			++nx;
		}
		if (nv<nx) {
			sprintf(msg,"lack of valid sats ns=%d",nv);
			break;
		}

		/* weight by variance */
		for (j=0;j<nv;++j) {
			sig=sqrt(var[j]);
			V[j] =v[j];
			v[j]/=sig;
			for (k=0;k<nx;++k) H_[j*nx+k]=H[j*NX_SPP+ix[k]]/sig;
			for (k=0;k<nx;++k) HH[j*nx+k]=H[j*NX_SPP+ix[k]];
		}

		/* least square estimation */
		if ((info=lsq(H_,v,nx,nv,dx,Q))) {
			sprintf(msg,"lsq error info=%d",info);
			break;
		}
		for(j=0;j< 3;++j) rr[j]       +=dx[j];
		for(j=3;j<nx;++j) dtr[ix[j]-3]+=dx[j];  
		if (i>1&&norm(dx,NX_SPP)<1E-4) {
			if(sol){
				for (j=0;j<6   ;++j) sol->rr[j]=j<3?rr[j]:0.0;
			    for (j=0;j<NSYS;++j) sol->dtr[j]=dtr[j]/CLIGHT;
			    sol->ns=(unsigned char)nv; 
			}
			/* validate solution GDOP */
			if ((stat=valsol(rscs,vsat,n,opt,v,nv,NX_SPP,msg,sol))) {
				if(sol) sol->stat=SOLQ_SINGLE;
			}
			return stat;
		}

	}
    if (i >= MAXITR) {
        sprintf(msg, "SPP iteration divergent i=%d", i);
    }
	return 0;
}
/* single-point positioning */
extern int pntpos(const obsd_t *obs,int n,const nav_t *nav,const prcopt_t *opt, sol_t *sol,rscs_t *rscs,ssat_t *ssat,char *msg)
{
    int i,stat,vsat[MAXRCV*MAXOBS]={0};
    if(sol) sol->stat=SOLQ_NONE;
    if (n<=0) {
        strcpy(msg,"no observation data");
        return 0;
    }
    msg[0]='\0';
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,nav,opt,rscs,sol,vsat,msg);
	for(i=0;i<MAXSAT;++i) ssat[i].vs=0;
    if (ssat) {
        for (i=0;i<n&&i<MAXRCV*MAXOBS;++i) if (vsat[i]) ssat[obs[i].sat-1].vs=1;
    }
    return stat;
}