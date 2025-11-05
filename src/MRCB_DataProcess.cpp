#include "MRCB.h"
#define THRES_MW_JUMP 10.0

/* initialize state and covariance */
static void initx(rcb_t*rcb,int i)
{
    int j;
    rcb->x[i]=0;
    for (j=0;j<rcb->nx;j++) rcb->P[i+j*rcb->nx]=rcb->P[j+i*rcb->nx]=0.0;
}
/* geometry-free phase measurement */
static double gfmeas(const obsd_t *obs,const nav_t *nav)
{
	
	if (obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;
	return obs->L[0]*nav->lam[obs->sat-1][0]-obs->L[1]*nav->lam[obs->sat-1][1];
}
/* Melbourne-Wubbena linear combination */
static double mwmeas(const obsd_t *obs,const nav_t *nav)
{
	double freq1,freq2;
	freq1=CLIGHT/nav->lam[obs->sat-1][0];
    freq2=CLIGHT/nav->lam[obs->sat-1][1];
    if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0||obs->P[0]==0.0||obs->P[1]==0.0) return 0.0;
    return (obs->L[0]-obs->L[1])*CLIGHT/(freq1-freq2)-(freq1*obs->P[0]+freq2*obs->P[1])/(freq1+freq2);
}
/* detect cycle slip by geometry free phase jump */
static void detslp_gf(const obsd_t* obs, int n, const nav_t* nav, int nf, rcb_t* rcb)
{
	int sys, prn;
	double g0, g1;
	int i, j, maxI = 0;
	double yuzhi = 0.05;
	double max = 0;
	int cntMAX = 0, cnt = 0;
	for (i = 0; i < n && i < MAXOBS; i++) {
		sys = satsys(obs[i].sat, &prn);
		if (sys == 4) cntMAX++;
	}
	for (i = 0; i < n && i < MAXOBS; i++) {
		sys = satsys(obs[i].sat, &prn);
		if ((g1 = gfmeas(obs + i, nav)) == 0.0) continue;
		g0 = rcb->ssat[obs[i].sat - 1].gf[0];
		rcb->ssat[obs[i].sat - 1].gf[0] = g1;
		if (sys == 4) {
			if (max < fabs(g1 - g0) && g0 != 0.0) {
				max = fabs(g1 - g0);
				maxI = obs[i].sat;
			}
			cnt++;
			if (cnt != cntMAX) continue;
		}
		if (sys == 4 && g0 != 0.0 && max > 0.5) {
			for (j = 0; j < rcb->opt.nf; j++) rcb->ssat[maxI - 1].slip[j] |= 1;
		}
		else if (g0 != 0.0 && fabs(g1 - g0) > yuzhi) {
			for (j = 0; j < rcb->opt.nf; j++) rcb->ssat[obs[i].sat - 1].slip[j] |= 1;
		}
	}
}
/* detect slip by Melbourne-Wubbena linear combination jump */
static void detslp_mw(const obsd_t *obs,int n,const nav_t *nav, int nf,rcb_t* rcb)
{
	double w0, w1;
	int i, j;
	int sys, prn;
	for (i = 0; i < n && i < MAXOBS; i++) {
		sys = satsys(obs[i].sat, &prn);
		if ((w1 = mwmeas(obs + i, nav)) == 0.0) continue;
		w0 = rcb->ssat[obs[i].sat - 1].mw[0];
		rcb->ssat[obs[i].sat - 1].mw[0] = w1;
		if (w0 != 0.0 && fabs(w1 - w0) > THRES_MW_JUMP) {
			for (j = 0; j < rcb->opt.nf; j++) rcb->ssat[obs[i].sat - 1].slip[j] |= 1;
		}
	}
}

/* cal rcb */
extern int calrcb(FILE**fp,rcb_t*rcb,obsd_t*obs,int n, nav_t *nav,sta_t*sta, const solopt_t* sopt)
{
    prcopt_t *opt=&rcb->opt;
	gtime_t time;
    rscs_t rscs[MAXRCV*MAXOBS];
	obsd_t odata[MAXRCV*MAXOBS]={0};
	double ep[6]={0.0};
	int i,j,t=0;
    char msg[128]="";

	for (j = 0; j < n && j < MAXRCV * MAXOBS; ++j) {
		odata[t++] = obs[j];
	}
    /* previous epoch */
    time=rcb->sol.time;
	/* epoch time */
	rcb->sol.time = odata[0].time;
    for (i=0;i<n&&i<MAXRCV*MAXOBS;++i) rscs[i].iode=-1;
	/* satellite positions and clocks */
	if (satposs(rcb->sol.time, opt, odata, n,  nav, opt->sateph, rscs) < 4) {
		printf("station: stellite clock lack!\n");
		return 0;
	}
	/* rover position by single point positioning */
	if (!pntpos(odata, n, nav, &rcb->opt, &rcb->sol, rscs, rcb->ssat, msg)) {
		printf("SPP point pos error:%s, Possible poor data quality. \n", msg);
		printf("You can try using a multi-system solution to increase the number of satellites and resolve this issue.");
		return 0;
	}
	/*cycle slip detection*/
	for(i=0;i<MAXSAT;++i) {
	    for(j=0;j < NFREQ;++j) rcb->ssat[i].slip[j]=0;
	}
	detslp_gf(odata, n, nav, NFREQ, rcb);
	detslp_mw(odata, n, nav, NFREQ, rcb);

    int stat=RCB(fp,rcb,odata, n,nav,rscs,sta,sopt);
    return stat;
}
