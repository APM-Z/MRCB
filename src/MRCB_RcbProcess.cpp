#include "MRCB.h"
#include "MRCB_SolBody.h"
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
/* measurement error variance */
static double varerr(double el,int type,int sys,const prcopt_t *opt,int frq)
{
    double fact=1.0,sinel=sin(el);
	if (type == 1) fact = 100.0;
	else           fact = 1.0;
    switch (sys) {
        case SYS_GPS: fact*=EFACT_GPS; break;
        case SYS_GLO: fact*=EFACT_GLO; break;
        case SYS_GAL: fact*=EFACT_GAL; break;
        case SYS_CMP: fact*=EFACT_CMP; break;
		case SYS_QZS: fact*=EFACT_QZS; break;
    }
	if (sys == SYS_GPS || sys == SYS_QZS) {
		if (frq == 2) {
			fact *= EFACT_GPS_L5;
		}
	}

    return SQR(fact * opt->err[1] /sinel);
}


/* initialize state and covariance */
static void initx(rcb_t* rcb, double xi, double var, int i)
{
    int j;
    rcb->x[i]=xi;
	for (j = 0; j < rcb->nx; j++) {
		rcb->P[i + j * rcb->nx] = rcb->P[j + i * rcb->nx] = i == j ? var : 0.0;
	}
}
/* temporal update of frequency-independent parameter */
static void uprex_mrcb(rcb_t* rcb, const obsd_t* obs, int n, unsigned char* F) {
	int i, j, sat;
	prcopt_t* opt = &rcb->opt;
	
	for (i = 0; i < n && i < MAXOBS; i++) {
		if ((sat = obs[i].sat) <= 0) {
			continue;
		}
		j = sat - 1;
		initx(rcb, 0, 0, j);
		F[j] = 0;
	}
}
/* temporal update of rcb */
static void urcb_mrcb(rcb_t* rcb, const obsd_t* obs, int n, unsigned char* F)
{
	int i, j, f;
	prcopt_t* opt = &rcb->opt;
	for (i = 0; i < NSYS; i++) {
		for (f = 0; f < opt->nf; f++) {
			j = IM(&rcb->opt, i, f);
			initx(rcb, 0, 0, j);
			F[j] = 0;
		}
	}
}
/* temporal update of scb */
static void uscb_mrcb(rcb_t* rcb, const obsd_t* obs, int n, unsigned char* F)
{
	int i, j, f, sat;
	prcopt_t* opt = &rcb->opt;

	for (f = 0; f < opt->nf; ++f) {
		if (f <= 1) continue;
		for (i = 0; i < n && i < MAXOBS; i++) {
			if ((sat = obs[i].sat) <= 0) {
				continue;
			}
			j = SCBN(&rcb->opt, sat, f);
			/*first epoch at a arc*/
			if (rcb->x[j] == 0 || !rcb->cvsat[sat - 1])
			{
				initx(rcb, 0, 0, j);
				F[j] = 0;
			}
			else
			{
				rcb->P[j + j * rcb->nx] += SQR(rcb->opt.prn[0]) * fabs(rcb->tt);
				F[j] = 1;
			}
		}
	}
}
/* temporal update of GLONASS IFB */
static void uIFB_mrcb(rcb_t* rcb, const obsd_t* obs, int n, unsigned char* F)
{
	int j, f;
	prcopt_t* opt = &rcb->opt;

	for (f = 0; f < opt->nf; f++) {
		j = IFB(&rcb->opt, f);
		initx(rcb, 0, 0, j);
		F[j] = 0;
	}
}
/* temporal update of ion parameters */
static void udion_mrcb(rcb_t* rcb, const obsd_t* obs, int n, unsigned char* F)
{
	int i, j, sat;
	for (i = 0; i < n && i < MAXOBS; i++) {
		if ((sat = obs[i].sat) <= 0) continue;
		j = II(&rcb->opt, sat);
		if (rcb->x[j] == 0.0 || !rcb->cvsat[sat - 1]||rcb->opt.Initial==0) {
			initx(rcb,0,0,j);
			F[j] = 0;
		}
		else {
			initx(rcb, 0, 0, j);
			F[j] = 0;
		}
	}
}

/* temporal update of phase biases */
static void udamb_mrcb(rcb_t *rcb,const obsd_t *obs,int n,unsigned char *F)
{
    int i,j,f,sat;
	prcopt_t* opt = &rcb->opt;

	for(f=0;f<opt->nf;++f){
		for(i=0;i<n&&i<MAXOBS;i++) {
			if ((sat = obs[i].sat) <= 0) {
				continue;
			}
	        j=IA(&rcb->opt,sat,f);
	        /*first epoch at a arc*/
	        if(rcb->x[j]==0||rcb->ssat[sat-1].slip[f]||!rcb->cvsat[sat-1]||rcb->opt.Initial==0)
			{
				initx(rcb, 0, 0, j);
		        F[j]=0;
			}
			else 
			{
				rcb->P[j + j * rcb->nx] += SQR(rcb->opt.prn[0]) * fabs(rcb->tt);
				F[j] = 1;
			}
        }
	}
}

/* temporal update of states */
static void udstate_rcb_sd(rcb_t*rcb,const obsd_t*obs,int n,unsigned char *F)
{
	int i;
	prcopt_t* opt = &rcb->opt;
	/* Initial */
	for (i = 0; i < rcb->nx; ++i) F[i] = 0;
	/*temporal update of P*/
	uprex_mrcb(rcb, obs, n, F);
	/*temporal update of ionoopt parameters*/
	if (rcb->opt.ionoopt == IONOOPT_EST) {
		udion_mrcb(rcb, obs, n, F);
	}
	/*temporal update of rcb*/
	urcb_mrcb(rcb, obs, n, F);
	/*temporal update of scb*/
	if (opt->nf >= 3) uscb_mrcb(rcb, obs, n, F);
	/*temporal update of ifb*/
	uIFB_mrcb(rcb, obs, n, F);
	/*temporal update of ambiguity*/
	udamb_mrcb(rcb, obs, n, F);
}

/* satellite antenna phase center variation */
static void satantpcv(const double *rs,const double *rr,const pcv_t *pcv,double *dant)
{
    double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
    int i;

    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;

    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    antmodel_s(pcv,nadir,dant);
}

/* phase and code residuals */
static int rcb_unimbined(FILE** fp, const obsd_t *obs,int n,const nav_t*nav,rcb_t*rcb,sta_t*sta,rscs_t*rscs,double *v,double *H,double *R,int *ih,int *ifrq)
{

	prcopt_t* opt = &rcb->opt;
	double var[MAXRCV * MAXOBS * 2 * NFREQ + MAXRCV * MAXOBS + MAXRCV] = { 0 };
	double rr[3] = { 0 }, e[3] = { 0 }, ep[6] = { 0 };
	double pos[3] = { 0 };
	double P[NFREQ] = { 0 }, L[NFREQ] = { 0 };
	double C, y, r, dion = 0, rcbs = 0, bias = 0, scb = 0, IFB = 0;
	int i, j, k, f, s[MAXRCV * MAXOBS] = { 0 }, sat, sys, prn, nv = 0, ns = 0;
	int nx = rcb->nx;
	unsigned char flag, stat = 1;/*0: phase        1:code*/
	char buff[8192], * p = buff;
	char timestr[64], id[4];
	gtime_t time=rcb->sol.time;
	/*time*/
	time2str3(time, timestr, 2);
	p = buff;
	p += sprintf(p, "%s", timestr);
	p += sprintf(p, "\n");
	fwrite(buff, (int)(p - buff), 1, fp[2]);

	for (i = 0; i < MAXSAT; ++i) {
		rcb->cvsat[i] = 0;
		for (j = 0; j < NFREQ; ++j) {
			rcb->ssat[i].vsat[j] = 0;
			rcb->sol.resp[i][j] = 0.0;
			rcb->sol.resc[i][j] = 0.0;
		}
	}
	for (i = 0; i < 3; i++) rr[i] = rcb->sol.rr[i];
	ecef2pos(rr, pos);

	for (i = 0; i < n && i < MAXOBS; ++i) {
		/* check */
		if ((sat = obs[i].sat) <= 0) continue;
		satno2id(sat, id);
		//s[i] G R E Q C -> 0 1 2 3 4
		if (!(sys = satsys(sat, &prn)) || (s[i] = sys2num(sys) - 1) < 0) continue;
		if (satexclude(sat, rscs[i].var, rscs[i].svh, opt)) continue;
		if ((r = geodist(rscs[i].rs, rr, e)) <= 0.0) continue;
		if (satazel(pos, e, rscs[i].azel) < opt->elmin) {
			OutLogElv(fp, opt->elmin * R2D, id, rscs[i].azel[1] * R2D); //LOG file
			continue;
		}
		if (!rcb->ssat[sat - 1].vs) continue;
		for (j = 0, f = 0; j < opt->nf; j++) {
			if (obs[i].L[j] == 0.0 || obs[i].P[j] == 0.0 || nav->lam[sat - 1][j] == 0.0) continue;
			if (obs[i].SNR[j] != 0 && testsnr(1, 0, rscs[i].azel[1], obs[i].SNR[j] * 0.25, &opt->snrmask)) continue;
			f++;
		}
		if (f < opt->nf) {
			OutLogFrq(fp, id, opt->nf, f); //LOG file
			continue;
		}

		rcb->cvsat[sat - 1] = 1;

		/* corrected phase and code measurements */
		for (j = 0; j < opt->nf; ++j) {
			L[j] = P[j] = 0.0;
			if (nav->lam[sat - 1][0] == 0.0 || nav->lam[sat - 1][j] == 0.0) continue;
			if (obs[i].L[j] == 0.0 || obs[i].P[j] == 0.0) continue;
			L[j] = obs[i].L[j] * nav->lam[sat - 1][j];
			P[j] = obs[i].P[j];
		}
		/* stack phase and code residuals {L1 P1 L2 P2} */
		for (j = 0; j < 2 * opt->nf; ++j) {
			IFB = dion = scb = rcbs = bias = 0.0;
			f = j / 2; flag = j % 2;//flag(0--L 1--P)
			if (nav->lam[sat - 1][0] == 0 || nav->lam[sat - 1][f] == 0) continue;
			if (L[f] == 0 || P[f] == 0) continue;
			if ((y = flag == 0 ? L[f] : P[f]) == 0.0) continue;
			C = SQR(nav->lam[sat - 1][f] / nav->lam[sat - 1][0]);
			for (k = 0; k < nx; k++) {
				if (k == (sat - 1)) {
					H[k + nx * nv] = 1.0;
				}
				else {
					H[k + nx * nv] = 0.0;
				}
			}
			if (opt->ionoopt == IONOOPT_EST) {
				H[II(opt, sat) + nx * nv] = C * (flag == 0 ? -1.0 : 1.0);
				dion = H[II(opt, sat) + nx * nv] * rcb->sol.ion[sat - 1];
			}
			if (opt->Initial != 0 && flag == 1) {
				H[IM(opt, s[i], f) + nx * nv] = 1.0;
				rcbs = 1 * rcb->sol.rcb[s[i]][f];
				if (sys == SYS_GLO) {
					for (k = 0; k < nav->ng; k++) {
						if (nav->geph[k].sat != sat) continue;
						H[IFB(opt, f) + nx * nv] = nav->geph[k].frq;
						IFB = H[IFB(opt, f) + nx * nv] * rcb->sol.IFB[f];
					}
				}
			}
			if (opt->nf >= 3 && flag == 1) {
				if (f <= 1) {
					H[SCBN(opt, sat, f) + nx * nv] = 0;
				}
				else {
					H[SCBN(opt, sat, f) + nx * nv] = -1.0;
					scb = rcb->sol.scb[sat - 1][f];
				}
			}
			if (flag == 0) {
				H[IA(opt, sat, f) + nx * nv] = nav->lam[sat - 1][f];
				bias = nav->lam[sat - 1][f] * rcb->sol.amb[sat - 1][f];
			}
			/*residuals*/
			v[nv] = y - (rcb->sol.prex[sat - 1] + dion + rcbs - scb + bias + IFB);
			if (!flag) rcb->sol.resc[sat - 1][f] = v[nv];
			else       rcb->sol.resp[sat - 1][f] = v[nv];
			var[nv] = varerr(rscs[i].azel[1], flag, sys, opt, f);
			rcb->ssat[sat - 1].vsat[f] = 1;
			if (!flag) {
				ih[nv] = (sat);   /*>0:phase*/
			}
			else {
				ih[nv] = -(sat);   /*<0:code*/
			}
			ifrq[nv] = f;
			++nv; 
		}
		++ns;
	}
	rcb->sol.ns = ns;
	for (i = 0; i < nv; ++i) for (j = 0; j < nv; ++j) R[i * nv + j] = i == j ? var[i] : 0.0;
	return nv;
}

static void update_pre(rcb_t*rcb,double*x,const prcopt_t*opt,int n, const obsd_t* obs)
{
	int i, j, k, f, sat;
	/*update P*/
	for (i = 0; i < n && i < 2 * MAXOBS; ++i) {
		sat = obs[i].sat;
		rcb->sol.prex[sat - 1] += x[sat - 1];
	}
	/*update ionosphere*/
	if (opt->ionoopt == IONOOPT_EST) {
		for (j = 0; j < MAXSAT; ++j)
		{
			k = II(opt, j + 1);
			rcb->sol.ion[j] += x[k];
		}
	}
	/*update rcb*/
	for (j = 0; j < NSYS; j++) {
		for (f = 0; f < opt->nf; f++) {
			rcb->sol.rcb[j][f] += x[IM(opt, j, f)];
		}
	}
	/*update scb*/
	if (opt->nf >= 3) {
		for (i = 0; i < n && i < 2 * MAXOBS; ++i) {
			sat = obs[i].sat;
			for (f = 0; f < opt->nf; ++f) {
				k = SCBN(opt, sat, f);
				rcb->sol.scb[sat - 1][f] += x[k];
			}
		}
	}
	/*update IFB*/
	for (f = 0; f < opt->nf; f++) {
		k = IFB(opt, f);
		rcb->sol.IFB[f] += x[k];
	}

	/*ambiguity parameters*/
	for (i = 0; i < n && i < 2 * MAXOBS; ++i) {
		sat = obs[i].sat;
		for (f = 0; f < opt->nf; ++f) {
			k = IA(opt, sat, f);
			rcb->sol.amb[sat - 1][f] += x[k];
		}
	}
}

/* sqrt of covariance */
static double sqvar(double covar)
{
    return covar < 0.0 ? -sqrt(-covar) : sqrt(covar);
}

static double getcri(int red)
{
	if     (red<1   ) return 0;
	else if(red<1500) return diacri[red-1];
	else              return 0.98;
}

/* LSF filter -------------------------------------------------------------*
** args   : double *Xp  (n x 1)  I  unknown parameters vector
**          double *Pp  (n x n)  I  covariance matrix of unknown parameters
**          double *H   (m x n)  I  
**          double *v   (m x 1)  I  measurement - model
**          double *R   (m x m)  I  covariance matrix of measurement error
**          int    *ih  (m x 1)  I  
**          char   *ifrq(m x 1)  I  
*			char   *F	(m x n)  I  
*			int	   *ix  (n x 1)  I	
*			double *TqLOM
**          int    n,m           I  number of unknown parameters and measurements
**----------------------------------------------------------------------------*/
static int LSFfilt(rcb_t*rcb,prcopt_t *opt,double *Xp,double *Pp,const double *H,
	               const double *v,const double *R,const int *ih,int *ifrq,const unsigned char *F,
	               double *TqLOM,int n,int m,int *ix,int*redundancy)
{
	int i,j,M,N,im,p,nx,np,nc,red;
	int *ip,*C,T=m;
	unsigned char stat=1;
	double d,LOM,kLOM,tqlom,TQLOM=1,maxwt;
	double *A1,*Y1,*R1,*iR1,*N1,*r1;
	double *A2,*Y2,*R2,*iR2,*N2,*r2;
	double *Am,*Ym,*Rm,*iRm,*X,*S,*V;
	double *Qv,*Qc,*Pk,*W,*wt;
	int exc_sat;
	/* obtain valid number of parameters */
	ip=imat(n,1);
	for (i=nx=np=0;i<n;++i) {
		for (j=0,d=0;j<m;++j) d+=fabs(H[j*n+i]);
		if (d>0) {
			if (F[i] == 1) {
				ip[np++] = i;
			}
			ix[nx++]=i;
		}
	}

	*redundancy=m-nx;
	/* LSF filter */
	if (opt->Initial!=0&&np>0)  {
		M=m+np;  /* number of measurements  */
		N=nx;    /* number of unknown parameters */
		red=M-N;
		if(red<=0) {
			free(ip);
			*TqLOM=99999.9;
			return -1;
		}
		/* N1=A1'*iR1*A1, r1=A1'*iR1*Y1 */
		A1=zeros(np,nx);
		Y1=zeros(np,1);
		R1=zeros(np,np);
		for (i=0;i<np;++i) {
			for (j=0;j<nx;++j) if (ip[i]==ix[j]) A1[i*nx+j]=1;
			Y1[i]=Xp[ip[i]];
			for (j=0;j<np;++j) R1[i*np+j]=Pp[ip[i]*n+ip[j]];
		}
		free(ip);

		// iR1==R1^-1
		iR1=zeros(np,np); matcpy(iR1,R1,np,np);
		if (matinv(iR1,np)==-1) {
			free(A1); free(Y1); free(R1); free(iR1);
			return -1;
		}
		W=zeros(nx,np); N1 =zeros(nx,nx); r1 =zeros(nx,1);
		matmul("NT",nx,np,np,1.0,A1,iR1,0.0,W); /* W=A1'*iR1 */
		matmul("NT",nx,nx,np,1.0,W,A1,0.0,N1);  /* N1=A1'*iR1*A1 */
		matmul("NN",nx,1 ,np,1.0,W,Y1,0.0,r1);  /* r1=A1'*iR1*Y1 */
		free(W);

		/* N2=A2'*iR2*A2, r2=A2'*iR2*Y2 */
		A2=zeros(m,nx);
		Y2=zeros(m,1);
		R2=zeros(m,m);
		for (i=0;i<m;++i) {
			for (j=0;j<nx;++j) A2[i*nx+j]=H[i*n+ix[j]];
			for (j=0;j<m ;++j) R2[i*m +j]=R[i*m+j];
			Y2[i]=v[i];
			if (fabs(R[i * m + i]) < 1e-12) {
				free(A1); free(Y1); free(R1); free(iR1); free(N1); free(r1);
				free(A2); free(Y2); free(R2);
				return -1;
			}
		}

		iR2=zeros(m,m); matcpy(iR2,R2,m,m);
		if (matinv(iR2,m)==-1) {
			free(A1); free(Y1); free(R1); free(iR1); free(N1); free(r1);
			free(A2); free(Y2); free(R2); free(iR2);
			return -1;
		}
		W=zeros(nx,m); N2 =zeros(nx,nx); r2 =zeros(nx,1);
		matmul("NT",nx,m ,m,1.0,A2,iR2,0.0,W); /* A2'*iR2 */
		matmul("NT",nx,nx,m,1.0,W,A2,0.0,N2);  /* A2'*iR2*A2 */
		matmul("NN",nx,1 ,m,1.0,W,Y2,0.0,r2);  /* A2'*iR2*Y2 */
		free(W);
		/* S=(N1+N2)^-1 */
		S=zeros(nx,nx);
		for (i=0;i<nx;++i) for (j=0;j<nx;++j) {
			 S[i*nx+j]=N1[i*nx+j]+N2[i*nx+j];
		}
		free(N1); free(N2);
		if (matinv(S,nx)==-1) { /* S=(N1+N2)^-1 */
			free(S);
			free(A1); free(Y1); free(R1); free(iR1); free(r1);
			free(A2); free(Y2); free(R2); free(iR2); free(r2);
			return -1;
		}
		/* X=S*(r1+r2) */
		X=zeros(nx,1);
		W=zeros(nx,1);
		for (i=0;i<nx;++i) W[i]=r1[i]+r2[i];
		matmul("TN",nx,1,nx,1.0,S,W,0.0,X); /* X=Q*(r1+r2) */
		free(r1); free(r2); free(W);
		/* V=[Y1-A1*X; Y2-A2*X] */
		V=zeros(M,1);
		W=zeros(np,1);
		matcpy(W,Y1,np,1);
		matmul("TN",np,1,nx,-1.0,A1,X,1.0,W);
		for (i=0;i<np;++i) V[i]=W[i];
		free(W);

		W=zeros(m,1);
		matcpy(W,Y2,m,1);
		matmul("TN",m,1,nx,-1.0,A2,X,1.0,W);
		for (i=0;i<m;++i) V[i+np]=W[i];
		free(W);

		/* iR=blkdiag(iR1,iR2), Rm=blkdiag(R1,R2) */
		//iR1=R1^-1,iR2=R2^-1
		iRm=zeros(M,M); Rm=zeros(M,M);
		for (i=0;i<np;++i) for (j=0;j<np;++j) {
			iRm[i*M+j]=iR1[i*np+j];
			Rm[i*M+j]=R1[i*np+j];
		}
		for (i=0;i<m;++i) for (j=0;j<m;++j) {
			iRm[(i+np)*M+(j+np)]=iR2[i*m+j];
			Rm[(i+np)*M+(j+np)]=R2[i*m+j];
		}
		free(R1); free(R2); free(iR1); free(iR2);


		/* Am=[A1; A2], Ym=[Y1; Y2] */
		Am=zeros(M,N); Ym=zeros(M,1);
		for (i=0;i<np;++i) {
			for (j=0;j<N;++j) Am[i*N+j]=A1[i*N+j];
			Ym[i]=Y1[i];
		}
		for (i=0;i<m;++i) {
			for (j=0;j<N;++j) Am[(i+np)*N+j]=A2[i*N+j];
			Ym[i+np]=Y2[i];
		}
		free(A1); free(A2); free(Y1); free(Y2);

		/* LOM = V'*iR*V/red */
		W = zeros(M, 1);
		matmul("NT", 1, M, M, 1.0, V, iRm, 0.0, W);
		matmul("NN", 1, 1, M, 1.0, W, V, 0.0, &LOM);
		free(W);
		LOM = LOM / red;
		kLOM = getcri(red);
		*TqLOM = tqlom = LOM / kLOM;
		Qv= mat(M,M); 
		Pk= mat(M,M);
		Qc= mat(M,M);
		wt= mat(m,1);
		C =imat(m,1); memset(C,0,sizeof(int)*m);
		/* QC */
		for (p=1,nc=0;p<min(red-1,m)&&fabs(tqlom)>TQLOM;++p) {
			memset(Qv,0,sizeof(double)*M*M);
            memset(Pk,0,sizeof(double)*M*M);
            memset(Qc,0,sizeof(double)*M*M);
            memset(wt,0,sizeof(double)*m);
			N=nx+nc;
			/* Qv = Rm - Am*S*Am' */
		    W=zeros(M,N);
		    matmul("TT",M,N,N,1.0,Am,S ,0.0,W);   /*      Am*S     */
		    matcpy(Qv,Rm,M,M);
		    matmul("NN",M,M,N,-1.0,W,Am,1.0,Qv); /* Rm - Am*S*Am' */
		    free(W);
			/* Pk = iRm*Qv*iRm */
		    matmul("NN",M,M,M,1.0,iRm,Qv,0.0,Qc); /* iRm*Qv */
		    matmul("NN",M,M,M,1.0,Qc,iRm,0.0,Pk); /* iRm*Qv*iRm */
			/* DIA: detect */
			im=-1; maxwt=0;
			for (i=0;i<T;++i) {
				if (C[i] == 1) continue; 
			    if ((d=Pk[(i+np)*M+(i+np)])<1e-4){ 
					wt[i]=0;
					continue;
				}
			    if (ih[i]>0) { /* phase */
				    for (j=0,wt[i]=0;j<M;++j) wt[i]+=iRm[(i+np)*M+j]*V[j];
				    wt[i]/=sqrt(d);
				    if (fabs(wt[i])>maxwt) {
						maxwt=fabs(wt[i]); 
						im=i;
					}
			    }
				else {         /* code */
					wt[i]=0;
				}
		    }
			if(im==-1||maxwt<5.0) break;
			/* DIA identify */
			if(ih[im]>0) { /* phase */
		        exc_sat=ih[im]%MAXSAT;
				/* cycle slip */
				rcb->ssat[exc_sat-1].slip[ifrq[im]]=1; 
				C[im]=1;
			}
			else {         /* code */
				C[im]=1;
			}
			nc+=1;
			/* DIA: adaptation */
			N=nx+nc; red=M-N;
			if(red<1){
				free(Rm); free(iRm); free(Qv); free(Pk); free(Qc); free(wt); free(C);
				free(Am); free(Ym);  free(X);  free(S);  free(V);
				return -1;
			}
			W=zeros(M,N);
			for (i=0;i<M;++i) {
				for (j = 0; j < N - 1; ++j) {
					W[i * N + j] = Am[i * (N - 1) + j];
				}
			}
			W[(im+np)*N+(N-1)]=1;
			free(Am); 
			Am=zeros(M,N);
			matcpy(Am,W,M,N);
			free(W);		
			/*lsq*/
			free(X); X=zeros(N,1);
			free(S); S=zeros(N,N);
			memset(V,0,sizeof(double)*M);
			if(lsq_(Am,iRm,Ym,N,M,X,S,V)==-1){
				free(Rm); free(iRm); free(Qv); free(Pk); free(Qc); free(wt); free(C);
				free(Am); free(Ym);  free(X);  free(S);  free(V);
				return -1;
			}
			/*computation of Local Overall Model (LOM) test*/	
			W=zeros(1,M);
			matmul("NT",1,M,M,1.0,V,iRm,0.0,W);     /* V'*P   */
			matmul("NN",1,1,M,1.0,W,V  ,0.0,&LOM);  /* V'*P*V */
			free(W);
			LOM=LOM/red; 
			kLOM=diacri[red-1]; 
			*TqLOM=tqlom=LOM/kLOM;
		}

		for(i=0;i<n;++i){
		    Xp[i]=0.0;
	        for(j=0;j<n;++j) Pp[i*n+j]=0.0;
	    }
	    for (i=0;i<nx;++i) {
		    Xp[ix[i]]=X[i];
		    for (j=0;j<nx;++j) Pp[ix[i]*n+ix[j]]=S[i*N+j];
	    }
		free(Rm); free(iRm); free(Qv); free(Pk); free(Qc); free(wt); free(C);free(Am); free(Ym);  free(X);  free(S);  free(V);
		return nx;
	}
	else { /* lsq */
		free(ip);
		M = m;  /* number of measurements */
		N = nx; /* number of unknown parameters */
		red = M - N;
		if (red < 0) return -1;

		Am = zeros(M, N); Ym = zeros(M, 1); Rm = zeros(M, M); X = zeros(N, 1); S = zeros(N, N); V = zeros(M, 1);
		for (i = 0; i < M; ++i) {
			for (j = 0; j < N; ++j)
			{
				Am[i * N + j] = H[i * n + ix[j]]; 
			}
			Ym[i] = v[i]; 
			for (j = 0; j < M; ++j) Rm[i * M + j] = R[i * M + j];
		}
		// Rm
		if (matinv(Rm, M) == -1) {
			free(Am); free(Ym); free(Rm); free(X); free(S); free(V);
			return -1;
		}
		if (lsq_(Am, Rm, Ym, N, M, X, S, V) != 0) {
			free(Am); free(Ym); free(Rm); free(X); free(S); free(V);
			return -1;
		}
		/* TqLOM */
		if (red > 0) {
			for (i = 0, d = 0; i < M; ++i) d += V[i] * V[i];
			LOM = d / red;
			kLOM = diacri[red - 1];
			*TqLOM = LOM / kLOM;
		}
		else {
			*TqLOM = 99999.99;
		}
		for(i=0;i<n;++i){
		    Xp[i]=0.0;
	        for(j=0;j<n;++j) Pp[i*n+j]=0.0;
	    }

	    for (i=0;i<nx;++i) {
		    Xp[ix[i]]=X[i];
		    for (j=0;j<nx;++j) Pp[ix[i]*n+ix[j]]=S[i*N+j];
	    }
		free(Am); free(Ym); free(Rm); free(X); free(S); free(V);
	    return nx;
	}
}


/* number of estimated states */
extern int rtknx(const prcopt_t *opt)
{
    return NX(opt);
}
/* RCB Calculate */
extern int RCB(FILE**fp,rcb_t* rcb,const obsd_t *obs, int n,const nav_t *nav,rscs_t *rscs,sta_t*sta, const solopt_t* sopt)
{
    prcopt_t *opt=&rcb->opt; 
	int stat= SOLQ_SINGLE,nix,nx=rcb->nx,ny=n*NFREQ*2+n+1,*ih,*ix,*ifrq;
	double* v, * H, * R, * Xp, * Pp, dr[3] = { 0 };
	int exc[MAXSAT][NFREQ]={0};
	rcb->sol.redundancy=0; 
	Xp = zeros(nx, 1);
	Pp = zeros(nx, nx);
	/* prefit residuals */ 
	v = zeros(ny, 1); H = zeros(nx, ny); R = zeros(ny, ny); ih = imat(ny, 1); ifrq = imat(ny, 1);
	if (!(ny = rcb_unimbined(fp, obs, n, nav, rcb, sta, rscs, v, H, R, ih, ifrq))) {
		free(v); free(H); free(R); free(ih); free(ifrq);free(Xp); free(Pp);
		return stat;
	}
	/* temporal update  */
	udstate_rcb_sd(rcb, obs, n, rcb->F);
	/*parameters to Xp Pp*/
	//matcpy(Xp, rcb->x, nx, 1); 
	matcpy(Pp, rcb->P, nx, nx);
	/* parameter estimation */
	ix = imat(nx, 1); nix = 0;
	if ((nix = LSFfilt(rcb, opt, Xp, Pp, H, v, R, ih, ifrq,rcb->F, &rcb->TqLOM, nx, ny, ix, &rcb->sol.redundancy)) != -1) {
		matcpy(rcb->x, Xp, nx, 1);
		matcpy(rcb->P, Pp, nx, nx);
		stat = SOLQ_FLOAT;
		update_pre(rcb, rcb->x, opt, n, obs);
	}
	else
	{
		if (rcb->sol.redundancy < 0) {
			printf("++++++++++++ LSF filter failure,redundancy=%d ++++++++++++\n", rcb->sol.redundancy);
			printf("Possible poor data quality. The num of sat, GPS:%d GLO:%d GAL:%d QZSS:%d BDS:%d  Total:%d\n",
				 rcb->sol.sys_ns[0], rcb->sol.sys_ns[1], rcb->sol.sys_ns[2], rcb->sol.sys_ns[3], rcb->sol.sys_ns[4], rcb->sol.ns);

		}
		else {
			printf("++++++++++++ LSF filter failure,redundancy ++++++++++++\n");
		}
		free(Xp); free(Pp); free(v); free(H); free(R); free(ih); free(ifrq); free(ix);
		return 0;
	}
	/* out solution */
	if (stat == SOLQ_FLOAT) 
	{
		rcb->opt.ratio = 1;
		rcb->sol.stat = stat;
		OutPut(fp, rcb, rscs, obs, n, stat, sta, sopt);
		opt->Initial = opt->estm;
	}
	else {
		OutPut(fp, rcb, rscs, obs, n, stat, sta,sopt);
	}
	free(Xp); free(Pp); free(v); free(H); free(R); free(ih); free(ifrq); free(ix);
    return stat;
}
