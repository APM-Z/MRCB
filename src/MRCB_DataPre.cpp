#include "MRCB.h"        
#include "MRCB_SolHead.h"


/* constants/global variables */
static pcvs_t pcvss={0};        /* satellite antenna parameters */
static pcvr_t pcvsr={0};        /* receiver antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sta_t stas;              /* station infomation */
static int aborts=0;            /* abort status */

/* set antenna parameters */
static void setpcv(gtime_t time, const prcopt_t* popt, nav_t* nav, const pcvs_t* pcvs, const pcvr_t* pcvr, sta_t* sta)
{
    pcv_t* pcv;
    int i, j;
    char id[64];
    for (i = 0; i < MAXSAT; i++) {
        if (!(satsys(i + 1, NULL) & popt->navsys)) continue;
        if (!(pcv = searchpcv(i + 1, time, pcvs))) {
            satno2id(i + 1, id);
            continue;
        }
        nav->pcvs[i] = *pcv;
    }
    for (i = 0; i < 1; ++i) {
        for (j = 0; j < pcvr->n; ++j) {
            if (!strcmp(sta[i].antsno, pcvr->pcv[j].type)) {
                sta[i].pcvr = pcvr->pcv[j];
            }
        }
    }
}

/* free obs and nav data */
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}

/* input obs data */
static int inputobs(FILE*fp,int index,obs_t *obs,double *ver,int *tsys,char tobs[][NUMSYS][MAXOBSTYPE][4], const prcopt_t* popt, char* file)
{
	gtime_t ts={0},te={0};
	int n;
	double tint=0.0;
	obs->n=0;
	tint=0.0;
    readrnxobs(fp, ts, te, tint, NULL, 1, ver[index], &tsys[index], tobs[index], obs, NULL, popt);
	n=sortobs(obs);
	return n;
}
/* initialize control */
static void init(rcb_t *rcb, const prcopt_t *opt)
{
    sol_t  sol0={{0}};
    ssat_t ssat0={0};
    int i;
	rcb->opt=*opt;
    for (i=0;i<MAXSAT*MAXRCV;i++) {
        rcb->ssat[i]=ssat0;
    }
    
    for (i=0;i<NSYS;++i) {
        rcb->rsat[i]=-1;
        rcb->basis_freq[i]=0;
    }
    rcb->sol=sol0; 
    rcb->nx=rtknx(opt); 
    rcb->tt=0.0;
    rcb->F=(unsigned char *)calloc(sizeof(unsigned char),rcb->nx); memset(rcb->F,0,(unsigned char)rcb->nx);
    rcb->x=zeros(rcb->nx,1);
    rcb->P=zeros(rcb->nx,rcb->nx);
}
/* free rtk control */
static void confree(rcb_t *rcb)
{
    rcb->nx=0;
    free(rcb->F);  rcb->F =NULL;
    free(rcb->x ); rcb->x =NULL;
    free(rcb->P ); rcb->P =NULL;
}

/* read nav data */
static int readnavdata(char* infile, nav_t* nav)
{
    FILE* fp;
    char type = ' ';
    if (!(fp = fopen(infile, "r"))) {
        printf("open nav file error:%s", infile);
        return 0;
    }
        
    int sys, tsys = TSYS_GPS;
    double ver;
    char tobs[NUMSYS][MAXOBSTYPE][4] = { {""} };
    /* read rinex header */
    if (!readrnxh(fp, &ver, &type, &sys, &tsys, NULL, nav, NULL)) return 0;
    /* read rinex body */
    switch (type) {
    case 'N': return readrnxnav(fp, NULL, ver, sys, nav);
    case 'G': return readrnxnav(fp, NULL, ver, SYS_GLO, nav);
    case 'H': return readrnxnav(fp, NULL, ver, SYS_SBS, nav);
    case 'J': return readrnxnav(fp, NULL, ver, SYS_QZS, nav);
    case 'L': return readrnxnav(fp, NULL, ver, SYS_GAL, nav);
    default:return 0;
    }
}
/* process nav */
static int readnav(char* infile, nav_t* nav, const prcopt_t* popt)
{
    int stat;
    stat = readnavdata(infile, nav);
    if (nav->n <= 0 && nav->ng <= 0 && nav->ns <= 0) {
        return 0;
    }
    uniqnav(nav, popt);
    return stat;
}
/* process */
static void process(const prcopt_t *popt, const solopt_t *sopt, Ifile infile,int len, const filopt_t* fopt)
{
    
    char type=' ',tobs[60][NUMSYS][MAXOBSTYPE][4]={{{""}}};
    sol_t sol={0};rcb_t rcb;obs_t obs={0};
	char str[25];
	double ver[60];
    int i, nobs, n, prn, sys, tsys[60] = { TSYS_GPS };
    int k,judge=0;
    double ep[6] = { 0.0 };
    FILE* fp[3] = { NULL };
    int first_run = 1;
    char outfilename[256], profilename[256], logfilename[256];


    char sysstr[16] = "";
    // Build system string based on navigation systems
    if (popt->navsys & 1) strcat(sysstr, "G");   // GPS
    if (popt->navsys & 4) strcat(sysstr, "R");   // GLONASS  
    if (popt->navsys & 8) strcat(sysstr, "E");   // GALILEO
    if (popt->navsys & 16) strcat(sysstr, "J");  // QZSS
    if (popt->navsys & 32) strcat(sysstr, "C");  // BDS
    // Use default value if no systems are selected
    if (strlen(sysstr) == 0) {
        strcpy(sysstr, "Nosystem");
    }
    /* Create output files */
    if (popt->CONTINUE == 0) {
        sprintf(outfilename, "%s_%s_%d_%s_RCB", popt->staname, sysstr, popt->nf, popt->DOY);
        sprintf(profilename, "%s_%s_%d_%s_Inf", popt->staname, sysstr, popt->nf, popt->DOY);
        sprintf(logfilename, "%s_%s_%d_%s_Sat", popt->staname, sysstr, popt->nf, popt->DOY);
        if (!creatFile(fp, fopt, popt, sopt, outfilename, profilename, logfilename)) {
            return;
        }
    }
    else {
        sprintf(outfilename, "%s_%s_%d_Continues_RCB", popt->staname, sysstr, popt->nf);
        sprintf(profilename, "%s_%s_%d_Continues_Inf", popt->staname, sysstr, popt->nf);
        sprintf(logfilename, "%s_%s_%d_Continues_Sat", popt->staname, sysstr, popt->nf);
        if (!creatFile(fp, fopt, popt, sopt, outfilename, profilename, logfilename)) {
            return;
        }
    }


    init(&rcb, popt);
    for (k = 0; k < len; k++) {
        /* read nav */
        navs.eph = NULL; navs.n = navs.nmax = 0;
        navs.geph = NULL; navs.ng = navs.ngmax = 0;
        navs.seph = NULL; navs.ns = navs.nsmax = 0;
        if (!readnav(infile.rnxfile[k], &navs, popt)) {
            printf("read obs file error:%s\n", infile.rnxfile[k]);
            return;
        }
        /* read obs header */
        FILE* fpobs;
        if (!(obs.data = (obsd_t*)malloc(sizeof(obsd_t) * MAXOBS))) return;
        if (!(fpobs = fopen(infile.obsfile[k], "r"))) {
            printf("open obs file error:%s\n", infile.obsfile[k]);
            free(obs.data);
            return;
        }
        if (!readrnxh(fpobs, &ver[k], &type, &sys, &tsys[k], tobs[k], &navs, NULL)) {
            printf("read obs header error:%s\n", infile.obsfile[k]);
            fclose(fpobs);
            free(obs.data);
            return;
        }


        while ((nobs = inputobs(fpobs, k, &obs, ver, tsys, tobs, popt, infile.obsfile[k])) > 0) {
            for (i = n = 0; i < nobs; i++) {
                if (((sys = satsys(obs.data[i].sat, &prn)) & popt->navsys) &&
                    popt->exsats[obs.data[i].sat - 1] != 1) {
                    obs.data[n++] = obs.data[i];
                }
            }
            if (n <= 0) continue;
            if (!calrcb(fp, &rcb, obs.data, n, &navs, &stas, sopt)) {
                judge = 1;
                break;
            }
            time2str(rcb.sol.time, str, 0);
            if (rcb.TqLOM > 1E5) {
                judge = 1;
                printf("%s TqLOM=%.2f   Possible small number of satellites, ensure the target system has ¡Ý3. The num of sat, GPS:%d GLO:%d GAL:%d QZSS:%d BDS:%d  Total:%d\n",
                    str, rcb.TqLOM, rcb.sol.sys_ns[0], rcb.sol.sys_ns[1], rcb.sol.sys_ns[2], rcb.sol.sys_ns[3], rcb.sol.sys_ns[4], rcb.sol.ns);
                break;
            }

            if (rcb.sol.stat == SOLQ_FLOAT) {
                printf("%s %s ", str, "FLOAT");
                printf("%6.3f %6.3f %6.3f %8.3f %6.3f %8.3f %6.3f %6.3f %6.3f %6.3f %8.3f %6.3f %6.3f %8.3f %6.3f %6.3f %6.3f %6.3f %6d %10.3f\n",
                    rcb.sol.rcb[0][0], rcb.sol.rcb[0][1], rcb.sol.rcb[0][2], 
                    rcb.sol.rcb[1][0], rcb.sol.rcb[1][1],
                    rcb.sol.rcb[2][0], rcb.sol.rcb[2][1], rcb.sol.rcb[2][2], rcb.sol.rcb[2][3], rcb.sol.rcb[2][4],
                    rcb.sol.rcb[3][0], rcb.sol.rcb[3][1], rcb.sol.rcb[3][2],
                    rcb.sol.rcb[4][0], rcb.sol.rcb[4][1], rcb.sol.rcb[4][2], rcb.sol.rcb[4][3], rcb.sol.rcb[4][4],
                    rcb.sol.ns, rcb.TqLOM);
            }
            else if (rcb.sol.stat == SOLQ_SINGLE) {
                printf("%s %s ", str, "SINGLE");
                printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %6d %10.3f\n",rcb.sol.rcb[0][0], rcb.sol.rcb[0][1], rcb.sol.rcb[0][2], 
                    rcb.sol.rcb[1][0], rcb.sol.rcb[1][1],
                    rcb.sol.rcb[2][0], rcb.sol.rcb[2][1], rcb.sol.rcb[2][2], rcb.sol.rcb[2][3], rcb.sol.rcb[2][4],
                    rcb.sol.rcb[3][0], rcb.sol.rcb[3][1], rcb.sol.rcb[3][2],
                    rcb.sol.rcb[4][0], rcb.sol.rcb[4][1], rcb.sol.rcb[4][2], rcb.sol.rcb[4][3], rcb.sol.rcb[4][4],
                    rcb.sol.ns, rcb.TqLOM);
            }
        }
        if (!judge) {
            fclose(fpobs);
            free(obs.data);
        }
        else {
            break;
        }
    }

    for (i = 0; i < 3; i++) {
        if (fp[i] != NULL) {
            fclose(fp[i]);
        }
    }
    confree(&rcb);
}
/* execute processing session */
static int execses(prcopt_t *popt, const solopt_t *sopt, const filopt_t *fopt, int len, Ifile infile)
{
    if (popt->soltype == SOLTYPE_FORWARD) {
        process(popt, sopt, infile, len, fopt);
        freeobsnav(&obss, &navs); freeIfile(&infile);
    }
    return aborts?1:0;
}

extern int postpos(prcopt_t *popt, const solopt_t *sopt, filopt_t filopt)
{
    int stat=0;
    Ifile ifiles;
    if (!initIfile(&ifiles)) return 0;
    if (!getpath(&ifiles,popt,filopt)) return 0; 
    stat=execses(popt,sopt,&filopt,ifiles.idx[0],ifiles);
    printf("ending!\n");
    return 1;
}