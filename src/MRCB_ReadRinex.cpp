#include "MRCB.h"
/* constants/macros ----------------------------------------------------------*/
#define NUMSYS      8                   /* number of systems */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */
#define MINFREQ_GLO -7                  /* min frequency number glonass */
#define MAXFREQ_GLO 13                  /* max frequency number glonass */
#define NINCOBS     262144              /* inclimental number of obs data */

static const int navsys[]={             /* satellite systems */
    SYS_GPS,SYS_GLO,SYS_GAL,SYS_QZS,SYS_SBS,SYS_CMP,SYS_IRN,0
};
static const char syscodes[]="GREJSCI"; /* satellite system codes */

static const char obscodes[]="CLDS";    /* obs type codes */

static const char frqcodes[]="1256789"; /* frequency codes */

static const double ura_eph[]={         /* ura values (ref [3] 20.3.3.3.1.1) */
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};
static const double ura_nominal[]={     /* ura nominal values */
    2.0,2.8,4.0,5.7,8.0,11.3,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,
    2048.0,4096.0,8192.0
};
/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    unsigned char pri [MAXOBSTYPE];     /* signal priority (15-0) */
    unsigned char type[MAXOBSTYPE];     /* type (0:C,1:L,2:D,3:S) */
    unsigned char code[MAXOBSTYPE];     /* obs code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

/* set string without tail space ---------------------------------------------*/
static void setstr(char *dst, const char *src, int n)
{
    char *p=dst;
    const char *q=src;
    while (*q&&q<src+n) *p++=*q++;
    *p--='\0';
    while (p>=dst&&*p==' ') *p--='\0';
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjweek(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-302400.0) return timeadd(t, 604800.0);
    if (tt> 302400.0) return timeadd(t,-604800.0);
    return t;
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjday(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-43200.0) return timeadd(t, 86400.0);
    if (tt> 43200.0) return timeadd(t,-86400.0);
    return t;
}
/* satellite to satellite code -----------------------------------------------*/
static int sat2code(int sat, char *code)
{
    int prn;
    switch (satsys(sat,&prn)) {
        case SYS_GPS: sprintf(code,"G%2d",prn-MINPRNGPS+1); break;
        case SYS_GLO: sprintf(code,"R%2d",prn-MINPRNGLO+1); break;
        case SYS_GAL: sprintf(code,"E%2d",prn-MINPRNGAL+1); break;
        case SYS_SBS: sprintf(code,"S%2d",prn-100); break;
        case SYS_QZS: sprintf(code,"J%2d",prn-MINPRNQZS+1); break;
        case SYS_CMP: sprintf(code,"C%2d",prn-MINPRNCMP+1); break;
        case SYS_IRN: sprintf(code,"I%2d",prn-MINPRNIRN+1); break;
        default: return 0;
    }
    return 1;
}
/* ura index to ura nominal value (m) ----------------------------------------*/
static double uravalue(int sva)
{
    return 0<=sva&&sva<15?ura_nominal[sva]:8192.0;
}
/* ura value (m) to ura index ------------------------------------------------*/
static int uraindex(double value)
{
    int i;
    for (i=0;i<15;i++) if (ura_eph[i]>=value) break;
    return i;
}
/* galileo sisa index to sisa nominal value (m) ------------------------------*/
static double sisa_value(int sisa)
{
    if (sisa<= 49) return sisa*0.01;
    if (sisa<= 74) return 0.5+(sisa- 50)*0.02;
    if (sisa<= 99) return 1.0+(sisa- 75)*0.04;
    if (sisa<=125) return 2.0+(sisa-100)*0.16;
    return -1.0; /* unknown or NAPA */
}
/* galileo sisa value (m) to sisa index --------------------------------------*/
static int sisa_index(double value)
{
    if (value<0.0 || value>6.0) return 255; /* unknown or NAPA */
    else if (value<=0.5) return (int)(value/0.01);
    else if (value<=1.0) return (int)((value-0.5)/0.02)+50;
    else if (value<=2.0) return (int)((value-1.0)/0.04)+75;
    return (int)(((int)(value-2.0)/0.16)+100);
}
/* initialize station parameter ----------------------------------------------*/
static void init_sta(sta_t *sta)
{
    int i;
    *sta->name   ='\0';
    *sta->marker ='\0';
    *sta->antdes ='\0';
    *sta->antsno ='\0';
    *sta->rectype='\0';
    *sta->recver ='\0';
    *sta->recsno ='\0';
    sta->antsetup=sta->itrf=sta->deltype=0;
    for (i=0;i<3;i++) sta->pos[i]=0.0;
    for (i=0;i<3;i++) sta->del[i]=0.0;
    sta->hgt=0.0;
}
/*------------------------------------------------------------------------------
* input rinex functions
*-----------------------------------------------------------------------------*/

/* convert rinex obs type ver.2 -> ver.3 -------------------------------------*/
static void convcode(double ver, int sys, const char *str, char *type)
{
    strcpy(type,"   ");
    
    if      (!strcmp(str,"P1")) { /* ver.2.11 GPS L1PY,GLO L2P */
        if      (sys==SYS_GPS) sprintf(type,"%c1W",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c1P",'C');
    }
    else if (!strcmp(str,"P2")) { /* ver.2.11 GPS L2PY,GLO L2P */
        if      (sys==SYS_GPS) sprintf(type,"%c2W",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c2P",'C');
    }
    else if (!strcmp(str,"C1")) { /* ver.2.11 GPS L1C,GLO L1C/A */
        if      (ver>=2.12) ; /* reject C1 for 2.12 */
        else if (sys==SYS_GPS) sprintf(type,"%c1C",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c1C",'C');
        else if (sys==SYS_GAL) sprintf(type,"%c1X",'C'); /* ver.2.12 */
        else if (sys==SYS_QZS) sprintf(type,"%c1C",'C');
        else if (sys==SYS_SBS) sprintf(type,"%c1C",'C');
    }
    else if (!strcmp(str,"C2")) {
        if (sys==SYS_GPS) {
            if (ver>=2.12) sprintf(type,"%c2W",'C'); /* L2P(Y) */
            else           sprintf(type,"%c2X",'C'); /* L2C */
        }
        else if (sys==SYS_GLO) sprintf(type,"%c2C",'C');
        else if (sys==SYS_QZS) sprintf(type,"%c2X",'C');
        else if (sys==SYS_CMP) sprintf(type,"%c1X",'C'); /* ver.2.12 B1 */
    }
    else if (ver>=2.12&&str[1]=='A') { /* ver.2.12 L1C/A */
        if      (sys==SYS_GPS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_SBS) sprintf(type,"%c1C",str[0]);
    }
    else if (ver>=2.12&&str[1]=='B') { /* ver.2.12 GPS L1C */
        if      (sys==SYS_GPS) sprintf(type,"%c1X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c1X",str[0]);
    }
    else if (ver>=2.12&&str[1]=='C') { /* ver.2.12 GPS L2C */
        if      (sys==SYS_GPS) sprintf(type,"%c2X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c2X",str[0]);
    }
    else if (ver>=2.12&&str[1]=='D') { /* ver.2.12 GLO L2C/A */
        if      (sys==SYS_GLO) sprintf(type,"%c2C",str[0]);
    }
    else if (ver>=2.12&&str[1]=='1') { /* ver.2.12 GPS L1PY,GLO L1P */
        if      (sys==SYS_GPS) sprintf(type,"%c1W",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1P",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c1X",str[0]); /* tentative */
        else if (sys==SYS_CMP) sprintf(type,"%c1X",str[0]); /* extension */
    }
    else if (ver<2.12&&str[1]=='1') {
        if      (sys==SYS_GPS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c1X",str[0]); /* tentative */
        else if (sys==SYS_QZS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_SBS) sprintf(type,"%c1C",str[0]);
    }
    else if (str[1]=='2') {
        if      (sys==SYS_GPS) sprintf(type,"%c2W",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c2P",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c2X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c1X",str[0]); /* ver.2.12 B1 */
    }
    else if (str[1]=='5') {
        if      (sys==SYS_GPS) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_SBS) sprintf(type,"%c5X",str[0]);
    }
    else if (str[1]=='6') {
        if      (sys==SYS_GAL) sprintf(type,"%c6X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c6X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c6X",str[0]); /* ver.2.12 B3 */
    }
    else if (str[1]=='7') {
        if      (sys==SYS_GAL) sprintf(type,"%c7X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c7X",str[0]); /* ver.2.12 B2 */
    }
    else if (str[1]=='8') {
        if      (sys==SYS_GAL) sprintf(type,"%c8X",str[0]);
    }
    
}
/* decode obs header ---------------------------------------------------------*/
static void decode_obsh(FILE *fp, char *buff, const double ver, int *tsys,
                        char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta)
{
    /* default codes for unknown code */
    const char *defcodes[]={
        "CWX    ",  /* GPS: L125____ */
        "CC     ",  /* GLO: L12_____ */
        "X XXXX ",  /* GAL: L1_5678_ */
        "CXXX   ",  /* QZS: L1256___ */
        "C X    ",  /* SBS: L1_5____ */
        "X  XX  ",  /* BDS: L1__67__ */
        "  A   A"   /* IRN: L__5___9 */
    };
    double del[3];
    int i,j,k,n,nt,prn,fcn;
    const char *p;
    char *label=buff+60,str[4];
    

    
    if      (strstr(label,"MARKER NAME"         )) {
        if (sta) setstr(sta->name,buff,60);
    }
    else if (strstr(label,"MARKER NUMBER"       )) { /* opt */
        if (sta) setstr(sta->marker,buff,20);
    }
    else if (strstr(label,"MARKER TYPE"         )) ; /* ver.3 */
    else if (strstr(label,"OBSERVER / AGENCY"   )) ;
    else if (strstr(label,"REC # / TYPE / VERS" )) {
        if (sta) {
            setstr(sta->recsno, buff,   20);
            setstr(sta->rectype,buff+20,20);
            setstr(sta->recver, buff+40,20);
        }
    }
    else if (strstr(label,"ANT # / TYPE"        )) {
        if (sta) {
            setstr(sta->antsno,buff   ,20);
            setstr(sta->antdes,buff+20,20);
        }
    }
    else if (strstr(label,"APPROX POSITION XYZ" )) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) sta->pos[i]=str2num(buff,j,14);
        }
    }
    else if (strstr(label,"ANTENNA: DELTA H/E/N")) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) del[i]=str2num(buff,j,14);
            sta->del[2]=del[0]; /* h */
            sta->del[0]=del[1]; /* e */
            sta->del[1]=del[2]; /* n */
        }
    }
    else if (strstr(label,"ANTENNA: DELTA X/Y/Z")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: PHASECENTER")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: B.SIGHT XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR AZI")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"CENTER OF MASS: XYZ" )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / # / OBS TYPES" )) { /* ver.3 */
        if (!(p=strchr(syscodes,buff[0]))) {
            
            return;
        }
        i=(int)(p-syscodes);
        n=(int)str2num(buff,3,3);
        for (j=nt=0,k=7;j<n;j++,k+=4) {
            if (k>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                k=7;
            }
            if (nt<MAXOBSTYPE-1) setstr(tobs[i][nt++],buff+k,3);
        }
        *tobs[i][nt]='\0';
       
        /* if unknown code in ver.3, set default code */
        for (j=0;j<nt;j++) {
            if (tobs[i][j][2]) continue;
            if (!(p=strchr(frqcodes,tobs[i][j][1]))) continue;
            tobs[i][j][2]=defcodes[i][(int)(p-frqcodes)];
        }
    }
    else if (strstr(label,"WAVELENGTH FACT L1/2")) ; /* opt ver.2 */
    else if (strstr(label,"# / TYPES OF OBSERV" )) { /* ver.2 */
        n=(int)str2num(buff,0,6);
        for (i=nt=0,j=10;i<n;i++,j+=6) {
            if (j>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                j=10;
            }
            if (nt>=MAXOBSTYPE-1) continue;
            if (ver<=2.99) {
                setstr(str,buff+j,2);
                convcode(ver,SYS_GPS,str,tobs[0][nt]);
                convcode(ver,SYS_GLO,str,tobs[1][nt]);
                convcode(ver,SYS_GAL,str,tobs[2][nt]);
                convcode(ver,SYS_QZS,str,tobs[3][nt]);
                convcode(ver,SYS_SBS,str,tobs[4][nt]);
                convcode(ver,SYS_CMP,str,tobs[5][nt]);
            }
            nt++;
        }
        *tobs[0][nt]='\0';
    }
    else if (strstr(label,"SIGNAL STRENGTH UNIT")) ; /* opt ver.3 */
    else if (strstr(label,"INTERVAL"            )) ; /* opt */
    else if (strstr(label,"TIME OF FIRST OBS"   )) {
        if      (!strncmp(buff+48,"GPS",3)) *tsys=TSYS_GPS;
        else if (!strncmp(buff+48,"GLO",3)) *tsys=TSYS_UTC;
        else if (!strncmp(buff+48,"GAL",3)) *tsys=TSYS_GAL;
        else if (!strncmp(buff+48,"QZS",3)) *tsys=TSYS_QZS; /* ver.3.02 */
        else if (!strncmp(buff+48,"BDT",3)) *tsys=TSYS_CMP; /* ver.3.02 */
        else if (!strncmp(buff+48,"IRN",3)) *tsys=TSYS_IRN; /* ver.3.03 */
    }
    else if (strstr(label,"TIME OF LAST OBS"    )) ; /* opt */
    else if (strstr(label,"RCV CLOCK OFFS APPL" )) ; /* opt */
    else if (strstr(label,"SYS / DCBS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PCVS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / SCALE FACTOR"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PHASE SHIFTS"  )) ; /* ver.3.01 */
    else if (strstr(label,"GLONASS SLOT / FRQ #")) { /* ver.3.02 */
        if (nav) {
            for (i=0,p=buff+4;i<8;i++,p+=8) {
                if (sscanf(p,"R%2d %2d",&prn,&fcn)<2) continue;
                if (1<=prn&&prn<=MAXPRNGLO) nav->glo_fcn[prn-1]=fcn+8;
            }
        }
    }
    else if (strstr(label,"GLONASS COD/PHS/BIS" )) { /* ver.3.02 */
        if (nav) {
            for (i=0,p=buff;i<4;i++,p+=13) {
                if      (strncmp(p+1,"C1C",3)) nav->glo_cpbias[0]=str2num(p,5,8);
                else if (strncmp(p+1,"C1P",3)) nav->glo_cpbias[1]=str2num(p,5,8);
                else if (strncmp(p+1,"C2C",3)) nav->glo_cpbias[2]=str2num(p,5,8);
                else if (strncmp(p+1,"C2P",3)) nav->glo_cpbias[3]=str2num(p,5,8);
            }
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
    else if (strstr(label,"# OF SALTELLITES"    )) { /* opt */
        /* skip */ ;
    }
    else if (strstr(label,"PRN / # OF OBS"      )) { /* opt */
        /* skip */ ;
    }
}
/* decode nav header ---------------------------------------------------------*/
static void decode_navh(char *buff, nav_t *nav)
{
    int i,j;
    char *label=buff+60;
    
    
    if      (strstr(label,"ION ALPHA"           )) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=2;i<4;i++,j+=12) nav->ion_gps[i]=str2num(buff,j,12);
        }
    }
    else if (strstr(label,"ION BETA"            )) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=2;i<4;i++,j+=12) nav->ion_gps[i+4]=str2num(buff,j,12);
        }
    }
    else if (strstr(label,"DELTA-UTC: A0,A1,T,W")) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=3;i<2;i++,j+=19) nav->utc_gps[i]=str2num(buff,j,19);
            for (;i<4;i++,j+=9) nav->utc_gps[i]=str2num(buff,j,9);
        }
    }
    else if (strstr(label,"IONOSPHERIC CORR"    )) { /* opt ver.3 */
        if (nav) {
            if (!strncmp(buff,"GPSA",4)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gps[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"GPSB",4)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gps[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"GAL",3)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gal[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"QZSA",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_qzs[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"QZSB",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_qzs[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"BDSA",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_cmp[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"BDSB",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_cmp[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"IRNA",4)) { /* v.3.03 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_irn[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"IRNB",4)) { /* v.3.03 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_irn[i+4]=str2num(buff,j,12);
            }
        }
    }
    else if (strstr(label,"TIME SYSTEM CORR"    )) { /* opt ver.3 */
        if (nav) {
            if (!strncmp(buff,"GPUT",4)) {
                nav->utc_gps[0]=str2num(buff, 5,17);
                nav->utc_gps[1]=str2num(buff,22,16);
                nav->utc_gps[2]=str2num(buff,38, 7);
                nav->utc_gps[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"GLUT",4)) {
                nav->utc_glo[0]=str2num(buff, 5,17);
                nav->utc_glo[1]=str2num(buff,22,16);
            }
            else if (!strncmp(buff,"GAUT",4)) { /* v.3.02 */
                nav->utc_gal[0]=str2num(buff, 5,17);
                nav->utc_gal[1]=str2num(buff,22,16);
                nav->utc_gal[2]=str2num(buff,38, 7);
                nav->utc_gal[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"QZUT",4)) { /* v.3.02 */
                nav->utc_qzs[0]=str2num(buff, 5,17);
                nav->utc_qzs[1]=str2num(buff,22,16);
                nav->utc_qzs[2]=str2num(buff,38, 7);
                nav->utc_qzs[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"BDUT",4)) { /* v.3.02 */
                nav->utc_cmp[0]=str2num(buff, 5,17);
                nav->utc_cmp[1]=str2num(buff,22,16);
                nav->utc_cmp[2]=str2num(buff,38, 7);
                nav->utc_cmp[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"SBUT",4)) { /* v.3.02 */
                nav->utc_cmp[0]=str2num(buff, 5,17);
                nav->utc_cmp[1]=str2num(buff,22,16);
                nav->utc_cmp[2]=str2num(buff,38, 7);
                nav->utc_cmp[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"IRUT",4)) { /* v.3.03 */
                nav->utc_irn[0]=str2num(buff, 5,17);
                nav->utc_irn[1]=str2num(buff,22,16);
                nav->utc_irn[2]=str2num(buff,38, 7);
                nav->utc_irn[3]=str2num(buff,45, 5);
            }
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
}
/* decode gnav header --------------------------------------------------------*/
static void decode_gnavh(char *buff, nav_t *nav)
{
    char *label=buff+60;
    if      (strstr(label,"CORR TO SYTEM TIME"  )) ; /* opt */
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
}
/* decode geo nav header -----------------------------------------------------*/
static void decode_hnavh(char *buff, nav_t *nav)
{
    char *label=buff+60;
    
    
    
    if      (strstr(label,"CORR TO SYTEM TIME"  )) ; /* opt */
    else if (strstr(label,"D-UTC A0,A1,T,W,S,U" )) ; /* opt */
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
}
/* read rinex header ---------------------------------------------------------*/
extern int readrnxh(FILE *fp, double *ver, char *type, int *sys, int *tsys,
                    char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta)
{
    char buff[MAXRNXLEN], * label = buff + 60;
    int i = 0, block = 0;


    *ver = 2.10; *type = ' '; *sys = SYS_GPS; *tsys = TSYS_GPS;


    while (fgets(buff, MAXRNXLEN, fp)) {
        if (strlen(buff) <= 60) continue;
        else if (strstr(label, "RINEX VERSION / TYPE")) {
            *ver = str2num(buff, 0, 9);
            *type = (*ver < 3.04) ? *(buff + 20) : *(buff + 21);
            if (*type != 'C') {
                *type = *(buff + 20);
            }
            /* satellite system */
            switch ((*ver < 3.04) ? *(buff + 40) : *(buff + 42)) {
            case ' ':
            case 'G': *sys = SYS_GPS;  *tsys = TSYS_GPS; break;
            case 'R': *sys = SYS_GLO;  *tsys = TSYS_UTC; break;
            case 'E': *sys = SYS_GAL;  *tsys = TSYS_GAL; break; /* v.2.12 */
            case 'S': *sys = SYS_SBS;  *tsys = TSYS_GPS; break;
            case 'J': *sys = SYS_QZS;  *tsys = TSYS_QZS; break; /* v.3.02 */
            case 'C': *sys = SYS_CMP;  *tsys = TSYS_CMP; break; /* v.2.12 */
            case 'I': *sys = SYS_IRN;  *tsys = TSYS_IRN; break; /* v.3.03 */
            case 'M': *sys = SYS_NONE; *tsys = TSYS_GPS; break; /* mixed */
            default: break;
            }
            continue;
        }
        else if (strstr(label, "PGM / RUN BY / DATE")) continue;
        else if (strstr(label, "COMMENT")) { /* opt */
            continue;
        }
        /* file type */
        switch (*type) {
            case 'O': decode_obsh(fp, buff, *ver, tsys, tobs, nav, sta); break;
            case 'N': decode_navh(buff, nav); break;
            case 'G': decode_gnavh(buff, nav); break;
            case 'H': decode_hnavh(buff, nav); break;
            case 'J': decode_navh(buff, nav); break; /* extension */
            case 'L': decode_navh(buff, nav); break; /* extension */
        }
        if (strstr(label, "END OF HEADER")) return 1;
        if (++i >= MAXPOSHEAD && *type == ' ') break; /* no rinex file */
    }
    return 0;
}

/* decode obs epoch ----------------------------------------------------------*/
static int decode_obsepoch(FILE *fp, char *buff, double ver, gtime_t *time,
                           int *flag, int *sats)
{
    int i,j,n;
    char satid[8]="";
    
    
    
    if (ver<=2.99) { /* ver.2 */
        if ((n=(int)str2num(buff,29,3))<=0) return 0;
        
        /* epoch flag: 3:new site,4:header info,5:external event */
        *flag=(int)str2num(buff,28,1);
        
        if (3<=*flag&&*flag<=5) return n;
        
        if (str2time(buff,0,26,time)) {
            return 0;
        }
        for (i=0,j=32;i<n;i++,j+=3) {
            if (j>=68) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                j=32;
            }
            if (i<MAXOBS) {
                strncpy(satid,buff+j,3);
                sats[i]=satid2no(satid);
            }
        }
    }
    else { /* ver.3 */
        if ((n=(int)str2num(buff,32,3))<=0) return 0;
        
        *flag=(int)str2num(buff,31,1);
        
        if (3<=*flag&&*flag<=5) return n;
        
        if (buff[0]!='>'||str2time(buff,1,28,time)) {
            return 0;
        }
    }
    return n;
}
/* decode obs data -----------------------------------------------------------*/
static int decode_obsdata(FILE *fp, char *buff, double ver, int mask,
                          sigind_t *index, obsd_t *obs)
{
    sigind_t *ind;
    double val[MAXOBSTYPE]={0},sumval=0.0;
    unsigned char lli[MAXOBSTYPE]={0};
    char satid[8]="";
    int i,j,n,m,stat=1,p[MAXOBSTYPE],k[16],l[16],pntype=0,lntype=0,pri[NFREQ]={0};
    
    
    if (ver>2.99) { /* ver.3 */
        strncpy(satid,buff,3);
        obs->sat=(unsigned char)satid2no(satid);
    }
    if (!obs->sat) {
        stat=0;
    }
    else if (!(satsys(obs->sat,NULL)&mask)) {
        stat=0;
    }
    /* read obs data fields */
    switch (satsys(obs->sat,NULL)) {
        case SYS_GLO: ind=index+1; break;
        case SYS_GAL: ind=index+2; break;
        case SYS_QZS: ind=index+3; break;
        case SYS_SBS: ind=index+4; break;
        case SYS_CMP: ind=index+5; break;
        default:      ind=index  ; break;
    }
    for (i=0,j=ver<=2.99?0:3;i<ind->n;i++,j+=16) {
        
        if (ver<=2.99&&j>=80) { /* ver.2 */
            if (!fgets(buff,MAXRNXLEN,fp)) break;
            j=0;
        }
        if (stat) {
            val[i]=str2num(buff,j,14)+ind->shift[i];
            lli[i]=(unsigned char)str2num(buff,j+14,1)&3;
        }
    }
    if (!stat) return 0;
    
    for (i=0;i<NFREQ+NEXOBS;i++) {
        obs->P[i]=obs->L[i]=0.0; obs->D[i]=0.0f;
        obs->SNR[i]=obs->LLI[i]=obs->code[i]=0;
        for (j=0;j<TYPEN;j++) obs->EP[i][j]=obs->EL[i][j]=0.0f;
    }
    /* assign position in obs data */
    for (i=n=m=0;i<ind->n;i++) {
        
        p[i]=ver<=2.11?ind->frq[i]-1:ind->pos[i];
        
        if (ind->type[i]==0&&p[i]==0) k[n++]=i; /* C1? index */
        if (ind->type[i]==0&&p[i]==1) l[m++]=i; /* C2? index */
    }
    if (ver<=2.11) {
        
        /* if multiple codes (C1/P1,C2/P2), select higher priority */
        if (n>=2) {
            if (val[k[0]]==0.0&&val[k[1]]==0.0) {
                p[k[0]]=-1; p[k[1]]=-1;
            }
            else if (val[k[0]]!=0.0&&val[k[1]]==0.0) {
                p[k[0]]=0; p[k[1]]=-1;
            }
            else if (val[k[0]]==0.0&&val[k[1]]!=0.0) {
                p[k[0]]=-1; p[k[1]]=0;
            }
            else if (ind->pri[k[1]]>ind->pri[k[0]]) {
                p[k[1]]=0; p[k[0]]=NEXOBS<1?-1:NFREQ;
            }
            else {
                p[k[0]]=0; p[k[1]]=NEXOBS<1?-1:NFREQ;
            }
        }
        if (m>=2) {
            if (val[l[0]]==0.0&&val[l[1]]==0.0) {
                p[l[0]]=-1; p[l[1]]=-1;
            }
            else if (val[l[0]]!=0.0&&val[l[1]]==0.0) {
                p[l[0]]=1; p[l[1]]=-1;
            }
            else if (val[l[0]]==0.0&&val[l[1]]!=0.0) {
                p[l[0]]=-1; p[l[1]]=1; 
            }
            else if (ind->pri[l[1]]>ind->pri[l[0]]) {
                p[l[1]]=1; p[l[0]]=NEXOBS<2?-1:NFREQ+1;
            }
            else {
                p[l[0]]=1; p[l[1]]=NEXOBS<2?-1:NFREQ+1;
            }
        }
    }
    for (i=0;i<ind->n;++i){
        if(p[i]==-1) val[i]=0.0;
    }
    /* save obs data */
    while(1){
        for (i=0;i<ind->n;i++) {
            if (val[i]==0.0||p[i]>=NFREQ) continue;
            switch (ind->type[i]) {
                case 0: obs->P[p[i]]=val[i]; obs->code[p[i]]=ind->code[i]; break;
                case 1: obs->L[p[i]]=val[i]; obs->LLI [p[i]]=lli[i];       break;
                case 2: obs->D[p[i]]=(float)val[i];                        break;
                case 3: obs->SNR[p[i]]=(unsigned char)(val[i]*4.0+0.5);    break;
            }
        }
        for(i=0;i<NFREQ;++i){
            if(obs->P[i]!=0&&obs->L[i]!=0){
                for(j=0;j<ind->n;++j){
                    if(p[j]%NFREQ==i) {
                        p[j]=-1;
                        val[j]=0.0;
                    }
                }
            }
            else{
                for(j=0;j<ind->n;++j){
                    if(p[j]==i) {
                        p[j]=-1;
                        val[j]=0.0;
                    }
                }
                pri[i]+=1;
                for(j=0;j<ind->n;++j){
                    if(p[j]==pri[i]*NFREQ+i) 
                        p[j]=i;
                }
            }
        }
        sumval=0.0;
        for(i=0;i<ind->n;++i) sumval+=val[i];
        if (sumval==0.0) break;
    }
    return 1;
}
/* save slips ----------------------------------------------------------------*/
static void saveslips(unsigned char slips[][NFREQ], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ;i++) {
        if (data->LLI[i]&1) slips[data->sat-1][i]|=LLI_SLIP;
    }
}
/* restore slips -------------------------------------------------------------*/
static void restslips(unsigned char slips[][NFREQ], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ;i++) {
        if (slips[data->sat-1][i]&1) data->LLI[i]|=LLI_SLIP;
        slips[data->sat-1][i]=0;
    }
}
/* add obs data --------------------------------------------------------------*/
static int addobsdata(obs_t *obs, const obsd_t *data)
{
    obsd_t *obs_data;
    
    if (obs->nmax<=obs->n) {
        if (obs->nmax<=0) obs->nmax=NINCOBS; else obs->nmax*=2;
        if (!(obs_data=(obsd_t *)realloc(obs->data,sizeof(obsd_t)*obs->nmax))) {
            free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
            return -1;
        }
        obs->data=obs_data;
    }
    obs->data[obs->n++]=*data;
    return 1;
}
/* set system mask -----------------------------------------------------------*/
static int set_sysmask(const char *opt)
{
    const char *p;
    int mask=SYS_NONE;
    
    if (!opt) return SYS_ALL;
    if (!(p=strstr(opt,"-SYS="))) return SYS_ALL;
    
    for (p+=5;*p&&*p!=' ';p++) {
        switch (*p) {
            case 'G': mask|=SYS_GPS; break;
            case 'R': mask|=SYS_GLO; break;
            case 'E': mask|=SYS_GAL; break;
            case 'J': mask|=SYS_QZS; break;
            case 'C': mask|=SYS_CMP; break;
            case 'I': mask|=SYS_IRN; break;
            case 'S': mask|=SYS_SBS; break;
        }
    }
    return mask;
}

static void set_freq(sigind_t* ind, const prcopt_t* popt,int sys,int i,double ver) 
{
    if (popt->bds2 == 0 || popt->bds2 == 2) {
        if (sys == SYS_CMP) {
            if (ind->frq[i] == 2) ind->frq[i] = 1; /* B1-2       */
            else if (ind->frq[i] == 7) ind->frq[i] = 3; /* B2b(BDS-3)         */
            else if (ind->frq[i] == 4) ind->frq[i] = 2; /* B3 */
            else if (ind->frq[i] == 1) ind->frq[i] = 4; /* B1 */
            else if (ind->frq[i] == 3) ind->frq[i] = 5; /* B2a(BDS-3) */
            else if (ind->frq[i] == 5) ind->frq[i] = 6; /* B2b(BDS-2) */
        }
        else if (sys == SYS_GAL) {
            if (ind->frq[i] == 3) ind->frq[i] = 2; /* E5a     */
            else if (ind->frq[i] == 5) ind->frq[i] = 3; /* E5b     */
            else if (ind->frq[i] == 6) ind->frq[i] = 5; /* E5      */
        }
    }
    else if (popt->bds2 == 1) {
        if (sys == SYS_CMP) {
            if (ind->frq[i] == 2) ind->frq[i] = 1; /* B1-2       */
            else if (ind->frq[i] == 5) ind->frq[i] = 3; /* B2b(BDS-2)       */
            else if (ind->frq[i] == 4) ind->frq[i] = 2; /* B3 */
            else if (ind->frq[i] == 3) ind->frq[i] = 5; /* B2a(BDS-3) */
            else if (ind->frq[i] == 7) ind->frq[i] = 6; /* B2b(BDS-3) */
            else if (ind->frq[i] == 1 && ver != 3.02) ind->frq[i] = 4; /* B1 (BDS-3) */
        }
        else if (sys == SYS_GAL) {
            if (ind->frq[i] == 3) ind->frq[i] = 2; /* E5a     */
            else if (ind->frq[i] == 5) ind->frq[i] = 3; /* E5b     */
            else if (ind->frq[i] == 6) ind->frq[i] = 5; /* E5      */
        }
    }
}



/* set signal index ----------------------------------------------------------*/
static void set_index(double ver, int sys, const char *opt,
                      char tobs[MAXOBSTYPE][4], sigind_t *ind, const prcopt_t* popt)
{
    const char *p;
    char str[8];
    const char *optstr = "";
    const char *tt = "";      
    double shift;
    int i,j,k,n,pri=0,freq_num=0;
    for (i=n=0;*tobs[i];i++,n++) {
        tt = *tobs;
        ind->code[i]=obs2code(tobs[i]+1,ind->frq+i);
        ind->type[i]=(p=strchr(obscodes,tobs[i][0]))?(int)(p-obscodes):255;
        ind->pri[i]=getcodepri(sys,ind->code[i],opt);
        ind->pos[i]=-1;
        set_freq(ind,popt,sys,i,ver);
    }
    /* parse phase shift options */
    switch (sys) {
        case SYS_GPS: optstr="-GL%2s=%lf"; break;
        case SYS_GLO: optstr="-RL%2s=%lf"; break;
        case SYS_GAL: optstr="-EL%2s=%lf"; break;
        case SYS_QZS: optstr="-JL%2s=%lf"; break;
        case SYS_SBS: optstr="-SL%2s=%lf"; break;
        case SYS_CMP: optstr="-CL%2s=%lf"; break;
        case SYS_IRN: optstr="-IL%2s=%lf"; break;
    }
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str,&shift)<2) continue;
        for (i=0;i<n;i++) {
            if (strcmp(code2obs(ind->code[i],NULL),str)) continue;
            ind->shift[i]=shift;
        }
    }
    /* assign index for highest priority code */
    while(1){
        freq_num=0;
        for (i=0;i<NFREQ;i++) {
            for (j=0,k=-1;j<n;j++) {
                if (ind->frq[j]==i+1&&ind->pri[j]&&(k<0||ind->pri[j]>ind->pri[k])) {
                    k=j;
                }
            }
            if (k<0) {
                freq_num++;
                if(freq_num>=NFREQ) break;
                continue;
            }
            for (j=0;j<n;j++) {
                if (ind->code[j]==ind->code[k]) {
                    ind->pos[j]=i+NFREQ*pri;
                    ind->pri[j]=0;
                }
            }
        }
        for (i=0;i<n;i++) {
            if (!ind->code[i]||!ind->pri[i]||ind->pos[i]>=0) continue;
        }
        if(freq_num>=NFREQ) break;
        pri+=1;
    }
    ind->n=n;
}
/* read rinex obs data body --------------------------------------------------*/
static int readrnxobsb(FILE *fp, const char *opt, double ver, int *tsys,
                       char tobs[][MAXOBSTYPE][4], int *flag, obsd_t *data,
                       sta_t *sta, const prcopt_t* popt)
{
    gtime_t time={0};
    sigind_t index[7]={{0}};
    char buff[MAXRNXLEN];
    int i=0,n=0,nsat=0,sats[MAXOBS]={0},mask;
    
    /* set system mask */
    mask=set_sysmask(opt);
    
    /* set signal index */
    set_index(ver,SYS_GPS,opt,tobs[0],index  , popt);
    set_index(ver,SYS_GLO,opt,tobs[1],index+1, popt);
    set_index(ver,SYS_GAL,opt,tobs[2],index+2, popt);
    set_index(ver,SYS_QZS,opt,tobs[3],index+3, popt);
    set_index(ver,SYS_SBS,opt,tobs[4],index+4, popt);
    set_index(ver,SYS_CMP,opt,tobs[5],index+5, popt);
    set_index(ver,SYS_IRN,opt,tobs[6],index+6, popt);
    
    /* read record */
    while (fgets(buff,MAXRNXLEN,fp)) {
        
        /* decode obs epoch */
        if (i==0) {
            if ((nsat=decode_obsepoch(fp,buff,ver,&time,flag,sats))<=0) {
                continue;
            }
        }
        else if (*flag<=2||*flag==6) {
            
            data[n].time=time;
            data[n].sat=(unsigned char)sats[i-1];
            
            /* decode obs data */
            if (decode_obsdata(fp,buff,ver,mask,index,data+n)&&n<MAXOBS) n++;
        }
        else if (*flag==3||*flag==4) { 
            /* decode obs header */
            decode_obsh(fp,buff,ver,tsys,tobs,NULL,sta);
        }
        if (++i>nsat) return n;
    }
    return -1;
}
/* read rinex obs ------------------------------------------------------------*/
extern int readrnxobs(FILE *fp, gtime_t ts, gtime_t te, double tint,
                      const char *opt, int rcv, double ver, int *tsys,
                      char tobs[][MAXOBSTYPE][4], obs_t *obs, sta_t *sta, const prcopt_t* popt)
{
 
    obsd_t *data;
    unsigned char slips[MAXSAT][NFREQ]={{0}};
    int i,n,flag=0;

    if (!obs||rcv>MAXRCV) return 0;
    
    if (!(data = (obsd_t*)malloc(sizeof(obsd_t) * MAXOBS))) {
        return 0;
    }
        
    
    /* read rinex obs data body */
    if ((n = readrnxobsb(fp, opt, ver, tsys, tobs, &flag, data, sta, popt)) <= 0) {
        return 0;
    }
        

    for (i=0;i<n;i++) {
        /* utc -> gpst */
        if (*tsys==TSYS_UTC) data[i].time=utc2gpst(data[i].time);
        /* save cycle-slip */
        saveslips(slips,data+i);
    }
    /* screen data by time */
    if (n > 0 && !screent(data[0].time, ts, te, tint)) {
        return 0;
    }
        
        
    for (i=0;i<n;i++) {
        /* restore cycle-slip */
        restslips(slips,data+i);
        data[i].rcv=(unsigned char)rcv;
        /* save obs data */
	    obs->data[obs->n++]=data[i];
    }
    free(data);
    return 1;
}
/* decode ephemeris ----------------------------------------------------------*/
static int decode_eph(double ver, int sat, gtime_t toc, const double *data,
                      eph_t *eph)
{
    eph_t eph0={0};
    int sys;
    
    sys=satsys(sat,NULL);
    
    if (!(sys&(SYS_GPS|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN))) {
        return 0;
    }
    *eph=eph0;
    
    eph->sat=sat;
    eph->toc=toc;
    
    eph->f0=data[0];
    eph->f1=data[1];
    eph->f2=data[2];
    
    eph->A=SQR(data[10]); eph->e=data[ 8]; eph->i0  =data[15]; eph->OMG0=data[13];
    eph->omg =data[17]; eph->M0 =data[ 6]; eph->deln=data[ 5]; eph->OMGd=data[18];
    eph->idot=data[19]; eph->crc=data[16]; eph->crs =data[ 4]; eph->cuc =data[ 7];
    eph->cus =data[ 9]; eph->cic=data[12]; eph->cis =data[14];
    
    if (sys==SYS_GPS||sys==SYS_QZS) {
        eph->iode=(int)data[ 3];      /* IODE */
        eph->iodc=(int)data[26];      /* IODC */
        eph->toes=     data[11];      /* toe (s) in gps week */
        eph->week=(int)data[21];      /* gps week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);
        
        eph->code=(int)data[20];      /* GPS: codes on L2 ch */
        eph->svh =(int)data[24];      /* sv health */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */
        eph->flag=(int)data[22];      /* GPS: L2 P data flag */
        
        eph->tgd[0]=   data[25];      /* TGD */
        if (sys==SYS_GPS) {
            eph->fit=data[28];        /* fit interval (h) */
        }
        else {
            eph->fit=data[28]==0.0?1.0:2.0; /* fit interval (0:1h,1:>2h) */
        }
    }
    else if (sys==SYS_GAL) { /* GAL ver.3 */
        eph->iode=(int)data[ 3];      /* IODnav */
        eph->toes=     data[11];      /* toe (s) in galileo week */
        eph->week=(int)data[21];      /* gal week = gps week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);
        
        eph->code=(int)data[20];      /* data sources */
                                      /* bit 0 set: I/NAV E1-B */
                                      /* bit 1 set: F/NAV E5a-I */
                                      /* bit 2 set: F/NAV E5b-I */
                                      /* bit 8 set: af0-af2 toc are for E5a.E1 */
                                      /* bit 9 set: af0-af2 toc are for E5b.E1 */
        eph->svh =(int)data[24];      /* sv health */
                                      /* bit     0: E1B DVS */
                                      /* bit   1-2: E1B HS */
                                      /* bit     3: E5a DVS */
                                      /* bit   4-5: E5a HS */
                                      /* bit     6: E5b DVS */
                                      /* bit   7-8: E5b HS */
        eph->sva =sisa_index(data[23]); /* sisa (m->index) */
        
        eph->tgd[0]=   data[25];      /* BGD E5a/E1 */
        eph->tgd[1]=   data[26];      /* BGD E5b/E1 */
    }
    else if (sys==SYS_CMP) { /* BeiDou v.3.02 */
        eph->toc=bdt2gpst(eph->toc);  /* bdt -> gpst */
        eph->iode=(int)data[ 3];      /* AODE */
        eph->iodc=(int)data[28];      /* AODC */
        eph->toes=     data[11];      /* toe (s) in bdt week */
        eph->week=(int)data[21];      /* bdt week */
        eph->toe=bdt2gpst(bdt2time(eph->week,data[11])); /* bdt -> gpst */
        eph->ttr=bdt2gpst(bdt2time(eph->week,data[27])); /* bdt -> gpst */
        eph->toe=adjweek(eph->toe,toc);
        eph->ttr=adjweek(eph->ttr,toc);
        
        eph->svh =(int)data[24];      /* satH1 */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */
        
        eph->tgd[0]=   data[25];      /* TGD1 B1/B3 */
        eph->tgd[1]=   data[26];      /* TGD2 B2/B3 */

        eph->iode=(int)((int)(eph->toes/720)%240);
    }
    else if (sys==SYS_IRN) { /* IRNSS v.3.03 */
        eph->iode=(int)data[ 3];      /* IODEC */
        eph->toes=     data[11];      /* toe (s) in irnss week */
        eph->week=(int)data[21];      /* irnss week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);
        eph->svh =(int)data[24];      /* sv health */
        eph->sva=uraindex(data[23]);  /* ura (m->index) */
        eph->tgd[0]=   data[25];      /* TGD */
    }
    if (eph->iode<0||1023<eph->iode) {}
    if (eph->iodc<0||1023<eph->iodc) {}
    return 1;
}
/* decode glonass ephemeris --------------------------------------------------*/
static int decode_geph(double ver, int sat, gtime_t toc, double *data, geph_t *geph)
{
    geph_t geph0={0};
    gtime_t tof;
    double tow,tod;
    int week,dow;
    
    
    if (satsys(sat,NULL)!=SYS_GLO) {
        return 0;
    }
    *geph=geph0;
    
    geph->sat=sat;
    
    /* toc rounded by 15 min in utc */
    tow=time2gpst(toc,&week);
    toc=gpst2time(week,floor((tow+450.0)/900.0)*900);
    dow=(int)floor(tow/86400.0);
    
    /* time of frame in utc */
    tod=ver<=2.99?data[2]:fmod(data[2],86400.0); /* tod (v.2), tow (v.3) in utc */
    tof=gpst2time(week,tod+dow*86400.0);
    tof=adjday(tof,toc);
    
    geph->toe=utc2gpst(toc);   /* toc (gpst) */
    geph->tof=utc2gpst(tof);   /* tof (gpst) */
    
    /* iode = tb (7bit), tb =index of UTC+3H within current day */
    geph->iode=(int)(fmod(tow+10800.0,86400.0)/900.0+0.5);
    
    geph->taun=-data[0];       /* -taun */
    geph->gamn= data[1];       /* +gamman */
    
    geph->pos[0]=data[3]*1E3; geph->pos[1]=data[7]*1E3; geph->pos[2]=data[11]*1E3;
    geph->vel[0]=data[4]*1E3; geph->vel[1]=data[8]*1E3; geph->vel[2]=data[12]*1E3;
    geph->acc[0]=data[5]*1E3; geph->acc[1]=data[9]*1E3; geph->acc[2]=data[13]*1E3;
    
    geph->svh=(int)data[ 6];
    geph->frq=(int)data[10];
    geph->age=(int)data[14];
    
    /* some receiver output >128 for minus frequency number */
    if (geph->frq>128) geph->frq-=256;
    
    if (geph->frq<MINFREQ_GLO||MAXFREQ_GLO<geph->frq) {}
    return 1;
}
/* decode geo ephemeris ------------------------------------------------------*/
static int decode_seph(double ver, int sat, gtime_t toc, double *data,
                       seph_t *seph)
{
    seph_t seph0={0};
    int week;
    
    
    if (satsys(sat,NULL)!=SYS_SBS) {
        return 0;
    }
    *seph=seph0;
    
    seph->sat=sat;
    seph->t0 =toc;
    
    time2gpst(toc,&week);
    seph->tof=adjweek(gpst2time(week,data[2]),toc);
    
    seph->af0=data[0];
    seph->af1=data[1];
    
    seph->pos[0]=data[3]*1E3; seph->pos[1]=data[7]*1E3; seph->pos[2]=data[11]*1E3;
    seph->vel[0]=data[4]*1E3; seph->vel[1]=data[8]*1E3; seph->vel[2]=data[12]*1E3;
    seph->acc[0]=data[5]*1E3; seph->acc[1]=data[9]*1E3; seph->acc[2]=data[13]*1E3;
    
    seph->svh=(int)data[6];
    seph->sva=uraindex(data[10]);
    
    return 1;
}
/* read rinex navigation data body -------------------------------------------*/
static int readrnxnavb(FILE *fp, const char *opt, double ver, int sys,
                       int *type, eph_t *eph, geph_t *geph, seph_t *seph)
{
    gtime_t toc;
    double data[64];
    int i=0,j,prn,sat=0,sp=3,mask;
    char buff[MAXRNXLEN],id[8]="",*p;
    
    
    /* set system mask */
    mask=set_sysmask(opt);
    
    while (fgets(buff,MAXRNXLEN,fp)) {
        
        if (i==0) {
            
            /* decode satellite field */
            if (ver>=3.0||sys==SYS_GAL||sys==SYS_QZS) { /* ver.3 or GAL/QZS */
                strncpy(id,buff,3);
                sat=satid2no(id);
                sp=4;
                if (ver>=3.0) sys=satsys(sat,NULL);
            }
            else {
                prn=(int)str2num(buff,0,2);
                
                if (sys==SYS_SBS) {
                    sat=satno(SYS_SBS,prn+100);
                }
                else if (sys==SYS_GLO) {
                    sat=satno(SYS_GLO,prn);
                }
                else if (93<=prn&&prn<=97) { /* extension */
                    sat=satno(SYS_QZS,prn+100);
                }
                else sat=satno(SYS_GPS,prn);
            }
            /* decode toc field */
            if (str2time(buff+sp,0,19,&toc)) {
                return 0;
            }
            /* decode data fields */
            for (j=0,p=buff+sp+19;j<3;j++,p+=19) {
                data[i++]=str2num(p,0,19);
            }
        }
        else {
            /* decode data fields */
            for (j=0,p=buff+sp;j<4;j++,p+=19) {
                data[i++]=str2num(p,0,19);
            }
            /* decode ephemeris */
            if (sys==SYS_GLO&&i>=15) {
                if (!(mask&sys)) return 0;
                *type=1;
                return decode_geph(ver,sat,toc,data,geph);
            }
            else if (sys==SYS_SBS&&i>=15) {
                if (!(mask&sys)) return 0;
                *type=2;
                return decode_seph(ver,sat,toc,data,seph);
            }
            else if (i>=31) {
                if (!(mask&sys)) return 0;
                *type=0;
                return decode_eph(ver,sat,toc,data,eph);
            }
        }
    }
    return -1;
}
/* add ephemeris to navigation data ------------------------------------------*/
static int add_eph(nav_t *nav, const eph_t *eph)
{
    eph_t *nav_eph;
    
    if (nav->nmax<=nav->n) {
        nav->nmax+=1024;
        if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->nmax))) {
            free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
            return 0;
        }
        nav->eph=nav_eph;
    }
    nav->eph[nav->n++]=*eph;
    return 1;
}
static int add_geph(nav_t *nav, const geph_t *geph)
{
    geph_t *nav_geph;
    
    if (nav->ngmax<=nav->ng) {
        nav->ngmax+=1024;
        if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ngmax))) {
            free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
            return 0;
        }
        nav->geph=nav_geph;
    }
    nav->geph[nav->ng++]=*geph;
    return 1;
}
static int add_seph(nav_t *nav, const seph_t *seph)
{
    seph_t *nav_seph;
    
    if (nav->nsmax<=nav->ns) {
        nav->nsmax+=1024;
        if (!(nav_seph=(seph_t *)realloc(nav->seph,sizeof(seph_t)*nav->nsmax))) {
            free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
            return 0;
        }
        nav->seph=nav_seph;
    }
    nav->seph[nav->ns++]=*seph;
    return 1;
}
/* read rinex nav/gnav/geo nav -----------------------------------------------*/
extern int readrnxnav(FILE *fp, const char *opt, double ver, int sys,
                      nav_t *nav)
{
    eph_t eph;
    geph_t geph;
    seph_t seph;
    int stat,type;
    
    
    if (!nav) return 0;
    
    /* read rinex navigation data body */
    while ((stat=readrnxnavb(fp,opt,ver,sys,&type,&eph,&geph,&seph))>=0) {
        
        /* add ephemeris to navigation data */
        if (stat) {
            switch (type) {
                case 1 : stat=add_geph(nav,&geph); break;
                case 2 : stat=add_seph(nav,&seph); break;
                default: stat=add_eph (nav,&eph ); break;
            }
            if (!stat) return 0;
        }
    }
    return nav->n>0||nav->ng>0||nav->ns>0;
}
/* read rinex clock ----------------------------------------------------------*/
static int readrnxclk(FILE *fp, const char *opt, double ver, int index, nav_t *nav)
{
    pclk_t *nav_pclk;
    gtime_t time;
    double data[2];
    int i, j, sat, mask, off;
    char buff[MAXRNXLEN],satid[8]="";
   
    
    if (!nav) return 0;
    
    /* set system mask */
    mask=set_sysmask(opt);
    off = ver >= 3.04 ? 5 : 0; /* format change for ver>=3.04 */
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (str2time(buff, 8 + off, 26, &time)) {
            continue;
        }
        strncpy(satid,buff+3,4);
        
        /* only read AS (satellite clock) record */
        if (strncmp(buff,"AS",2)||!(sat=satid2no(satid))) continue;
        
        if (!(satsys(sat,NULL)&mask)) continue;
        
        for (i = 0, j = 40 + off; i < 2; i++, j += 20) {
            data[i] = str2num(buff, j, 19);  
        }

        if (nav->nc >= nav->ncmax) {
            nav->ncmax += 1024;
            if (!(nav_pclk = (pclk_t*)realloc(nav->pclk, sizeof(pclk_t) * (nav->ncmax)))) {
                free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
                return -1;
            }
            nav->pclk = nav_pclk;
        }

        if (nav->nc<=0||fabs(timediff(time,nav->pclk[nav->nc-1].time))>1E-9) {
            nav->nc++;
            nav->pclk[nav->nc-1].time =time;
            nav->pclk[nav->nc-1].index=index;
            for (i=0;i<MAXSAT;i++) {
                nav->pclk[nav->nc-1].clk[i][0]=0.0;
                nav->pclk[nav->nc-1].std[i][0]=0.0f;
            }
        }
        nav->pclk[nav->nc-1].clk[sat-1][0]=data[0];
        nav->pclk[nav->nc-1].std[sat-1][0]=(float)data[1];
    }
    return nav->nc>0;
}
/* read rinex file -----------------------------------------------------------*/
static int readrnxfp(FILE *fp, gtime_t ts, gtime_t te, double tint,
                     const char *opt, int flag, int index, char *type,
                     obs_t *obs, nav_t *nav, sta_t *sta)
{
    double ver;
    int sys,tsys=TSYS_GPS;
    char tobs[NUMSYS][MAXOBSTYPE][4]={{""}};
    
    
    /* read rinex header */
    if (!readrnxh(fp,&ver,type,&sys,&tsys,tobs,nav,sta)) return 0;
    
    /* flag=0:except for clock,1:clock */
    if ((!flag&&*type=='C')||(flag&&*type!='C')) return 0;
    
    /* read rinex body */
    switch (*type) {
        case 'O': return readrnxobs(fp,ts,te,tint,opt,index,ver,&tsys,tobs,obs,
                                    sta,NULL);
        case 'N': return readrnxnav(fp,opt,ver,sys    ,nav);
        case 'G': return readrnxnav(fp,opt,ver,SYS_GLO,nav);
        case 'H': return readrnxnav(fp,opt,ver,SYS_SBS,nav);
        case 'J': return readrnxnav(fp,opt,ver,SYS_QZS,nav); /* extension */
        case 'L': return readrnxnav(fp,opt,ver,SYS_GAL,nav); /* extension */
        case 'C': return readrnxclk(fp,opt,ver,index,nav);
    }
    return 0;
}
/* uncompress and read rinex file --------------------------------------------*/
static int readrnxfile(const char *file, gtime_t ts, gtime_t te, double tint,
                       const char *opt, int flag, int index, char *type,
                       obs_t *obs, nav_t *nav, sta_t *sta)
{
    FILE *fp;
    int cstat,stat;
    char tmpfile[1024];
    
    if (sta) init_sta(sta);
    
    /* uncompress file */
    if ((cstat=rcb_uncompress(file,tmpfile))<0) {
        return 0;
    }
    if (!(fp=fopen(cstat?tmpfile:file,"r"))) {
        return 0;
    }
    /* read rinex file */
    stat=readrnxfp(fp,ts,te,tint,opt,flag,index,type,obs,nav,sta);
    
    fclose(fp);
    
    /* delete temporary file */
    if (cstat) remove(tmpfile);
    
    return stat;
}

/* compare precise clock -----------------------------------------------------*/
static int cmppclk(const void *p1, const void *p2)
{
    pclk_t *q1=(pclk_t *)p1,*q2=(pclk_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise clock -----------------------------------------------------*/
static void combpclk(nav_t *nav)
{
    pclk_t *nav_pclk;
    int i,j,k;
    
    if (nav->nc<=0) return;
    
    qsort(nav->pclk,nav->nc,sizeof(pclk_t),cmppclk);
    
    for (i=0,j=1;j<nav->nc;j++) {
        if (fabs(timediff(nav->pclk[i].time,nav->pclk[j].time))<1E-9) {
            for (k=0;k<MAXSAT;k++) {
                if (nav->pclk[j].clk[k][0]==0.0) continue;
                nav->pclk[i].clk[k][0]=nav->pclk[j].clk[k][0];
                if (nav->pclk[i].clk[k][0] == 0) {
                    int a = 0;
                }
                nav->pclk[i].std[k][0]=nav->pclk[j].std[k][0];
            }
        }
        else if (++i<j) nav->pclk[i]=nav->pclk[j];
    }
    nav->nc=i+1;
    
    if (!(nav_pclk=(pclk_t *)realloc(nav->pclk,sizeof(pclk_t)*nav->nc))) {
        free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
        return;
    }
    nav->pclk=nav_pclk;
    nav->ncmax=nav->nc;

}
