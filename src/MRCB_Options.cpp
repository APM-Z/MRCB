#include "MRCB.h"

/* system options buffer */
static prcopt_t prcopt_;
static solopt_t solopt_;
static filopt_t filopt_;
static int antpostype_[2];
static double elmask_;
static double antpos_[2][3];
static char exsats_[1024];
static char snrmask_[NFREQ][1024];
extern opt_t sysoptz[];

opt_t sysopts[]={
    {"OBS_PATH",          2,  (void*)&filopt_.obs,          ""     },
    {"NAV_PATH",          2,  (void*)&filopt_.nav,          ""     },

    {"OUT_PRO",         2,  (void*)&filopt_.outFile_PRO,    ""     },
    {"OUT_PATH",        2,  (void*)&filopt_.outFile_OUT,     ""     },
    {"OUT_LOG",         2,  (void*)&filopt_.outFile_LOG,     ""     },

    {"Ephopt",          0,  (void*)&prcopt_.sateph,     "" },
    {"FiltMode",        0,  (void*)&prcopt_.soltype,     "" },
    {"Frequency",       0,  (void*)&prcopt_.nf,         "" },
    {"Elmask",          1,  (void*)&elmask_,            "" },
    {"Exclsats",        2,  (void*)exsats_,             ""},
    {"System",          0,  (void*)&prcopt_.navsys,     "" },
    {"UseBDS2",         0,  (void*)&prcopt_.bds2,       ""},
    {"DOY",             2,  (void*)&prcopt_.DOY,       ""},



    {"Continue",        0,  (void*)&prcopt_.CONTINUE,   ""        },


    {"OBS_DIR",       2,  (void*)&filopt_.ObsDir,          ""     },
    {"NAV_DIR",       2,  (void*)&filopt_.NavDir,          ""     },
    {"StaName",       2,  (void*)&prcopt_.staname,         ""     },
   




    {"",0,NULL,""} /* terminator */
};

/* discard space characters at tail */
static void chop(char *str)
{
    char *p;
    if ((p=strchr(str,'#'))) *p='\0'; /* comment */
    for (p=str+strlen(str)-1;p>=str&&!isgraph((int)*p);p--) *p='\0';
}
/* enum to string */
static int enum2str(char *s,const char *comment,int val)
{
    char str[32];
    const char* p, * q; 
    int n;

    n=sprintf(str,"%d:",val);
    if (!(p=strstr(comment,str))) {
        return sprintf(s,"%d",val);
    }
    if (!(q=strchr(p+n,','))&&!(q=strchr(p+n,')'))) {
        strcpy(s,p+n);
        return (int)strlen(p+n);
    }
    strncpy(s,p+n,q-p-n); s[q-p-n]='\0';
    return (int)(q-p-n);
}
/* string to enum */
static int str2enum(const char *str,const char *comment,int *val)
{
    const char *p;
    char s[32];

    for (p=comment;;p++) {
       if (!(p=strstr(p,str))) break;
       if (*(p-1)!=':') continue;
       for (p-=2;'0'<=*p&&*p<='9';p--) ;
       return sscanf(p+1,"%d",val)==1;
    }
    sprintf(s,"%30.30s:",str);
    if ((p=strstr(comment,s))) { /* number */
        return sscanf(p,"%d",val)==1;
    }
    return 0;
}
/* search option */
extern opt_t *searchopt(const char *name,const opt_t *opts)
{
    int i;

    for (i=0;*opts[i].name;i++) {
        if (strstr(opts[i].name,name)) return (opt_t *)(opts+i);
    }
    return NULL;
}
/* string to option value */
extern int str2opt(opt_t *opt,const char *str)
{
    switch (opt->format) {
        case 0: *(int    *)opt->var=atoi(str); break;
        case 1: *(double *)opt->var=atof(str); break;
        case 2: strcpy((char *)opt->var,str);  break;
        case 3: return str2enum(str,opt->comment,(int *)opt->var);
        default: return 0;
    }
    return 1;
}
/* option value to string */
extern int opt2str(const opt_t *opt,char *str)
{
    char *p=str;
    switch (opt->format) {
        case 0: p+=sprintf(p,"%d"   ,*(int   *)opt->var); break;
        case 1: p+=sprintf(p,"%.15g",*(double*)opt->var); break;
        case 2: p+=sprintf(p,"%s"   , (char  *)opt->var); break;
        case 3: p+=enum2str(p,opt->comment,*(int *)opt->var); break;
    }
    return (int)(p-str);
}
/* load options */
extern int loadopts(const char *file,opt_t *opts)
{
    FILE *fp;
    opt_t *opt;
    char buff[2048],*p;
    int n=0;

    if (!(fp=fopen(file,"r"))) return 0;

    while (fgets(buff,sizeof(buff),fp)) {
        n++;
        chop(buff);

        if (buff[0]=='\0'||!(p=strstr(buff,"="))) continue;

        *p++='\0';
        chop(buff);
        if (!(opt=searchopt(buff,opts))) continue;

        if (!str2opt(opt,p)) continue;
    }
    fclose(fp);

    return 1;
}
/* system options buffer to options */
static void buff2sysopts(void)
{
    double pos[3];
    char buff[1024],*p,*id;
    int i,j,sat,*ps;
    prcopt_.elmin=elmask_*D2R;
    for (i=0;i<1;i++) {
        ps=&prcopt_.rovpos;
        if (antpostype_[i]==0) { /* lat/lon/hgt */
            *ps=0;
            pos[0]=antpos_[i][0]*D2R;
            pos[1]=antpos_[i][1]*D2R;
            pos[2]=antpos_[i][2];
        }
        else if (antpostype_[i]==1) { /* xyz-ecef */
            *ps=0;
        }
        else *ps=antpostype_[i]-1;
    }
    /* excluded satellites */
    for (i=0;i<MAXSAT;i++) prcopt_.exsats[i]=0;
    if (exsats_[0]!='\0') {
        strcpy(buff,exsats_);
        for (p=strtok(buff," ");p;p=strtok(NULL," ")) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        for (j=0;j<9;j++) prcopt_.snrmask.mask[i][j]=0.0;
        strcpy(buff,snrmask_[i]);
        for (p=strtok(buff,","),j=0;p&&j<9;p=strtok(NULL,",")) {
            prcopt_.snrmask.mask[i][j++]=atof(p);
        }
    }
}
/* options to system options buffer */
static void sysopts2buff(void)
{
    double pos[3],*rr=NULL;
    char id[32],*p;
    int i,j,sat,*ps;

    elmask_=prcopt_.elmin*R2D;

    for (i=0;i<1;i++) {
        ps=&prcopt_.rovpos;
        if (*ps==0) {
            antpostype_[i]=0;
            ecef2pos(rr,pos);
            antpos_[i][0]=pos[0]*R2D;
            antpos_[i][1]=pos[1]*R2D;
            antpos_[i][2]=pos[2];
        }
        else antpostype_[i]=*ps+1;
    }
    /* excluded satellites */
    exsats_[0]='\0';
    for (sat=1,p=exsats_;sat<=MAXSAT&&p-exsats_<(int)sizeof(exsats_)-32;sat++) {
        if (prcopt_.exsats[sat-1]) {
            satno2id(sat,id);
            p+=sprintf(p,"%s%s%s",p==exsats_?"":" ",prcopt_.exsats[sat-1]==2?"+":"",id);
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        snrmask_[i][0]='\0';
        p=snrmask_[i];
        for (j=0;j<9;j++) {
            p+=sprintf(p,"%s%.0f",j>0?",":"",prcopt_.snrmask.mask[i][j]);
        }
    }
}
/* reset system options to default */
extern void resetsysopts(void)
{
    int i,j;
    prcopt_=prcopt_default;
    solopt_=solopt_default;
    filopt_.satantp[0]='\0';
    filopt_.rcvantp[0]='\0';
    filopt_.stapos [0]='\0';
    filopt_.solstat[0]='\0';
    for (i=0;i<2;i++) antpostype_[i]=0;
    elmask_=15.0;
    for (i=0;i<2;i++) for (j=0;j<3;j++) {
        antpos_[i][j]=0.0;
    }
    exsats_[0] ='\0';
}
/* get system options */
extern void getsysopts(prcopt_t *popt,solopt_t *sopt,filopt_t *fopt)
{
    buff2sysopts();
    if (popt) *popt=prcopt_;
    if (sopt) *sopt=solopt_;
    if (fopt) *fopt=filopt_;
}
/* set system options */
extern void setsysopts(const prcopt_t *prcopt,const solopt_t *solopt,const filopt_t *filopt)
{
    resetsysopts();
    if (prcopt) prcopt_=*prcopt;
    if (solopt) solopt_=*solopt;
    if (filopt) filopt_=*filopt;
    sysopts2buff();
}
