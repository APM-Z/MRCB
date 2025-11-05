#include "MRCB.h"
/* solution option to field separator */
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* output processing options */
static int outprcopts(unsigned char* buff, const prcopt_t* opt)
{
    const int sys[] = { SYS_GPS,SYS_GLO,SYS_GAL,SYS_QZS,SYS_CMP,SYS_IRN,SYS_SBS,0 };
    const char* s3[] = { "forward","backward","combined" };
    const char* s4[] = { "off","broadcast","sbas","iono-free","estimation","" };
    const char* s6[] = { "broadcast","precise","" };
    const char* s7[] = { "GPS","GLONASS","GALILEO","QZSS","BDS","IRNSS","SBAS","" };
    int i;
    char* p = (char*)buff;
    p += sprintf(p, "%s Station Name: %s\n", COMMENTH, opt->staname);
    p += sprintf(p, "%s Ephemeris   : %s\n", COMMENTH, s6[opt->sateph]);
    p += sprintf(p, "%s Elev mask   : %.1f deg\n", COMMENTH, opt->elmin * R2D);
    p += sprintf(p, "%s Freq number : %d\n", COMMENTH, opt->nf);
    p += sprintf(p, "%s Sol  Type   : %s\n", COMMENTH, s3[opt->soltype]);
    p += sprintf(p, "%s Ionos opt   : %s\n", COMMENTH, s4[opt->ionoopt]);
    p += sprintf(p, "%s Navi sys    :", COMMENTH);
    for (i = 0; sys[i]; i++) {
        if (opt->navsys & sys[i]) p += sprintf(p, " %s", s7[i]);
    }
    p += sprintf(p, "\n");
    return (int)(p - (char*)buff);
}
/* output solution header */
static int outsolheads(unsigned char* buff, const solopt_t* opt, const prcopt_t* popt)
{
    const char* sep = opt2sep(opt);
    char* p = (char*)buff;
    p += sprintf(p, "#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    p += sprintf(p, "#%15s%s%20s%s%14s%s%14s%s%4s%s%6s", "Time(GPST)", sep, "x-ecef(m)", sep, "y-ecef(m)", sep, "z-ecef(m)", sep, "ns", sep, "GDOP");
    if (popt->bds2 == 0) {
        p += sprintf(p, "%s%s%10s%s%8s%s%8s%s", sep, sep, "G_P1(m)", sep, "G_P2(m)", sep, "G_P5(m)", sep);
        p += sprintf(p, "%s%9s%s%8s%s", sep, "R_R1(m)", sep, "R_R2(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s%8s%s%8s%s%8s%s", sep, "GAL_E1(m)", sep, "GAL_E5a(m)", sep, "GAL_E5b(m)", sep, "GAL_E6(m)", sep, "GAL_E5(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s%8s%s", sep, "J_P1(m)", sep, "J_P2(m)", sep, "J_P5(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s%8s%s%8s%s%8s", sep, "BDS3_B1I(m)", sep, "BDS3_B3I(m)", sep, "BDS3_B2b(m)", sep, "BDS3_B1C(m)", sep, "BDS3_B2a(m)");
    }
    else if (popt->bds2 == 1) {
        p += sprintf(p, "%s%s%10s%s%8s%s%8s%s", sep, sep, "G_P1(m)", sep, "G_P2(m)", sep, "G_P5(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s", sep, "R_R1(m)", sep, "R_R2(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s%8s%s%8s%s%8s%s", sep, "GAL_E1(m)", sep, "GAL_E5a(m)", sep, "GAL_E5b(m)", sep, "GAL_E6(m)", sep, "GAL_E5(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s%8s%s", sep, "J_P1(m)", sep, "J_P2(m)", sep, "J_P5(m)", sep);
        p += sprintf(p, "%s%10s%s%8s%s%8s", sep, "BDS2_B1I(m)", sep, "BDS2_B3I(m)", sep, "BDS2_B2b(m)");
    }
    p += sprintf(p, "\n");
    p += sprintf(p, "#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    return (int)(p - (char*)buff);
}
/* output solution header */
static void outsolhead(FILE* fp, const solopt_t* opt, const prcopt_t* popt)
{
    unsigned char buff[MAXSOLMSG + 1];
    int n;

    if ((n = outsolheads(buff, opt, popt)) > 0) {
        fwrite(buff, n, 1, fp);
    }
}
/* output processing option */
static void outprcopt(FILE* fp, const prcopt_t* opt)
{
    unsigned char buff[MAXSOLMSG + 1];
    int n;
    if ((n = outprcopts(buff, opt)) > 0) {
        fwrite(buff, n, 1, fp);
    }
}

/* output header */
static void outheader(FILE* fp, const prcopt_t* popt, const solopt_t* sopt)
{
    const char* s1[] = { "GPST","UTC","JST" };
    outprcopt(fp, popt);
    fprintf(fp, "%s\n", COMMENTH);
    outsolhead(fp, sopt, popt);
}
/* open output file for append */
static FILE* openfile(const char* outfile)
{
    return fopen(outfile, "w");
}
///* creat out file */
int creatFile(FILE** fp, const filopt_t* fopt, const prcopt_t* popt, const solopt_t* sopt, char* filename, char* profilename, char* logfilename) {
    char outfile[3][1024] = { "","" };
    char buff[1024], * p = buff;
#ifdef _WIN32
    snprintf(outfile[0], 1024, "%s%s.%s", fopt->outFile_PRO, profilename, "pro");
    snprintf(outfile[1], 1024, "%s%s.%s", fopt->outFile_OUT, filename, "out");
    snprintf(outfile[2], 1024, "%s%s.%s", fopt->outFile_LOG, logfilename, "log");
#else
    snprintf(outfile[0], 1024, "%s%s.%s", fopt->outFile_PRO, profilename, "pro");
    snprintf(outfile[1], 1024, "%s%s.%s", fopt->outFile_OUT, filename, "out");
    snprintf(outfile[2], 1024, "%s%s.%s", fopt->outFile_LOG, logfilename, "log");
#endif
    if ((fp[0] = openfile(outfile[0])) && (fp[1] = openfile(outfile[1])) && (fp[2] = openfile(outfile[2]))) {
        outheader(fp[1], popt, sopt);
        p += sprintf(p, "%s", "sat    sat-X(m)        sat-Y(m)         sat-Z(m)        sat-clock(m)");
        p += sprintf(p, "%s", "         ion(m)");
        p += sprintf(p, "%s", "       amb-f1(cycle)  amb-f2(cycle)    amb-f3(cycle)  amb-f4(cycle)   amb-f5(cycle)");
        p += sprintf(p, "%s", "     elevation      azimuth       resL-f1(m)       resL-f2(m)       resL-f3(m)       resL-f4(m)       resL-f5(m)");
        p += sprintf(p, "%s", "    resP-f1(m)      resP-f2(m)        resP-f3(m)      resP-f4(m)      resP-f5(m)");
        p += sprintf(p, "%s\n", "       SCB-f3(m)       SCB-f4(m)       SCB-f5(m)");
        fwrite(buff, (int)(p - buff), 1, fp[0]);
        return 1;
    }
    else {
        printf("File creat error,please check the out path!\n");
        return 0;
    }
}

