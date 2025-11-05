#include "MRCB_SolBody.h"

/* solution option to field separator */
static const char* opt2sep(const solopt_t* opt)
{
	if (!*opt->sep) return " ";
	else if (!strcmp(opt->sep, "\\t")) return "\t";
	return opt->sep;
}

/* output solution as the form of x/y/z-ecef */
static int outecef_float(uint8_t* buff, const char* s, const sol_t* sol, const prcopt_t* opt, sta_t* sta, const solopt_t* sopt)
{
	const char* sep = opt2sep(sopt);
	char* p = (char*)buff;
	double x = 0.0, y = 0.0, z = 0.0;
	x = sol->rr[0];
	y = sol->rr[1];
	z = sol->rr[2];
	if (opt->bds2 == 0) {
		p += sprintf(p, "%s%s%14.4f%s%14.4f%s%14.4f%s%4d%s%8.4f%s%9.4f%s%8.4f%s%8.4f%s%10.4f%s%8.4f%s%11.4f%s%9.4f%s%10.4f%s%9.4f%s%9.4f%s%12.4f%s%8.4f%s%8.4f%s%9.4f%s%11.4f%s%11.4f%s%12.4f%s%10.4f",
			s, sep, x, sep, y, sep, z, sep, sol->ns, sep, sol->dop[0], sep,
			sol->rcb[0][0], sep, sol->rcb[0][1], sep, sol->rcb[0][2], sep,//gps
			sol->rcb[1][0], sep, sol->rcb[1][1], sep,//glo
			sol->rcb[2][0], sep, sol->rcb[2][1], sep, sol->rcb[2][2], sep, sol->rcb[2][3], sep, sol->rcb[2][4], sep,//gal
			sol->rcb[3][0], sep, sol->rcb[3][1], sep, sol->rcb[3][2], sep, //qzss
			sol->rcb[4][0], sep, sol->rcb[4][1], sep, sol->rcb[4][2], sep, sol->rcb[4][3], sep, sol->rcb[4][4]);//bds
	}
	else if (opt->bds2 == 1) {
		p += sprintf(p, "%s%s%14.4f%s%14.4f%s%14.4f%s%4d%s%8.4f%s%9.4f%s%8.4f%s%8.4f%s%10.4f%s%8.4f%s%11.4f%s%9.4f%s%10.4f%s%9.4f%s%9.4f%s%12.4f%s%8.4f%s%8.4f%s%9.4f%s%11.4f%s%11.4f",
			s, sep, x, sep, y, sep, z, sep, sol->ns, sep, sol->dop[0], sep,
			sol->rcb[0][0], sep, sol->rcb[0][1], sep, sol->rcb[0][2], sep,
			sol->rcb[1][0], sep, sol->rcb[1][1], sep,//glo
			sol->rcb[2][0], sep, sol->rcb[2][1], sep, sol->rcb[2][2], sep, sol->rcb[2][3], sep, sol->rcb[2][4], sep,//gal
			sol->rcb[3][0], sep, sol->rcb[3][1], sep, sol->rcb[3][2], sep, //qzss
			sol->rcb[4][0], sep, sol->rcb[4][1], sep, sol->rcb[4][2]);
	}
	p += sprintf(p, "\n");
	return (int)((p - (char*)buff));
}

void OutPut(FILE** fp, rcb_t* rcb, rscs_t* rscs, const obsd_t* obs, int n, int stat, sta_t* sta, const solopt_t* sopt)
{
	double ep[6] = { 0.0 };
	const prcopt_t* opt = &rcb->opt;
	int i, j, sat = 0, sys = 0;
	rcb->sol.time_filter = rcb->sol.time;
	/* valid satellites */
	rcb->sol.ns = 0;
	for (int i = 0; i < NSYS; i++) {
		rcb->sol.sys_ns[i] = 0;
	}
	for (i = 0; i < n && i < MAXOBS; ++i) {
		if ((sat = obs[i].sat) <= 0) continue;
		if (rcb->cvsat[sat - 1] != 0) {
			sys = satsys(sat, NULL);
			++rcb->sol.ns;
			++rcb->sol.sys_ns[sys2num(sys) - 1];
		}
	}
	/*initialize*/
	gtime_t time;
	int timeu = 2;
	char s[64], id[4];
	uint8_t buff[8192];
	uint8_t* p = buff;
	const char* sep = opt2sep(sopt);
	time = rcb->sol.time;
	time2str(time, s, timeu);
	/*time*/
	p = buff;
	p += sprintf((char*)p, "%s", s);
	p += sprintf((char*)p, "\n");
	fwrite(buff, (int)(p - buff), 1, fp[0]);
	/*output float position data*/
	p = buff;
	int bytes_written = outecef_float(p, s, &rcb->sol, opt, sta, sopt);
	p += bytes_written;
	fwrite(buff, (int)(p - buff), 1, fp[1]);
	/*satellie and receiver info*/
	for (i = 0; i < n && i < MAXOBS; ++i) {
		/*check*/
		sat = obs[i].sat;
		if (!rcb->cvsat[sat - 1]) continue;
		sys = sys2num(satsys(sat, NULL)) - 1;
		/*out satellite*/
		/*satellite position and satellite clock*/
		p = buff;
		satno2id(sat, id);
		p += sprintf((char*)p, "%s%s%15.4f%s%15.4f%s%15.4f%s%15.4f%s", id, sep, rscs[i].rs[0], sep, rscs[i].rs[1], sep, rscs[i].rs[2], sep, CLIGHT * rscs[i].dts[0], sep);
		fwrite(buff, (int)(p - buff), 1, fp[0]);
		/*ionosphere*/
		p = buff;
		p += sprintf((char*)p, "%15.4f%s", rcb->sol.ion[sat - 1], sep);
		fwrite(buff, (int)(p - buff), 1, fp[0]);
		/*ambiguity*/
		p = buff;
		for (j = 0; j < NFREQ; ++j) p += sprintf((char*)p, "%15.4f%s", rcb->sol.amb[sat - 1][j], sep);
		fwrite(buff, (int)(p - buff), 1, fp[0]);
		/*elevation*/
		p = buff;
		p += sprintf((char*)p, "%15.4f%s%15.4f%s", rscs[i].azel[1] * R2D, sep, rscs[i].azel[0] * R2D, sep);
		fwrite(buff, (int)(p - buff), 1, fp[0]);
		/*residuals*/
		p = buff;
		for (j = 0; j < NFREQ; ++j) {
			p += sprintf((char*)p, "%15.4f%s", rcb->sol.resc[sat - 1][j], sep);
		}
		p += sprintf((char*)p, "%s", sep);
		for (j = 0; j < NFREQ; ++j) {
			p += sprintf((char*)p, "%15.4f%s", rcb->sol.resp[sat - 1][j], sep);
		}
		fwrite(buff, (int)(p - buff), 1, fp[0]);
		/*SCB*/
		p = buff;
		for (j = 2; j < NFREQ; ++j) {
			p += sprintf((char*)p, "%15.4f%s", rcb->sol.scb[sat - 1][j], sep);
		}
		p += sprintf((char*)p, "\n");
		fwrite(buff, (int)(p - buff), 1, fp[0]);
	}
}

void OutLogElv(FILE** fp, double maskElv, char id[4], double rsElv) {
	char buff[8192], * p = buff;
	p = buff;
	p += sprintf(p, "Deleter Sat for lower %.1f degree:%s(Ele:%.1f)\n", maskElv, id, rsElv);
	fwrite(buff, (int)(p - buff), 1, fp[2]);
}

void OutLogFrq(FILE** fp, char id[4], int nf, int f) {
	char buff[8192], * p = buff;
	p = buff;
	p += sprintf(p, "Deleter Sat for lack observations, need freq %d:%s(Frq:%d)\n", nf, id, f);
	fwrite(buff, (int)(p - buff), 1, fp[2]);
}
