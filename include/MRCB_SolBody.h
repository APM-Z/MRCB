#ifndef MRCB_SolBody_H
#define MRCB_SolBody_H

#include "MRCB.h"

void OutPut(FILE** fp, rcb_t* rcb, rscs_t* rscs, const obsd_t* obs, int n, int stat, sta_t* sta, const solopt_t* sopt);
void OutLogElv(FILE** fp, double maskElv, char id[4], double rsElv);
void OutLogFrq(FILE** fp, char id[4], int nf, int f);
#endif
